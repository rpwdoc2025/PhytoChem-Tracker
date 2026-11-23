"""Wikidata SPARQL: taxon resolution, LOTUS occurrences, cross-references.

The LOTUS initiative curates >750,000 referenced structure-organism pairs and
publishes them on Wikidata as 'found in taxon' (P703) statements, each with a
literature reference. Wikidata also links compounds to KEGG, ChEBI, KNApSAcK
and PubChem identifiers, which makes it the natural hub for entity resolution
(names <-> InChIKey <-> taxon).
"""

from __future__ import annotations

from typing import Any

from .. import config
from ..http import Client, HTTPError

SOURCE = "Wikidata SPARQL (LOTUS)"

PREFIXES = """PREFIX wd: <http://www.wikidata.org/entity/>
PREFIX wdt: <http://www.wikidata.org/prop/direct/>
PREFIX p: <http://www.wikidata.org/prop/>
PREFIX ps: <http://www.wikidata.org/prop/statement/>
PREFIX pr: <http://www.wikidata.org/prop/reference/>
PREFIX prov: <http://www.w3.org/ns/prov#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
"""


def _esc(s: str) -> str:
    return s.replace("\\", "\\\\").replace('"', '\\"')


def sparql(client: Client, query: str) -> list[dict[str, str]]:
    """Run a query on Wikidata; fall back to the QLever mirror on failure."""
    q = PREFIXES + query
    headers = {"Accept": "application/sparql-results+json"}
    last_exc: Exception | None = None
    for endpoint in (config.WIKIDATA_SPARQL, config.WIKIDATA_SPARQL_FALLBACK):
        if not endpoint:
            continue
        try:
            data = client.get_json(endpoint, source=SOURCE, params={"query": q, "format": "json"},
                                   headers=headers, timeout=60)
            rows = []
            for b in data.get("results", {}).get("bindings", []):
                rows.append({k: v.get("value", "") for k, v in b.items()})
            return rows
        except (HTTPError, ValueError) as exc:
            last_exc = exc
    raise HTTPError(f"Wikidata query failed on all endpoints: {last_exc}")


def qid(uri: str) -> str:
    return uri.rsplit("/", 1)[-1] if uri else ""


def resolve_taxon(client: Client, name: str) -> dict[str, Any] | None:
    """Find a plant taxon by scientific name, or by English common name."""
    name = " ".join(name.split())
    variants = {name, name.replace(" x ", " × ")}
    sci_values = " ".join(f'"{_esc(v)}"' for v in variants)
    query = f"""SELECT ?t ?name ?ncbi WHERE {{
  VALUES ?name {{ {sci_values} }}
  ?t wdt:P225 ?name .
  OPTIONAL {{ ?t wdt:P685 ?ncbi }}
}} LIMIT 5"""
    rows = sparql(client, query)
    if not rows:
        commons = {name, name.lower(), name.title(), name.capitalize()}
        values = " ".join(f'"{_esc(v)}"@en' for v in commons)
        query = f"""SELECT ?t ?name ?ncbi WHERE {{
  VALUES ?common {{ {values} }}
  ?t wdt:P1843 ?common ; wdt:P225 ?name .
  OPTIONAL {{ ?t wdt:P685 ?ncbi }}
}} LIMIT 5"""
        rows = sparql(client, query)
    if not rows:
        return None
    # Prefer the record that has an NCBI taxonomy ID (usually the accepted name)
    rows.sort(key=lambda r: (not r.get("ncbi"), r.get("name", "")))
    best = rows[0]
    t = qid(best["t"])
    commons = sparql(client, f"""SELECT ?c WHERE {{ wd:{t} wdt:P1843 ?c . FILTER(LANG(?c) = "en") }} LIMIT 10""")
    return {
        "qid": t, "name": best.get("name", name), "ncbi_taxid": best.get("ncbi", ""),
        "common_names": sorted({r["c"] for r in commons}),
    }


def compounds_in_taxon(client: Client, taxon_qid: str, include_infraspecific: bool = True,
                       limit: int = 5000) -> list[dict[str, Any]]:
    """All compounds with a 'found in taxon' statement for the taxon (and its varieties)."""
    taxon_clause = (f"?t wdt:P171* wd:{taxon_qid} ." if include_infraspecific
                    else f"VALUES ?t {{ wd:{taxon_qid} }}")
    query = f"""SELECT ?c ?ik (SAMPLE(?lbl) AS ?label) (SAMPLE(?smi1) AS ?smiles_a) (SAMPLE(?smi2) AS ?smiles_b)
       (SAMPLE(?pc) AS ?cid) (COUNT(DISTINCT ?ref) AS ?refs) WHERE {{
  {taxon_clause}
  ?c p:P703 ?st . ?st ps:P703 ?t .
  ?c wdt:P235 ?ik .
  OPTIONAL {{ ?st prov:wasDerivedFrom ?d . ?d pr:P248 ?ref . }}
  OPTIONAL {{ ?c wdt:P233 ?smi1 }}
  OPTIONAL {{ ?c wdt:P2017 ?smi2 }}
  OPTIONAL {{ ?c wdt:P662 ?pc }}
  OPTIONAL {{ ?c rdfs:label ?lbl . FILTER(LANG(?lbl) = "en") }}
}} GROUP BY ?c ?ik LIMIT {int(limit)}"""
    out = []
    for r in sparql(client, query):
        out.append({
            "qid": qid(r.get("c", "")),
            "name": r.get("label") or r.get("ik", ""),
            "has_label": bool(r.get("label")),
            "inchikey": r.get("ik", ""),
            "smiles": r.get("smiles_b") or r.get("smiles_a") or "",
            "cid": int(r["cid"]) if r.get("cid", "").isdigit() else None,
            "references": int(float(r.get("refs") or 0)),
        })
    return out


XREF_PROPS = {"kegg": "P665", "chebi": "P683", "knapsack": "P2064", "pubchem_cid": "P662",
              "chembl": "P592", "hmdb": "P2057"}


def compound_xrefs(client: Client, inchikeys: list[str]) -> dict[str, dict[str, str]]:
    """Cross-references for a batch of InChIKeys -> {inchikey: {kegg: C00389, ...}}."""
    keys = [k for k in inchikeys if k]
    if not keys:
        return {}
    values = " ".join(f'"{_esc(k)}"' for k in keys)
    selects = " ".join(f"(SAMPLE(?{k}) AS ?{k}_)" for k in XREF_PROPS)
    optionals = "\n".join(f"  OPTIONAL {{ ?c wdt:{p} ?{k} }}" for k, p in XREF_PROPS.items())
    query = f"""SELECT ?ik (SAMPLE(?c) AS ?item) {selects} WHERE {{
  VALUES ?ik {{ {values} }}
  ?c wdt:P235 ?ik .
{optionals}
}} GROUP BY ?ik"""
    out: dict[str, dict[str, str]] = {}
    for r in sparql(client, query):
        d = {k: r.get(f"{k}_", "") for k in XREF_PROPS}
        d["wikidata"] = qid(r.get("item", ""))
        out[r["ik"]] = {k: v for k, v in d.items() if v}
    return out


def taxa_for_compound(client: Client, inchikey: str, limit: int = 50) -> list[dict[str, Any]]:
    """Reverse lookup: organisms in which a compound has been reported (compound mode, tile 1)."""
    query = f"""SELECT ?t ?name (COUNT(DISTINCT ?ref) AS ?refs) WHERE {{
  ?c wdt:P235 "{_esc(inchikey)}" ; p:P703 ?st .
  ?st ps:P703 ?t . ?t wdt:P225 ?name .
  OPTIONAL {{ ?st prov:wasDerivedFrom ?d . ?d pr:P248 ?ref . }}
}} GROUP BY ?t ?name ORDER BY DESC(?refs) ?name LIMIT {int(limit)}"""
    return [{"taxon": r.get("name", ""), "wikidata": qid(r.get("t", "")),
             "references": int(float(r.get("refs") or 0))} for r in sparql(client, query)]
