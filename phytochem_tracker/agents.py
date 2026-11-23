"""The specialised agents behind the seven tiles.

Each agent is a small, deterministic function that calls one or more
connectors, records provenance, and never raises: failures are stored in
``details.errors`` so one slow or broken database cannot break the report.
"""

from __future__ import annotations

from typing import Any, Callable

from . import chem, config
from .connectors import kegg, links, mona, pubchem, pubmed, uniprot, wikidata
from .config import Settings
from .http import Client
from .models import (
    Compound, CompoundDetails, GeneMap, PathwayRef, PlantTaxon, RankedCompound, SpectrumRef, UseRecord,
)

Progress = Callable[[str, str], None]


def _safe(details: CompoundDetails, tile: str, fn: Callable[[], Any]) -> Any:
    try:
        return fn()
    except Exception as exc:  # noqa: BLE001 - agents must not crash the report
        details.errors[tile] = f"{type(exc).__name__}: {exc}"
        return None


# ==========================================================================
# Tile 1 - Literature agent: LOTUS candidates ranked by PubMed paper counts
# ==========================================================================
class LiteratureAgent:
    name = "Literature agent"

    def __init__(self, client: Client, settings: Settings, progress: Progress):
        self.client, self.settings, self.progress = client, settings, progress
        self.n_ranked = 0

    def candidates(self, taxon: PlantTaxon) -> list[dict[str, Any]]:
        self.progress(self.name, f"Querying LOTUS (Wikidata) for compounds reported in {taxon.name}…")
        rows = wikidata.compounds_in_taxon(self.client, taxon.qid, self.settings.include_infraspecific)
        self.progress(self.name, f"LOTUS returned {len(rows)} structure–organism pairs.")
        return rows

    def rank(self, taxon: PlantTaxon, rows: list[dict[str, Any]]) -> tuple[list[RankedCompound], dict[int, dict]]:
        s = self.settings
        # Pre-filter: compounds documented by more LOTUS references first. Entries with neither an
        # English label nor a PubChem CID cannot be searched in PubMed by name, so they are skipped.
        usable = [r for r in rows if r["has_label"] or r["cid"]]
        if len(usable) < len(rows):
            self.progress(self.name, f"Skipped {len(rows) - len(usable)} LOTUS entries with no name or PubChem CID.")
        pool = sorted(usable, key=lambda r: (-r["references"], not r["has_label"], r["name"]))[: s.n_candidates]

        # Resolve PubChem CIDs (Wikidata usually has them; otherwise look up by InChIKey).
        for r in pool:
            if not r["cid"] and r["inchikey"]:
                cids = pubchem.cids_for(self.client, r["inchikey"], "inchikey")
                r["cid"] = cids[0] if cids else None
        cids = [r["cid"] for r in pool if r["cid"]]
        props = pubchem.properties(self.client, cids) if cids else {}

        self.progress(self.name, f"Counting PubMed papers for {len(pool)} candidates × '{taxon.name}'…")
        ranked: list[RankedCompound] = []
        for r in pool:
            rec = props.get(r["cid"], {}) if r["cid"] else {}
            comp = pubchem.compound_from_props(rec, fallback_name=r["name"]) if rec else Compound(
                name=r["name"], inchikey=r["inchikey"], smiles=r["smiles"])
            comp.qid = r["qid"]
            if not comp.smiles:
                comp.smiles = r["smiles"]
            names = [r["name"]] if r["has_label"] else []
            if comp.name and comp.name.lower() not in {n.lower() for n in names}:
                names.append(comp.name)
            if s.synonyms_per_compound and comp.cid:
                syns = pubchem.synonyms(self.client, comp.cid)
                extra = pubchem.literature_synonyms([x for x in syns if x.lower() not in {n.lower() for n in names}],
                                                    s.synonyms_per_compound)
                names += extra
            names = pubchem.literature_synonyms(names, 1 + s.synonyms_per_compound)
            if not names:  # only registry-style names available: not searchable in PubMed
                continue
            comp.synonyms = names
            query = pubmed.build_query(names, taxon.name)
            try:
                n = pubmed.count(self.client, query, s.ncbi_api_key, s.ncbi_email)
            except Exception as exc:  # noqa: BLE001
                self.progress(self.name, f"PubMed count failed for {comp.name}: {exc}")
                n = -1
            ranked.append(RankedCompound(rank=0, compound=comp, paper_count=n, lotus_references=r["references"],
                                         pubmed_query=query, sources=["LOTUS/Wikidata", "PubMed"]))
        self.n_ranked = len(ranked)
        ranked.sort(key=lambda x: (-x.paper_count, -x.lotus_references, x.compound.name.lower()))
        top = ranked[: s.n_top]
        for i, rc in enumerate(top, 1):
            rc.rank = i
        return top, props


# ==========================================================================
# Tile 2 - Physico-chemical agent
# ==========================================================================
def physchem_agent(client: Client, d: CompoundDetails, pubchem_rec: dict[str, Any] | None) -> None:
    c = d.compound
    if pubchem_rec:
        d.computed["PubChem (computed)"] = pubchem.computed_table(pubchem_rec)
    rd = chem.descriptors(c.smiles, c.inchi)
    if rd:
        d.computed["RDKit (local)"] = rd
    if c.cid:
        rec = pubchem.pug_view(client, c.cid, "Experimental Properties")
        d.experimental = pubchem.experimental_properties(rec)
        d.solid_state = pubchem.solid_state(d.experimental)
        d.add_prov("physchem", client.provenance(pubchem.VIEW_SOURCE,
                   f"pug_view/data/compound/{c.cid}?heading=Experimental Properties", n=len(d.experimental)))


# ==========================================================================
# Tile 3 - Spectral agent
# ==========================================================================
def spectral_agent(client: Client, d: CompoundDetails, max_mona: int = 5) -> None:
    c = d.compound
    d.links.update({f"spectra::{k}": v for k, v in links.spectra_links(c).items()})
    if c.cid:
        rec = pubchem.pug_view(client, c.cid, "Spectral Information")
        for g in pubchem.spectral_annotations(rec):
            d.spectra.append(SpectrumRef(
                type=g["type"], source=", ".join(g["sources"][:4]) or "PubChem",
                url=(g["urls"][0] if g["urls"] else links.spectra_links(c).get("PubChem spectral information", "")),
                description=f"{g['records']} record(s) under '{g['heading']}'" + (f" — {'; '.join(g['text'])}" if g["text"] else ""),
            ))
        d.add_prov("spectra", client.provenance(pubchem.VIEW_SOURCE,
                   f"pug_view/data/compound/{c.cid}?heading=Spectral Information", n=len(d.spectra)))
    if c.inchikey:
        try:
            recs = mona.spectra_for_inchikey(client, c.inchikey, size=20)
        except Exception as exc:  # noqa: BLE001
            d.errors["spectra_mona"] = str(exc)
            recs = []
        recs.sort(key=lambda r: (r["ms_level"] != "MS2", r["ionization"] != "positive"))
        for r in recs[:max_mona]:
            desc = " · ".join(x for x in [r["ms_level"], r["precursor_type"], r["ionization"], r["instrument"],
                                         f"CE {r['collision_energy']}" if r["collision_energy"] else ""] if x)
            level = {"MS2": "MS/MS", "MS3": "MS3", "MS1": "MS"}.get(r["ms_level"], r["ms_level"] or "MS")
            d.spectra.append(SpectrumRef(type=f"{level} (reference spectrum)", source="MoNA", url=r["url"], accession=r["id"],
                                         description=desc, peaks=r["peaks"]))
        d.add_prov("spectra", client.provenance(mona.SOURCE, f"spectra/search InChIKey={c.inchikey}", n=len(recs)))


# ==========================================================================
# Tile 4 - Use-case agent
# ==========================================================================
def uses_agent(client: Client, d: CompoundDetails) -> None:
    c = d.compound
    d.links.update({f"uses::{k}": v for k, v in links.uses_links(c).items()})
    if not c.cid:
        return
    records = [pubchem.pug_view(client, c.cid, h) for h in pubchem.USE_HEADINGS]
    for u in pubchem.uses_from_records(records):
        d.uses.append(UseRecord(category=u["category"], text=u["text"], source=u["source"], url=u["url"]))
    d.add_prov("uses", client.provenance(pubchem.VIEW_SOURCE,
               f"pug_view/data/compound/{c.cid}?heading=" + " | ".join(pubchem.USE_HEADINGS), n=len(d.uses)))


# ==========================================================================
# Tile 5 - Structure agent
# ==========================================================================
def structure_agent(client: Client, d: CompoundDetails) -> None:
    c = d.compound
    s: dict[str, Any] = {}
    if c.smiles:
        s["svg"] = chem.svg_2d(c.smiles)
        s["mol2d"] = chem.molblock_2d(c.smiles, c.name)
        s["sdf2d"] = chem.sdf_from_molblock(s["mol2d"], {"NAME": c.name, "INCHIKEY": c.inchikey,
                                                         "PUBCHEM_CID": c.cid, "SMILES": c.smiles})
    if c.cid:
        s["png_url"] = config.LINKS["pubchem_png"].format(cid=c.cid)
        sdf3d = pubchem.sdf(client, c.cid, three_d=True)
        if sdf3d.strip():
            s["sdf3d"], s["sdf3d_source"] = sdf3d, "PubChem 3D conformer"
    if "sdf3d" not in s and c.smiles:
        mb = chem.molblock_3d(c.smiles, c.name)
        if mb:
            s["sdf3d"], s["sdf3d_source"] = chem.sdf_from_molblock(mb, {"NAME": c.name}), "RDKit ETKDGv3 + MMFF (computed)"
    s["identifiers"] = {"Name": c.name, "PubChem CID": c.cid, "InChIKey": c.inchikey, "InChI": c.inchi,
                        "SMILES": c.smiles, "Formula": c.formula}
    d.structure = s
    d.add_prov("structure", client.provenance("RDKit + PubChem", f"CID {c.cid} / SMILES", note=s.get("sdf3d_source", "")))


# ==========================================================================
# Tile 6 - Pathway agent
# ==========================================================================
def pathway_agent(client: Client, d: CompoundDetails, taxon: PlantTaxon | None, max_pathway_enzymes: int = 3) -> None:
    c = d.compound
    d.links.update({f"pathways::{k}": v for k, v in links.pathway_links(c).items()})
    kid = c.xrefs.get("kegg") or kegg.find_compound(client, [c.name] + c.synonyms, c.formula)
    if not kid:
        d.add_prov("pathways", client.provenance(kegg.SOURCE, f"find/compound {c.formula} / {c.name}", n=0,
                                                 note="no KEGG compound entry found"))
        return
    c.xrefs["kegg"] = kid
    maps = kegg.pathways_for_compound(client, kid)
    names = kegg.pathway_names(client) if maps else {}
    org_ids: set[str] = set()
    org = taxon.kegg_org if taxon else ""
    if org and maps:
        org_ids = set(kegg.pathway_names(client, org))
    ecs = kegg.enzymes_for_compound(client, kid)
    d.enzymes_ec = ecs
    for i, m in enumerate(maps):
        org_pid = m.replace("map", org) if org and m.replace("map", org) in org_ids else ""
        pid = org_pid or m
        enz: list[str] = []
        if i < max_pathway_enzymes and m not in kegg.GLOBAL_MAPS:
            enz = sorted(set(kegg.enzymes_for_pathway(client, m)) & set(ecs))
        d.pathways.append(PathwayRef(
            name=names.get(m, m) + (f" ({taxon.name})" if org_pid and taxon else ""),
            db="KEGG", pathway_id=pid, url=config.LINKS["kegg_pathway"].format(pathway=pid, kegg=kid),
            enzymes=enz,
        ))
    d.add_prov("pathways", client.provenance(kegg.SOURCE, f"link/pathway/cpd:{kid}; link/enzyme/cpd:{kid}",
                                             n=len(maps), note=f"organism code: {org or 'none (reference maps)'}"))


# ==========================================================================
# Tile 7 - Gene agent: EC -> species-specific gene/protein identifiers
# ==========================================================================
def gene_agent(client: Client, d: CompoundDetails, taxon: PlantTaxon | None, settings: Settings) -> None:
    if not taxon:
        d.genes_status = "no species context"
        return
    ecs = list(d.enzymes_ec)
    pathway_ecs = [e for p in d.pathways for e in p.enzymes]
    ecs = sorted(set(ecs), key=lambda e: (e not in pathway_ecs, e))[: settings.max_ec_per_compound]
    if not ecs:
        d.genes_status = "no EC numbers linked to this compound in KEGG"
        return
    enames = kegg.enzyme_names(client, ecs)
    for ec in ecs:
        ename = enames.get(ec, "")
        if taxon.kegg_org:
            for g in kegg.genes_for_ec(client, taxon.kegg_org, ec)[:5]:
                d.genes.append(GeneMap(ec=ec, enzyme_name=ename, species=taxon.name, gene_id=g, db="KEGG GENES",
                                       url=f"https://www.kegg.jp/entry/{g}", evidence="KEGG orthology annotation"))
        for p in uniprot.proteins_for_ec(client, ec, taxon.ncbi_taxid, taxon.name, size=5):
            d.genes.append(GeneMap(
                ec=ec, enzyme_name=ename, species=p["organism"] or taxon.name,
                gene_id=p["gene_names"].split(" ")[0] if p["gene_names"] else p["accession"],
                db="UniProtKB/" + ("Swiss-Prot" if p["reviewed"] == "reviewed" else "TrEMBL"),
                url=config.LINKS["uniprot_entry"].format(acc=p["accession"]),
                protein=f"{p['accession']} — {p['protein']}"[:160],
                evidence=("Phytozome: " + p["phytozome"]) if p.get("phytozome") else "",
            ))
    d.genes_status = f"{len(d.genes)} gene/protein records for {len(ecs)} EC numbers"
    d.add_prov("genes", client.provenance(uniprot.SOURCE, f"ec:<EC> AND taxonomy_id:{taxon.ncbi_taxid or taxon.name}",
                                          n=len(d.genes), note=f"KEGG organism: {taxon.kegg_org or 'not in KEGG'}"))
