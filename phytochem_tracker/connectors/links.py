"""Deep links to resources that have no open, stable API.

KNApSAcK, IMPPAT, Dr. Duke's, PlantCyc, Plant Reactome, Phytozome, NP-MRD,
GNPS and MassBank are linked rather than scraped: scraping is fragile and
some sites forbid it. Each link opens the right record or a pre-filled search.
"""

from __future__ import annotations

from urllib.parse import quote, quote_plus

from .. import config
from ..models import Compound


def _fmt(key: str, **kw: str) -> str:
    return config.LINKS[key].format(**kw)


def compound_links(c: Compound, species: str = "") -> dict[str, str]:
    name_q = quote_plus(c.name)
    links: dict[str, str] = {}
    if c.cid:
        links["PubChem"] = _fmt("pubchem_compound", cid=str(c.cid))
    if c.xrefs.get("wikidata") or c.qid:
        links["Wikidata / LOTUS"] = _fmt("lotus_wikidata", qid=c.xrefs.get("wikidata") or c.qid)
    if c.xrefs.get("knapsack"):
        links["KNApSAcK"] = _fmt("knapsack_record", knapsack=c.xrefs["knapsack"])
    else:
        links["KNApSAcK (search)"] = _fmt("knapsack_search", name=name_q)
    if c.inchikey:
        links["COCONUT"] = _fmt("coconut_search", inchikey=c.inchikey)
    if c.xrefs.get("chebi"):
        links["ChEBI"] = f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:{c.xrefs['chebi']}"
    if c.xrefs.get("kegg"):
        links["KEGG"] = _fmt("kegg_compound", kegg=c.xrefs["kegg"])
    term = f'"{c.name}"' + (f' AND "{species}"' if species else "")
    links["PubMed"] = _fmt("pubmed_search", query=quote_plus(term))
    return links


def spectra_links(c: Compound) -> dict[str, str]:
    out = {}
    if c.cid:
        out["PubChem spectral information"] = _fmt("pubchem_spectra", cid=str(c.cid))
    if c.inchikey:
        q = f"exists(compound.metaData.name:'InChIKey' and compound.metaData.value:'{c.inchikey}')"
        out["MoNA (MS/MS, all records)"] = _fmt("mona_inchikey", query=quote(q, safe=""))
        out["MassBank Europe"] = _fmt("massbank_inchikey", inchikey=c.inchikey)
        out["NIST Chemistry WebBook (IR, UV-Vis, EI-MS)"] = _fmt("nist_webbook", inchikey=c.inchikey)
    out["NP-MRD (NMR)"] = _fmt("npmrd_search", name=quote_plus(c.name))
    out["GNPS spectral library"] = _fmt("gnps_library")
    return out


def uses_links(c: Compound) -> dict[str, str]:
    return {
        "IMPPAT (Indian medicinal plants)": _fmt("imppat"),
        "Dr. Duke's Phytochemical & Ethnobotanical DB": _fmt("dr_duke"),
        "PubMed reviews": _fmt("pubmed_search", query=quote_plus(f'"{c.name}" AND review[pt]')),
    }


def pathway_links(c: Compound) -> dict[str, str]:
    return {
        "PlantCyc / PMN": _fmt("plantcyc_search", name=quote(c.name)),
        "Plant Reactome": _fmt("plant_reactome_search", name=quote_plus(c.name)),
    }


def gene_links(ec: str) -> dict[str, str]:
    return {
        "KEGG enzyme": _fmt("kegg_ec", ec=f"ec:{ec}"),
        "BRENDA": _fmt("brenda_ec", ec=ec),
        "Phytozome": _fmt("phytozome"),
    }
