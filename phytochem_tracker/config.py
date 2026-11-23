"""Central configuration.

All endpoints and link templates live here so that, when a database moves an
endpoint (it happens), the fix is a one-line change. Every value can be
overridden with an environment variable of the same name prefixed ``PCT_``
(e.g. ``PCT_PUBCHEM_BASE``).
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path


def _env(name: str, default: str) -> str:
    return os.getenv(f"PCT_{name}", default)


# --------------------------------------------------------------------------
# API endpoints (live, machine-readable)
# --------------------------------------------------------------------------
PUBCHEM_BASE = _env("PUBCHEM_BASE", "https://pubchem.ncbi.nlm.nih.gov/rest/pug")
PUBCHEM_VIEW_BASE = _env("PUBCHEM_VIEW_BASE", "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view")
EUTILS_BASE = _env("EUTILS_BASE", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils")
WIKIDATA_SPARQL = _env("WIKIDATA_SPARQL", "https://query.wikidata.org/sparql")
# QLever is a fast public Wikidata mirror used as an automatic fallback.
WIKIDATA_SPARQL_FALLBACK = _env("WIKIDATA_SPARQL_FALLBACK", "https://qlever.cs.uni-freiburg.de/api/wikidata")
KEGG_BASE = _env("KEGG_BASE", "https://rest.kegg.jp")
UNIPROT_BASE = _env("UNIPROT_BASE", "https://rest.uniprot.org")
MONA_BASE = _env("MONA_BASE", "https://mona.fiehnlab.ucdavis.edu/rest")

# --------------------------------------------------------------------------
# Deep-link templates for resources without an open API.
# Placeholders: {name} {inchikey} {cid} {species} {kegg} {knapsack} {ec}
# {gene} {taxid} {query}
# --------------------------------------------------------------------------
LINKS = {
    "pubchem_compound": "https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
    "pubchem_spectra": "https://pubchem.ncbi.nlm.nih.gov/compound/{cid}#section=Spectral-Information",
    "pubchem_png": PUBCHEM_BASE + "/compound/cid/{cid}/PNG?image_size=400x300",
    "mona_inchikey": "https://mona.fiehnlab.ucdavis.edu/spectra/browse?query={query}",
    "massbank_inchikey": "https://massbank.eu/MassBank/search?inchikey={inchikey}",
    "gnps_library": "https://library.gnps2.org/",
    "gnps_dashboard": "https://dashboard.gnps2.org/",
    "npmrd_search": "https://np-mrd.org/unearth/q?searcher=natural_products&query={name}",
    "nist_webbook": "https://webbook.nist.gov/cgi/cbook.cgi?InChI={inchikey}&Units=SI",
    "knapsack_record": "http://www.knapsackfamily.com/knapsack_core/information.php?word={knapsack}",
    "knapsack_search": "http://www.knapsackfamily.com/knapsack_core/result.php?sname=all&word={name}",
    "imppat": "https://cb.imsc.res.in/imppat/",
    "dr_duke": "https://phytochem.nal.usda.gov/",
    "coconut_search": "https://coconut.naturalproducts.net/search?q={inchikey}",
    "lotus_wikidata": "https://www.wikidata.org/wiki/{qid}",
    "kegg_compound": "https://www.kegg.jp/entry/{kegg}",
    "kegg_pathway": "https://www.kegg.jp/pathway/{pathway}+{kegg}",
    "kegg_ec": "https://www.kegg.jp/entry/{ec}",
    "plantcyc_search": "https://pmn.plantcyc.org/PLANT/substring-search?type=NIL&object={name}",
    "plant_reactome_search": "https://plantreactome.gramene.org/content/query?q={name}",
    "phytozome": "https://phytozome-next.jgi.doe.gov/",
    "uniprot_entry": "https://www.uniprot.org/uniprotkb/{acc}/entry",
    "brenda_ec": "https://www.brenda-enzymes.org/enzyme.php?ecno={ec}",
    "pubmed_search": "https://pubmed.ncbi.nlm.nih.gov/?term={query}",
}

# --------------------------------------------------------------------------
# Polite per-host request budgets (requests per second)
# --------------------------------------------------------------------------
RATE_LIMITS = {
    "pubchem.ncbi.nlm.nih.gov": 4.0,  # documented limit is 5/s
    "eutils.ncbi.nlm.nih.gov": 2.5,   # 3/s without key, raised to 9 with key
    "query.wikidata.org": 1.0,
    "qlever.cs.uni-freiburg.de": 1.0,
    "rest.kegg.jp": 3.0,
    "rest.uniprot.org": 5.0,
    "mona.fiehnlab.ucdavis.edu": 2.0,
}

USER_AGENT = _env(
    "USER_AGENT",
    "PhytoChem-Tracker/3.0 (https://github.com/rpwdoc2025/PhytoChem-Tracker; research use)",
)


@dataclass
class Settings:
    """Run-time knobs exposed in the sidebar."""

    n_top: int = 10                  # size of the ranked list
    n_candidates: int = 30           # LOTUS candidates sent to PubMed ranking
    synonyms_per_compound: int = 2   # extra PubChem synonyms in the PubMed query
    gene_compounds: int = 3          # compounds whose genes are mapped eagerly
    max_ec_per_compound: int = 8
    include_infraspecific: bool = True
    ppm_tolerance: float = 10.0      # spectrum dereplication
    fragment_tolerance: float = 0.02 # Da, for cosine similarity
    ncbi_api_key: str = field(default_factory=lambda: os.getenv("NCBI_API_KEY", ""))
    ncbi_email: str = field(default_factory=lambda: os.getenv("NCBI_EMAIL", ""))
    llm_provider: str = field(default_factory=lambda: os.getenv("PCT_LLM_PROVIDER", "none"))
    llm_api_key: str = ""
    llm_model: str = ""
    use_cache: bool = True
    cache_hours: int = 72
    offline_demo: bool = False


def cache_dir() -> Path:
    p = Path(_env("CACHE_DIR", str(Path.home() / ".cache" / "phytochem_tracker")))
    p.mkdir(parents=True, exist_ok=True)
    return p


def log_dir() -> Path:
    p = Path(_env("LOG_DIR", str(cache_dir() / "audit")))
    p.mkdir(parents=True, exist_ok=True)
    return p
