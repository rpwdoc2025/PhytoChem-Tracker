"""Check every live connector against the real databases.

Run this on a machine with internet access after installing requirements:

    python scripts/live_smoke_test.py            # quick check of each source
    python scripts/live_smoke_test.py --full     # also runs a full species query

Each line prints PASS/FAIL with the reason, so an API change at any one
source is easy to spot. No API keys are needed.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from phytochem_tracker.audit import AuditLog  # noqa: E402
from phytochem_tracker.config import Settings  # noqa: E402
from phytochem_tracker.connectors import kegg, mona, pubchem, pubmed, uniprot, wikidata  # noqa: E402
from phytochem_tracker.http import Client  # noqa: E402

QUERCETIN_KEY = "REFJWTPEDVJJIY-UHFFFAOYSA-N"


def check(name, fn, expect):
    t0 = time.perf_counter()
    try:
        value = fn()
        ok = expect(value)
        msg = repr(value)[:110]
    except Exception as exc:  # noqa: BLE001
        ok, msg = False, f"{type(exc).__name__}: {exc}"[:160]
    print(f"{'PASS' if ok else 'FAIL'}  {name:<44} {time.perf_counter() - t0:5.1f}s  {msg}")
    return ok


def main() -> int:
    c = Client(AuditLog(write_file=False), use_cache=False)
    results = [
        check("PubChem: name -> CID", lambda: pubchem.cids_for(c, "quercetin", "name"), lambda v: 5280343 in v),
        check("PubChem: properties (2025 SMILES names)", lambda: pubchem.properties(c, [5280343])[5280343],
              lambda v: v.get("InChIKey") == QUERCETIN_KEY and bool(pubchem._smiles(v))),
        check("PubChem PUG-View: melting point", lambda: [p.value for p in pubchem.experimental_properties(
            pubchem.pug_view(c, 5280343, "Experimental Properties")) if p.name == "Melting Point"], bool),
        check("PubChem PUG-View: spectral info", lambda: [g["type"] for g in pubchem.spectral_annotations(
            pubchem.pug_view(c, 5280343, "Spectral Information"))], bool),
        check("PubMed: count", lambda: pubmed.count(c, pubmed.build_query(["caffeine"], "Camellia sinensis")),
              lambda v: v > 100),
        check("Wikidata: resolve 'Camellia sinensis'", lambda: wikidata.resolve_taxon(c, "Camellia sinensis"),
              lambda v: v and v["ncbi_taxid"] == "4442"),
        check("Wikidata: resolve common name 'peppermint'", lambda: wikidata.resolve_taxon(c, "peppermint"),
              lambda v: v and "Mentha" in v["name"]),
        check("Wikidata/LOTUS: compounds in C. sinensis", lambda: len(wikidata.compounds_in_taxon(
            c, wikidata.resolve_taxon(c, "Camellia sinensis")["qid"])), lambda v: v > 50),
        check("Wikidata: xrefs for quercetin", lambda: wikidata.compound_xrefs(c, [QUERCETIN_KEY]).get(QUERCETIN_KEY),
              lambda v: v and v.get("kegg") == "C00389"),
        check("KEGG: pathways for C00389", lambda: kegg.pathways_for_compound(c, "C00389"),
              lambda v: "map00941" in v or "map00944" in v),
        check("KEGG: enzymes for C00389", lambda: kegg.enzymes_for_compound(c, "C00389"), bool),
        check("KEGG: organism code for A. thaliana", lambda: kegg.organism_code(c, "Arabidopsis thaliana"),
              lambda v: v == "ath"),
        check("KEGG: genes ath / EC 1.14.20.6", lambda: kegg.genes_for_ec(c, "ath", "1.14.20.6"), bool),
        check("UniProt: FLS in Arabidopsis", lambda: uniprot.proteins_for_ec(c, "1.14.20.6", "3702"), bool),
        check("MoNA: spectra for quercetin", lambda: len(mona.spectra_for_inchikey(c, QUERCETIN_KEY, 5)),
              lambda v: v > 0),
    ]
    if "--full" in sys.argv:
        from phytochem_tracker.orchestrator import Orchestrator
        from phytochem_tracker.router import route

        def full():
            r = Orchestrator(Settings(n_candidates=15), AuditLog(write_file=False)).run(route("Camellia sinensis"))
            return [(rc.compound.name, rc.paper_count) for rc in r.top], r.errors

        results.append(check("Full species run (15 candidates)", full, lambda v: len(v[0]) >= 5 and not v[1]))
    print(f"\n{sum(results)}/{len(results)} checks passed.")
    return 0 if all(results) else 1


if __name__ == "__main__":
    sys.exit(main())
