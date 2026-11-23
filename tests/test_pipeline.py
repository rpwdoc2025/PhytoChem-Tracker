"""End-to-end tests of the orchestrator and exports, using the offline DemoSession."""

import io
import zipfile

import pytest

from phytochem_tracker import report as rep
from phytochem_tracker.audit import AuditLog
from phytochem_tracker.config import Settings
from phytochem_tracker.demo import DEMO_MGF
from phytochem_tracker.orchestrator import Orchestrator
from phytochem_tracker.router import route


def orch():
    return Orchestrator(Settings(offline_demo=True), AuditLog(write_file=False))


@pytest.fixture(scope="module")
def species_report():
    o = orch()
    return o.run(route("Camellia sinensis")), o.audit


def test_species_ranking(species_report):
    r, audit = species_report
    assert not r.errors
    assert r.taxon.name == "Camellia sinensis" and r.taxon.ncbi_taxid == "4442"
    names = [rc.compound.name for rc in r.top]
    assert names[0] == "Caffeine"                      # highest (demo) PubMed count
    assert [rc.paper_count for rc in r.top] == sorted([rc.paper_count for rc in r.top], reverse=True)
    assert all("[tiab]" in rc.pubmed_query for rc in r.top)
    assert not any(n.endswith("-N") for n in names)    # unnamed LOTUS entries are skipped
    assert audit.stats()["errors"] == 0


def test_species_tiles(species_report):
    r, _ = species_report
    q = next(d for d in r.details.values() if d.compound.name == "Quercetin")
    assert q.computed["PubChem (computed)"]["Molecular formula"] == "C15H10O7"
    assert any(p.name == "Melting Point" and "316.5" in p.value for p in q.experimental)
    assert {"1H NMR", "13C NMR", "IR", "UV-Vis"} <= {s.type for s in q.spectra}
    assert any(s.peaks for s in q.spectra)             # MoNA reference spectrum
    assert "Pharmacological class" in {u.category for u in q.uses}
    assert q.structure["svg"] and q.structure["sdf3d"]
    flav = next(p for p in q.pathways if p.name == "Flavonoid biosynthesis")
    assert "1.14.20.6" in flav.enzymes
    assert q.genes and all(g.ec for g in q.genes)
    assert q.provenance["physchem"][0].source == "PubChem PUG-View"


def test_compound_mode_and_occurrences():
    r = orch().run(route("quercetin"))
    assert not r.errors and len(r.details) == 1
    assert r.occurrences[0]["taxon"] == "Camellia sinensis"
    assert r.notes  # explains missing species context


def test_inchikey_with_species_context():
    r = orch().run(route("REFJWTPEDVJJIY-UHFFFAOYSA-N"), context_species="Camellia sinensis")
    d = next(iter(r.details.values()))
    assert d.compound.cid == 5280343 and d.genes


def test_common_name_resolves():
    r = orch().run(route("tea", mode="species"))
    assert r.taxon.name == "Camellia sinensis" and r.taxon.resolved_from == "tea"


def test_unknown_species_is_reported_not_raised():
    r = orch().run(route("Nonexistent plantus"))
    assert r.errors and not r.details


def test_spectrum_dereplication():
    r = orch().run(route(filename="demo.mgf", mode="spectrum"), DEMO_MGF.encode(), "Camellia sinensis")
    assert not r.errors
    best = r.spectrum_matches[0]
    assert best["name"] == "quercetin" and best["adduct"] == "[M+H]+" and abs(best["ppm_error"]) < 2
    assert best["cosine"] > 0.9
    assert any(m["name"] == "caffeine" for m in r.spectrum_matches)


def test_spectrum_needs_species():
    r = orch().run(route(filename="demo.mgf", mode="spectrum"), DEMO_MGF.encode(), "")
    assert any("species" in e for e in r.errors)


def test_exports(species_report):
    r, audit = species_report
    assert rep.to_pdf(r)[:5] == b"%PDF-"
    assert b'"audit"' in rep.to_json(r, audit)
    names = zipfile.ZipFile(io.BytesIO(rep.to_csv_zip(r))).namelist()
    assert {"top_compounds.csv", "properties.csv", "genes.csv", "structures_2d.sdf"} <= set(names)
