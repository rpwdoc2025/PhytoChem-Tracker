"""Unit tests: router, chemistry, spectra and API-response parsers. No network needed."""

import pytest

from phytochem_tracker import chem, router, spectra
from phytochem_tracker.audit import redact
from phytochem_tracker.connectors import kegg, links, mona, pubchem, pubmed, uniprot
from phytochem_tracker.models import Compound

QUERCETIN = "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O"


# ---------------------------------------------------------------- router
@pytest.mark.parametrize("value,itype,idtype", [
    ("Camellia sinensis", "species", "binomial"),
    ("Mentha x piperita", "species", "binomial"),
    ("Camellia sinensis var. assamica", "species", "binomial"),
    ("quercetin", "compound", "name"),
    ("REFJWTPEDVJJIY-UHFFFAOYSA-N", "compound", "inchikey"),
    ("InChI=1S/CH4/h1H4", "compound", "inchi"),
    (QUERCETIN, "compound", "smiles"),
    ("green tea", "compound", "name"),
])
def test_router_autodetect(value, itype, idtype):
    r = router.route(value)
    assert (r.input_type, r.id_type) == (itype, idtype)
    assert r.reason


def test_router_modes():
    assert router.route("peppermint", mode="species").input_type == "species"
    assert router.route("Camellia sinensis", mode="compound").input_type == "compound"
    assert router.route(filename="a.mzML").id_type == "mzml"
    with pytest.raises(ValueError):
        router.route("   ")


# ---------------------------------------------------------------- chemistry
def test_masses_and_descriptors():
    assert chem.monoisotopic_mass(QUERCETIN) == pytest.approx(302.042653, abs=1e-5)
    assert chem.inchikey_from_smiles(QUERCETIN) == "REFJWTPEDVJJIY-UHFFFAOYSA-N"
    d = chem.descriptors(QUERCETIN)
    assert d["Formula (RDKit)"] == "C15H10O7"
    assert d["H-bond donors"] == 5 and d["Lipinski violations"] == 0


def test_structure_files():
    assert "<svg" in chem.svg_2d(QUERCETIN)
    mb3 = chem.molblock_3d("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "caffeine")
    assert "V2000" in mb3 and mb3.startswith("caffeine")
    sdf = chem.sdf_from_molblock(chem.molblock_2d(QUERCETIN), {"CID": 5280343})
    assert sdf.rstrip().endswith("$$$$") and "<CID>" in sdf


def test_adducts():
    m = 302.042653
    assert chem.adduct_mz(m, chem.ADDUCTS_POS[0]) == pytest.approx(303.04993, abs=1e-4)  # [M+H]+
    assert chem.adduct_mz(m, chem.ADDUCTS_NEG[0]) == pytest.approx(301.03538, abs=1e-4)  # [M-H]-


# ---------------------------------------------------------------- spectra
MGF = """BEGIN IONS
TITLE=feature 1
PEPMASS=303.0500 1000
CHARGE=1+
153.018 100
229.050 30
END IONS
BEGIN IONS
TITLE=neg
PEPMASS=301.0354
CHARGE=1-
151.003 100
END IONS
"""


def test_parse_mgf():
    sp = spectra.parse_mgf(MGF)
    assert len(sp) == 2
    assert sp[0].precursor_mz == pytest.approx(303.05) and sp[0].polarity == "positive"
    assert sp[1].polarity == "negative" and sp[1].peaks == [(151.003, 100.0)]


def test_parse_msp():
    msp = "NAME: Quercetin\nPRECURSORMZ: 303.05\nION_MODE: P\nNum Peaks: 2\n153.018 100\n229.05 30\n\n"
    sp = spectra.parse_msp(msp)
    assert sp[0].title == "Quercetin" and len(sp[0].peaks) == 2 and sp[0].polarity == "positive"


def test_precursor_matching_and_cosine():
    cands = [{"name": "quercetin", "monoisotopic_mass": 302.042653}, {"name": "other", "monoisotopic_mass": 290.079}]
    hits = spectra.match_precursor(303.0500, cands, "positive", ppm=5)
    assert hits[0]["name"] == "quercetin" and hits[0]["adduct"] == "[M+H]+" and abs(hits[0]["ppm_error"]) < 1
    assert spectra.match_precursor(301.0354, cands, "negative", 5)[0]["adduct"] == "[M-H]-"
    a = [(153.018, 100), (229.05, 40), (257.045, 30)]
    assert spectra.cosine_similarity(a, a)[0] == pytest.approx(1.0)
    assert spectra.cosine_similarity(a, [(400.0, 100)])[0] == 0.0
    assert 0 < spectra.cosine_similarity(a, [(153.02, 100), (500.0, 100)])[0] < 1


def test_peak_string():
    assert spectra.parse_peak_string("153.018:100 229.05:30 junk") == [(153.018, 100.0), (229.05, 30.0)]


# ---------------------------------------------------------------- API parsers
PUGVIEW = {"Record": {
    "Reference": [{"ReferenceNumber": 7, "SourceName": "HSDB", "URL": "https://example.org/hsdb"}],
    "Section": [{"TOCHeading": "Chemical and Physical Properties", "Section": [
        {"TOCHeading": "Experimental Properties", "Section": [
            {"TOCHeading": "Melting Point", "Information": [
                {"ReferenceNumber": 7, "Value": {"StringWithMarkup": [{"String": "316.5 °C"}]}},
                {"ReferenceNumber": 7, "Value": {"Number": [316.5], "Unit": "°C"}}]},
            {"TOCHeading": "Physical Description", "Information": [
                {"ReferenceNumber": 7, "Value": {"StringWithMarkup": [{"String": "Yellow needles"}]}}]},
        ]}]}]}}


def test_pugview_experimental():
    props = pubchem.experimental_properties(PUGVIEW)
    assert props[0].name == "Physical Description"  # display order
    mp = [p for p in props if p.name == "Melting Point"]
    assert {p.value for p in mp} == {"316.5 °C"} and mp[0].source == "HSDB"  # numeric duplicate collapsed
    assert pubchem.solid_state(props).startswith("crystalline")


def test_pubchem_smiles_rename_compat():
    """PubChem renamed SMILES properties in July 2025; both shapes must parse."""
    legacy = {"CID": 1, "CanonicalSMILES": "CCO", "Title": "ethanol"}
    current = {"CID": 1, "SMILES": "CCO", "ConnectivitySMILES": "CCO", "Title": "ethanol"}
    assert pubchem.compound_from_props(legacy).smiles == pubchem.compound_from_props(current).smiles == "CCO"


def test_use_categories():
    assert pubchem.categorise_use("Therapeutic Uses", "x") == "Therapeutic / clinical"
    assert pubchem.categorise_use("Uses", "Used as an insecticide") == "Agricultural / veterinary"
    assert pubchem.categorise_use("Uses", "Food additive and flavoring") == "Industrial / commercial"
    assert pubchem.categorise_use("MeSH Pharmacological Classification", "Antioxidants") == "Pharmacological class"


def test_literature_synonyms_filter():
    syns = ["Quercetin", "117-39-5", "quercetin", "CHEBI:16243", "Sophoretin",
            "2-(3,4-dihydroxyphenyl)-3,5,7-trihydroxychromen-4-one", "UNII-9IKM0I5T1E"]
    assert pubchem.literature_synonyms(syns, 5) == ["Quercetin", "Sophoretin"]


def test_pubmed_query():
    q = pubmed.build_query(["quercetin", "sophoretin"], "Camellia sinensis")
    assert q == '("quercetin"[tiab] OR "sophoretin"[tiab]) AND "Camellia sinensis"[tiab]'
    assert pubmed.build_query(["caffeine"], "Coffea arabica") == '"caffeine"[tiab] AND "Coffea arabica"[tiab]'


def test_uniprot_query():
    assert uniprot.build_query("1.14.20.6", taxid="4442") == "(ec:1.14.20.6) AND (taxonomy_id:4442)"
    assert 'organism_name:"Camellia sinensis"' in uniprot.build_query("1.1.1.1", species="Camellia sinensis")


def test_mona_parse():
    recs = mona.parse_records([{"id": "X1", "spectrum": "153.0:100 229.0:30",
                                "compound": [{"names": [{"name": "Quercetin"}]}],
                                "metaData": [{"name": "ms level", "value": "MS2"},
                                             {"name": "precursor m/z", "value": 303.05}]}])
    assert recs[0]["ms_level"] == "MS2" and recs[0]["precursor_mz"] == "303.05" and len(recs[0]["peaks"]) == 2
    assert mona.parse_records({"content": []}) == []


def test_links_are_url_safe():
    c = Compound(name="(-)-Epicatechin", inchikey="PFTAWBLQPZVEMU-UKRRQHHQSA-N", cid=72276)
    for url in links.spectra_links(c).values():
        assert " " not in url and "'" not in url


def test_kegg_global_maps():
    assert "map01100" in kegg.GLOBAL_MAPS and "map00941" not in kegg.GLOBAL_MAPS


def test_redact():
    assert "SECRET" not in redact("https://x.org/a?term=q&api_key=SECRET&b=1")
