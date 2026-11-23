"""Offline demo mode and test fixtures.

``DemoSession`` is a drop-in ``requests.Session`` that answers every external
API call with a small fixture shaped exactly like the real response (PubChem
JSON, Wikidata SPARQL JSON, KEGG tab-separated text, UniProt TSV, MoNA JSON).
It powers the "Offline demo" switch in the UI and the automated tests.

Identifiers, formulas and masses below are real; PubMed counts, UniProt
accessions and the MoNA record are ILLUSTRATIVE placeholders (prefixed DEMO
where they are identifiers) and the UI labels demo results accordingly.
"""

from __future__ import annotations

import json
import re
from typing import Any
from urllib.parse import parse_qs, unquote, urlparse

import requests

DEMO_SPECIES = "Camellia sinensis"

COMPOUNDS: dict[str, dict[str, Any]] = {
    "caffeine": {
        "CID": 2519, "Title": "Caffeine", "MolecularFormula": "C8H10N4O2", "MolecularWeight": "194.19",
        "MonoisotopicMass": "194.08037557", "SMILES": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "ConnectivitySMILES": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "InChI": "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
        "InChIKey": "RYYVLZVUVIJVGH-UHFFFAOYSA-N", "IUPACName": "1,3,7-trimethylpurine-2,6-dione",
        "XLogP": -0.1, "TPSA": 58.4, "HBondDonorCount": 0, "HBondAcceptorCount": 3, "RotatableBondCount": 0,
        "HeavyAtomCount": 14, "Complexity": 293, "Charge": 0,
        "qid": "Q60235", "kegg": "C07481", "chebi": "27732", "refs": 41, "papers": 2310,
        "synonyms": ["caffeine", "Caffeine", "1,3,7-Trimethylxanthine", "Guaranine", "Theine", "58-08-2"],
    },
    "quercetin": {
        "CID": 5280343, "Title": "Quercetin", "MolecularFormula": "C15H10O7", "MolecularWeight": "302.23",
        "MonoisotopicMass": "302.04265265", "SMILES": "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O",
        "ConnectivitySMILES": "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O",
        "InChI": "InChI=1S/C15H10O7/c16-7-4-10(19)12-11(5-7)22-15(14(21)13(12)20)6-1-2-8(17)9(18)3-6/h1-5,16-19,21H",
        "InChIKey": "REFJWTPEDVJJIY-UHFFFAOYSA-N",
        "IUPACName": "2-(3,4-dihydroxyphenyl)-3,5,7-trihydroxychromen-4-one",
        "XLogP": 1.5, "TPSA": 127, "HBondDonorCount": 5, "HBondAcceptorCount": 7, "RotatableBondCount": 1,
        "HeavyAtomCount": 22, "Complexity": 488, "Charge": 0,
        "qid": "Q409478", "kegg": "C00389", "chebi": "16243", "refs": 35, "papers": 1460,
        "synonyms": ["quercetin", "Quercetin", "Sophoretin", "Meletin", "117-39-5", "Quercetine"],
    },
    "epicatechin": {
        "CID": 72276, "Title": "(-)-Epicatechin", "MolecularFormula": "C15H14O6", "MolecularWeight": "290.27",
        "MonoisotopicMass": "290.07903816", "SMILES": "C1[C@H]([C@H](OC2=CC(=CC(=C21)O)O)C3=CC(=C(C=C3)O)O)O",
        "ConnectivitySMILES": "C1C(C(OC2=CC(=CC(=C21)O)O)C3=CC(=C(C=C3)O)O)O",
        "InChI": "InChI=1S/C15H14O6/c16-8-4-11(18)9-6-13(20)15(21-14(9)5-8)7-1-2-10(17)12(19)3-7/h1-5,13,15-20H,6H2/t13-,15-/m1/s1",
        "InChIKey": "PFTAWBLQPZVEMU-UKRRQHHQSA-N",
        "IUPACName": "(2R,3R)-2-(3,4-dihydroxyphenyl)-3,4-dihydro-2H-chromene-3,5,7-triol",
        "XLogP": 0.4, "TPSA": 110, "HBondDonorCount": 5, "HBondAcceptorCount": 6, "RotatableBondCount": 1,
        "HeavyAtomCount": 21, "Complexity": 364, "Charge": 0,
        "qid": "", "kegg": "", "chebi": "", "refs": 28, "papers": 1890,
        "synonyms": ["(-)-Epicatechin", "Epicatechin", "L-Epicatechin", "490-46-0"],
    },
    "theobromine": {
        "CID": 5429, "Title": "Theobromine", "MolecularFormula": "C7H8N4O2", "MolecularWeight": "180.16",
        "MonoisotopicMass": "180.06472551", "SMILES": "CN1C=NC2=C1C(=O)NC(=O)N2C",
        "ConnectivitySMILES": "CN1C=NC2=C1C(=O)NC(=O)N2C",
        "InChI": "InChI=1S/C7H8N4O2/c1-10-3-8-5-4(10)6(12)9-7(13)11(5)2/h3H,1-2H3,(H,9,12,13)",
        "InChIKey": "YAPQBXQYLJRXSA-UHFFFAOYSA-N", "IUPACName": "3,7-dimethylpurine-2,6-dione",
        "XLogP": -0.8, "TPSA": 67.2, "HBondDonorCount": 1, "HBondAcceptorCount": 3, "RotatableBondCount": 0,
        "HeavyAtomCount": 13, "Complexity": 268, "Charge": 0,
        "qid": "", "kegg": "C07480", "chebi": "", "refs": 22, "papers": 270,
        "synonyms": ["Theobromine", "theobromine", "3,7-Dimethylxanthine", "Xantheose", "83-67-0"],
    },
}
BY_CID = {v["CID"]: k for k, v in COMPOUNDS.items()}
BY_KEY = {v["InChIKey"]: k for k, v in COMPOUNDS.items()}
BY_NAME = {s.lower(): k for k, v in COMPOUNDS.items() for s in v["synonyms"] + [v["Title"]]}

EXPERIMENTAL = {
    "quercetin": [
        ("Physical Description", "Yellow needles or yellow powder.", 1),
        ("Color/Form", "Yellow needles from dilute alcohol", 2),
        ("Melting Point", "316.5 °C", 1),
        ("Solubility", "Practically insoluble in water; soluble in aqueous alkaline solutions", 2),
    ],
    "caffeine": [
        ("Physical Description", "Odorless white powder or white glistening needles, usually melted together.", 1),
        ("Color/Form", "Hexagonal prisms", 2),
        ("Melting Point", "238 °C", 1),
        ("Boiling Point", "178 °C (sublimes)", 2),
        ("Solubility", "In water, 2.16X10+4 mg/L at 25 °C", 2),
        ("LogP", "log Kow = -0.07", 2),
    ],
}
SPECTRAL = {
    "quercetin": [("1H NMR Spectra", 2), ("13C NMR Spectra", 2), ("LC-MS", 6), ("MS-MS", 4), ("FTIR Spectra", 1),
                  ("UV Spectra", 1)],
    "caffeine": [("1H NMR Spectra", 3), ("13C NMR Spectra", 3), ("GC-MS", 8), ("LC-MS", 10), ("FTIR Spectra", 2),
                 ("UV Spectra", 2)],
}
USES = {
    "caffeine": [
        ("Therapeutic Uses", "Central Nervous System Stimulants; used as an adjunct analgesic", 2),
        ("MeSH Pharmacological Classification", ("Central Nervous System Stimulants",
         "A loosely defined group of drugs that tend to increase behavioral alertness"), 3),
        ("Uses", "Flavoring agent in cola beverages; food additive", 2),
        ("Uses", "Insecticide and molluscicide (experimental)", 2),
    ],
    "quercetin": [
        ("MeSH Pharmacological Classification", ("Antioxidants",
         "Naturally occurring or synthetic substances that inhibit or retard oxidation"), 3),
        ("Uses", "Dietary supplement; antioxidant in food research", 2),
    ],
}
REFS = [
    {"ReferenceNumber": 1, "SourceName": "NTP", "URL": "https://ntp.niehs.nih.gov/"},
    {"ReferenceNumber": 2, "SourceName": "Hazardous Substances Data Bank (HSDB)", "URL": "https://pubchem.ncbi.nlm.nih.gov/source/hsdb"},
    {"ReferenceNumber": 3, "SourceName": "MeSH", "URL": "https://www.ncbi.nlm.nih.gov/mesh"},
]
KEGG = {
    "C00389": {"pathways": ["map00941", "map00944", "map01100", "map01110"],
               "enzymes": ["1.14.14.82", "1.14.20.6", "2.4.1.91"]},
    "C07481": {"pathways": ["map00232", "map01100"], "enzymes": ["2.1.1.160"]},
    "C07480": {"pathways": ["map00232", "map01100"], "enzymes": ["2.1.1.159", "2.1.1.160"]},
}
KEGG_PATHWAYS = {"map00941": "Flavonoid biosynthesis", "map00944": "Flavone and flavonol biosynthesis",
                 "map00232": "Caffeine metabolism", "map01100": "Metabolic pathways",
                 "map01110": "Biosynthesis of secondary metabolites"}
KEGG_PATHWAY_ENZYMES = {"map00941": ["1.14.20.6", "1.14.14.82"], "map00944": ["1.14.14.82", "2.4.1.91"],
                        "map00232": ["2.1.1.159", "2.1.1.160"]}
KEGG_EC = {"1.14.14.82": "flavonoid 3'-monooxygenase", "1.14.20.6": "flavonol synthase",
           "2.4.1.91": "flavonol 3-O-glucosyltransferase", "2.1.1.160": "caffeine synthase",
           "2.1.1.159": "theobromine synthase"}

# Illustrative MS/MS of quercetin [M+H]+ (characteristic RDA fragments 153, 137, 229, 257).
QUERCETIN_MS2 = "121.028:8 137.023:22 153.018:100 165.018:18 201.055:6 229.050:35 247.060:9 257.045:30 285.039:12 303.050:40"


def _resp(url: str, status: int = 200, body: Any = "", content_type: str = "application/json") -> requests.Response:
    r = requests.Response()
    r.status_code = status
    r.url = url
    r.encoding = "utf-8"
    r._content = (json.dumps(body) if not isinstance(body, (str, bytes)) else body)
    if isinstance(r._content, str):
        r._content = r._content.encode()
    r.headers["Content-Type"] = content_type
    r.from_cache = False  # type: ignore[attr-defined]
    return r


def _bindings(rows: list[dict[str, str]]) -> dict[str, Any]:
    vars_ = sorted({k for r in rows for k in r})
    return {"head": {"vars": vars_}, "results": {"bindings": [
        {k: {"type": "literal", "value": str(v)} for k, v in r.items() if v not in (None, "")} for r in rows]}}


def _pug_view_record(cid: int, sections: list[tuple[str, Any, int]], name_field: bool = False) -> dict[str, Any]:
    subs = []
    for heading, value, ref in sections:
        info: dict[str, Any] = {"ReferenceNumber": ref}
        if isinstance(value, tuple):
            info["Name"], text = value
        else:
            text = value
        info["Value"] = {"StringWithMarkup": [{"String": text}]}
        subs.append({"TOCHeading": heading, "Information": [info]})
    return {"Record": {"RecordType": "CID", "RecordNumber": cid, "Reference": REFS,
                       "Section": [{"TOCHeading": "Top", "Section": subs}]}}


class DemoSession(requests.Session):
    """Routes requests to fixtures. Unknown endpoints get a 404."""

    def request(self, method, url, params=None, data=None, headers=None, **kwargs):  # noqa: D401
        full = requests.Request(method, url, params=params).prepare().url or url
        p = urlparse(full)
        qs = {k: v[0] for k, v in parse_qs(p.query).items()}
        host, path = p.netloc, unquote(p.path)
        try:
            if "pubchem" in host:
                return self._pubchem(full, path, qs, data or {})
            if "eutils" in host:
                return self._eutils(full, qs)
            if "wikidata" in host or "qlever" in host:
                return self._sparql(full, qs.get("query", ""))
            if "kegg" in host:
                return self._kegg(full, path)
            if "uniprot" in host:
                return self._uniprot(full, qs)
            if "mona" in host:
                return self._mona(full, qs)
        except KeyError:
            pass
        return _resp(full, 404, {"Fault": {"Code": "PUGREST.NotFound"}})

    # ------------------------------------------------------------------
    def _pubchem(self, url, path, qs, data):
        m = re.search(r"/compound/(name|inchikey)/([^/]+)/cids/JSON", path)
        if m:
            key = BY_NAME.get(m.group(2).lower()) if m.group(1) == "name" else BY_KEY.get(m.group(2))
            if not key:
                return _resp(url, 404, {"Fault": {"Code": "PUGREST.NotFound"}})
            return _resp(url, 200, {"IdentifierList": {"CID": [COMPOUNDS[key]["CID"]]}})
        if re.search(r"/compound/(smiles|inchi)/cids/JSON", path):
            value = (data or {}).get("smiles") or (data or {}).get("inchi") or ""
            for v in COMPOUNDS.values():
                if value in (v["SMILES"], v["ConnectivitySMILES"], v["InChI"]):
                    return _resp(url, 200, {"IdentifierList": {"CID": [v["CID"]]}})
            return _resp(url, 404, {"Fault": {"Code": "PUGREST.NotFound"}})
        m = re.search(r"/compound/cid/([\d,]+)/property/", path)
        if m:
            props = []
            for cid in m.group(1).split(","):
                key = BY_CID.get(int(cid))
                if key:
                    props.append({k: v for k, v in COMPOUNDS[key].items()
                                  if k not in ("qid", "kegg", "chebi", "refs", "papers", "synonyms")})
            return _resp(url, 200, {"PropertyTable": {"Properties": props}})
        m = re.search(r"/compound/cid/(\d+)/synonyms/JSON", path)
        if m:
            key = BY_CID[int(m.group(1))]
            return _resp(url, 200, {"InformationList": {"Information": [
                {"CID": int(m.group(1)), "Synonym": COMPOUNDS[key]["synonyms"]}]}})
        m = re.search(r"/compound/cid/(\d+)/SDF", path)
        if m:
            from . import chem

            key = BY_CID[int(m.group(1))]
            mb = chem.molblock_3d(COMPOUNDS[key]["SMILES"], COMPOUNDS[key]["Title"])
            return _resp(url, 200, chem.sdf_from_molblock(mb), "chemical/x-mdl-sdfile")
        m = re.search(r"/pug_view/data/compound/(\d+)/JSON", path)
        if m:
            cid = int(m.group(1))
            key = BY_CID.get(cid, "")
            heading = qs.get("heading", "")
            if heading == "Experimental Properties" and key in EXPERIMENTAL:
                return _resp(url, 200, _pug_view_record(cid, EXPERIMENTAL[key]))
            if heading == "Spectral Information" and key in SPECTRAL:
                secs = []
                for h, n in SPECTRAL[key]:
                    for i in range(n):
                        secs.append((h, f"Spectrum record {i + 1}", 2))
                return _resp(url, 200, _pug_view_record(cid, secs))
            if heading in ("Drug and Medication Information", "Use and Manufacturing", "Pharmacology and Biochemistry") and key in USES:
                wanted = {"Drug and Medication Information": ("Therapeutic Uses",),
                          "Use and Manufacturing": ("Uses",),
                          "Pharmacology and Biochemistry": ("MeSH Pharmacological Classification",)}[heading]
                secs = [u for u in USES[key] if u[0] in wanted]
                if secs:
                    return _resp(url, 200, _pug_view_record(cid, secs))
            return _resp(url, 404, {"Fault": {"Code": "PUGVIEW.NotFound"}})
        return _resp(url, 404, {"Fault": {"Code": "PUGREST.NotFound"}})

    def _eutils(self, url, qs):
        term = qs.get("term", "").lower()
        n = 0
        for v in COMPOUNDS.values():
            if any(f'"{s.lower()}"' in term for s in v["synonyms"] + [v["Title"]]):
                n = v["papers"]
                break
        return _resp(url, 200, {"header": {"type": "esearch", "version": "0.3"},
                                "esearchresult": {"count": str(n), "retmax": "0", "retstart": "0", "idlist": []}})

    def _sparql(self, url, q):
        if "wdt:P225 ?name" in q and "VALUES ?name" in q:
            if "Camellia sinensis" in q:
                return _resp(url, 200, _bindings([{"t": "http://www.wikidata.org/entity/Q101815",
                                                   "name": "Camellia sinensis", "ncbi": "4442"}]))
            if "Arabidopsis thaliana" in q:
                return _resp(url, 200, _bindings([{"t": "http://www.wikidata.org/entity/Q158695",
                                                   "name": "Arabidopsis thaliana", "ncbi": "3702"}]))
            return _resp(url, 200, _bindings([]))
        if "wdt:P1843 ?common" in q:
            if '"tea"@en' in q:
                return _resp(url, 200, _bindings([{"t": "http://www.wikidata.org/entity/Q101815",
                                                   "name": "Camellia sinensis", "ncbi": "4442"}]))
            return _resp(url, 200, _bindings([]))
        if "wdt:P1843 ?c" in q:
            return _resp(url, 200, _bindings([{"c": "tea plant"}, {"c": "tea"}] if "Q101815" in q else []))
        if "GROUP BY ?c ?ik" in q:
            if "Q101815" not in q:
                return _resp(url, 200, _bindings([]))
            rows = []
            for v in COMPOUNDS.values():
                rows.append({"c": f"http://www.wikidata.org/entity/{v['qid'] or 'Q0'}", "ik": v["InChIKey"],
                             "label": v["Title"].lower() if v["Title"] != "(-)-Epicatechin" else "(−)-epicatechin",
                             "smiles_a": v["ConnectivitySMILES"], "smiles_b": v["SMILES"], "cid": v["CID"],
                             "refs": v["refs"]})
            rows.append({"c": "http://www.wikidata.org/entity/Q0", "ik": "XXXXXXXXXXXXXX-UHFFFAOYSA-N",
                         "label": "", "refs": 1})  # unlabeled, no structure: must be handled gracefully
            return _resp(url, 200, _bindings(rows))
        if "VALUES ?ik" in q:
            rows = []
            for v in COMPOUNDS.values():
                if v["InChIKey"] in q:
                    rows.append({"ik": v["InChIKey"], "item": f"http://www.wikidata.org/entity/{v['qid'] or 'Q0'}",
                                 "kegg_": v["kegg"], "chebi_": v["chebi"], "pubchem_cid_": v["CID"]})
            return _resp(url, 200, _bindings(rows))
        if "ORDER BY DESC(?refs)" in q:
            return _resp(url, 200, _bindings([
                {"t": "http://www.wikidata.org/entity/Q101815", "name": "Camellia sinensis", "refs": 12},
                {"t": "http://www.wikidata.org/entity/Q41500", "name": "Coffea arabica", "refs": 9},
                {"t": "http://www.wikidata.org/entity/Q1426208", "name": "Ilex paraguariensis", "refs": 4},
            ]))
        return _resp(url, 200, _bindings([]))

    def _kegg(self, url, path):
        parts = path.strip("/").split("/")
        op = parts[0]
        text = ""
        if op == "list" and parts[1] == "organism":
            text = "T00041\tath\tArabidopsis thaliana (thale cress)\tEukaryotes;Plants;Eudicots;Mustard family\n"
        elif op == "list" and parts[1] == "pathway":
            text = "".join(f"path:{k}\t{v}\n" for k, v in KEGG_PATHWAYS.items())
        elif op == "list" and parts[1].startswith("ec:"):
            for e in parts[1].split("+"):
                e = e.replace("ec:", "")
                if e in KEGG_EC:
                    text += f"ec:{e}\t{KEGG_EC[e]}\n"
        elif op == "list" and parts[1].startswith("cpd:"):
            for c in parts[1].split("+"):
                c = c.replace("cpd:", "")
                for v in COMPOUNDS.values():
                    if v["kegg"] == c:
                        text += f"cpd:{c}\t{v['Title']}; {v['synonyms'][2]}\n"
        elif op == "find" and len(parts) >= 4 and parts[3] == "formula":
            for v in COMPOUNDS.values():
                if v["kegg"] and v["MolecularFormula"] == parts[2]:
                    text += f"cpd:{v['kegg']}\t{v['MolecularFormula']}\n"
        elif op == "link" and parts[1] == "pathway":
            c = parts[2].replace("cpd:", "")
            text = "".join(f"cpd:{c}\tpath:{m}\n" for m in KEGG.get(c, {}).get("pathways", []))
        elif op == "link" and parts[1] == "enzyme" and parts[2].startswith("cpd:"):
            c = parts[2].replace("cpd:", "")
            text = "".join(f"cpd:{c}\tec:{e}\n" for e in KEGG.get(c, {}).get("enzymes", []))
        elif op == "link" and parts[1] == "enzyme" and parts[2].startswith("path:"):
            m = parts[2].replace("path:", "")
            text = "".join(f"path:{m}\tec:{e}\n" for e in KEGG_PATHWAY_ENZYMES.get(m, []))
        elif op == "link" and parts[1] == "ath":
            e = parts[2].replace("ec:", "")
            if e == "1.14.20.6":
                text = "ec:1.14.20.6\tath:AT5G08640\n"
        if not text:
            return _resp(url, 404, "", "text/plain")
        return _resp(url, 200, text, "text/plain")

    def _uniprot(self, url, qs):
        q = qs.get("query", "")
        m = re.search(r"ec:([\d.]+)", q)
        ec = m.group(1) if m else ""
        header = "Entry\tGene Names\tProtein names\tOrganism\tReviewed\n"
        name = KEGG_EC.get(ec)
        if not name:
            return _resp(url, 200, header, "text/plain")
        body = header + (f"DEMO{ec.replace('.', '')}\tDEMO-{ec}\t{name} (demo record)\t"
                         f"Camellia sinensis (Tea plant)\tunreviewed\n")
        return _resp(url, 200, body, "text/plain")

    def _mona(self, url, qs):
        q = qs.get("query", "")
        if COMPOUNDS["quercetin"]["InChIKey"] in q and q.startswith("exists("):
            return _resp(url, 200, [{
                "id": "DEMO-QUERCETIN-MS2", "spectrum": QUERCETIN_MS2,
                "compound": [{"names": [{"name": "Quercetin"}]}],
                "metaData": [{"name": "ms level", "value": "MS2"}, {"name": "precursor type", "value": "[M+H]+"},
                             {"name": "precursor m/z", "value": 303.0499}, {"name": "ionization mode", "value": "positive"},
                             {"name": "instrument type", "value": "LC-ESI-QTOF (demo)"},
                             {"name": "collision energy", "value": "20 eV"}],
            }])
        return _resp(url, 200, [])


DEMO_MGF = """BEGIN IONS
TITLE=Tea leaf extract feature 1 (demo: quercetin-like)
PEPMASS=303.0500 15000
CHARGE=1+
IONMODE=positive
RTINSECONDS=512.4
121.0290 7
137.0234 20
153.0180 100
165.0182 15
229.0499 33
257.0448 27
285.0390 10
END IONS
BEGIN IONS
TITLE=Tea leaf extract feature 2 (demo: caffeine-like)
PEPMASS=195.0877 42000
CHARGE=1+
IONMODE=positive
138.0662 100
110.0713 18
END IONS
"""
