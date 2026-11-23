"""PubChem PUG-REST (identifiers, computed properties, synonyms, files) and
PUG-View (experimental properties, spectral annotations, uses).

Note: in July 2025 PubChem renamed SMILES properties
(CanonicalSMILES -> ConnectivitySMILES, IsomericSMILES -> SMILES). The parser
accepts both the new and legacy names, and ``properties`` falls back to the
legacy names if the server rejects the new ones.
"""

from __future__ import annotations

import re
from typing import Any
from urllib.parse import quote

from .. import config
from ..http import Client, HTTPError
from ..models import Compound, PropertyValue

SOURCE = "PubChem PUG-REST"
VIEW_SOURCE = "PubChem PUG-View"

PROPS = [
    "Title", "MolecularFormula", "MolecularWeight", "MonoisotopicMass", "ExactMass",
    "SMILES", "ConnectivitySMILES", "InChI", "InChIKey", "IUPACName", "XLogP", "TPSA",
    "HBondDonorCount", "HBondAcceptorCount", "RotatableBondCount", "HeavyAtomCount",
    "Complexity", "Charge", "AtomStereoCount",
]
LEGACY_PROPS = [p for p in PROPS if p not in ("SMILES", "ConnectivitySMILES")] + ["IsomericSMILES", "CanonicalSMILES"]

# Display labels for the computed-property table
PROP_LABELS = {
    "MolecularFormula": "Molecular formula",
    "MolecularWeight": "Molecular weight (g/mol)",
    "MonoisotopicMass": "Monoisotopic mass (Da)",
    "XLogP": "XLogP3 (lipophilicity)",
    "TPSA": "TPSA (Å²)",
    "HBondDonorCount": "H-bond donors",
    "HBondAcceptorCount": "H-bond acceptors",
    "RotatableBondCount": "Rotatable bonds",
    "HeavyAtomCount": "Heavy atoms",
    "Complexity": "Complexity",
    "Charge": "Formal charge",
    "AtomStereoCount": "Stereocentres",
    "IUPACName": "IUPAC name",
    "InChIKey": "InChIKey",
}

# PUG-View headings shown in the physico-chemical tile, in display order
EXPERIMENTAL_HEADINGS = [
    "Physical Description", "Color/Form", "Melting Point", "Boiling Point", "Solubility",
    "Density", "LogP", "Vapor Pressure", "Decomposition", "Dissociation Constants",
    "Optical Rotation", "Stability/Shelf Life", "Odor", "Taste", "Kovats Retention Index",
    "Collision Cross Section",
]


def _smiles(rec: dict[str, Any]) -> str:
    for k in ("SMILES", "IsomericSMILES", "ConnectivitySMILES", "CanonicalSMILES"):
        if rec.get(k):
            return rec[k]
    return ""


def _float(v: Any) -> float | None:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


# --------------------------------------------------------------------------
# PUG-REST
# --------------------------------------------------------------------------
def cids_for(client: Client, value: str, id_type: str) -> list[int]:
    """Resolve name / smiles / inchi / inchikey to PubChem CIDs."""
    ns = {"name": "name", "inchikey": "inchikey", "smiles": "smiles", "inchi": "inchi"}[id_type]
    try:
        if ns in ("smiles", "inchi"):  # POST avoids URL-encoding problems
            data = client.request(
                "POST", f"{config.PUBCHEM_BASE}/compound/{ns}/cids/JSON", source=SOURCE,
                data={ns: value}, ok_statuses=(200, 404),
            ).json()
        else:
            data = client.request(
                "GET", f"{config.PUBCHEM_BASE}/compound/{ns}/{quote(value, safe='')}/cids/JSON",
                source=SOURCE, ok_statuses=(200, 404),
            ).json()
    except ValueError:
        return []
    cids = data.get("IdentifierList", {}).get("CID", [])
    return [int(c) for c in cids if int(c) > 0]


def properties(client: Client, cids: list[int]) -> dict[int, dict[str, Any]]:
    """Batch-fetch computed properties for up to ~100 CIDs in one call."""
    if not cids:
        return {}
    ids = ",".join(str(c) for c in cids)
    url = f"{config.PUBCHEM_BASE}/compound/cid/{ids}/property/{','.join(PROPS)}/JSON"
    try:
        data = client.get_json(url, source=SOURCE)
    except HTTPError as exc:
        if exc.status != 400:
            raise
        url = f"{config.PUBCHEM_BASE}/compound/cid/{ids}/property/{','.join(LEGACY_PROPS)}/JSON"
        data = client.get_json(url, source=SOURCE)
    out = {}
    for rec in data.get("PropertyTable", {}).get("Properties", []):
        out[int(rec["CID"])] = rec
    return out


def compound_from_props(rec: dict[str, Any], fallback_name: str = "") -> Compound:
    return Compound(
        name=rec.get("Title") or fallback_name or rec.get("IUPACName", ""),
        inchikey=rec.get("InChIKey", ""),
        cid=int(rec["CID"]) if rec.get("CID") else None,
        smiles=_smiles(rec),
        inchi=rec.get("InChI", ""),
        formula=rec.get("MolecularFormula", ""),
        mw=_float(rec.get("MolecularWeight")),
        monoisotopic_mass=_float(rec.get("MonoisotopicMass")),
    )


def computed_table(rec: dict[str, Any]) -> dict[str, Any]:
    out = {}
    for key, label in PROP_LABELS.items():
        if rec.get(key) not in (None, ""):
            out[label] = rec[key]
    if _smiles(rec):
        out["SMILES"] = _smiles(rec)
    return out


_JUNK_SYNONYM = re.compile(r"^(\d+-\d+-\d|[A-Z]{2,}\d|UNII|CHEBI|CHEMBL|DTXSID|NSC|SCHEMBL|ZINC|BRN|EINECS|HSDB|CCRIS|AKOS|MFCD)", re.I)


def synonyms(client: Client, cid: int, limit: int = 25) -> list[str]:
    try:
        data = client.get_json(f"{config.PUBCHEM_BASE}/compound/cid/{cid}/synonyms/JSON", source=SOURCE)
    except HTTPError:
        return []
    info = data.get("InformationList", {}).get("Information", [])
    syns = info[0].get("Synonym", []) if info else []
    return syns[:limit]


def literature_synonyms(names: list[str], n: int) -> list[str]:
    """Pick names suitable for a PubMed title/abstract query.

    Drops registry numbers and database codes, long systematic names, and
    duplicates (case-insensitive).
    """
    out: list[str] = []
    seen: set[str] = set()
    for s in names:
        s2 = s.strip()
        low = s2.lower()
        if (not s2 or low in seen or len(s2) > 40 or _JUNK_SYNONYM.match(s2)
                or sum(ch.isdigit() for ch in s2) > 3 or '"' in s2):
            continue
        seen.add(low)
        out.append(s2)
        if len(out) >= n:
            break
    return out


def sdf_url(cid: int, three_d: bool = False) -> str:
    suffix = "?record_type=3d" if three_d else ""
    return f"{config.PUBCHEM_BASE}/compound/cid/{cid}/SDF{suffix}"


def sdf(client: Client, cid: int, three_d: bool = False) -> str:
    try:
        return client.request("GET", sdf_url(cid, three_d), source=SOURCE, ok_statuses=(200,)).text
    except HTTPError:
        return ""


# --------------------------------------------------------------------------
# PUG-View
# --------------------------------------------------------------------------
def pug_view(client: Client, cid: int, heading: str) -> dict[str, Any] | None:
    url = f"{config.PUBCHEM_VIEW_BASE}/data/compound/{cid}/JSON?heading={quote(heading)}"
    resp = client.request("GET", url, source=VIEW_SOURCE, ok_statuses=(200, 404))
    if resp.status_code == 404:
        return None
    return resp.json()


def _value_text(value: dict[str, Any]) -> str:
    if "StringWithMarkup" in value:
        return " ".join(s.get("String", "") for s in value["StringWithMarkup"]).strip()
    if "Number" in value:
        nums = ", ".join(f"{n:g}" if isinstance(n, (int, float)) else str(n) for n in value["Number"])
        return f"{nums} {value.get('Unit', '')}".strip()
    if "Boolean" in value:
        return ", ".join(str(b) for b in value["Boolean"])
    if "DateISO8601" in value:
        return ", ".join(value["DateISO8601"])
    return ""


def iter_information(record: dict[str, Any] | None):
    """Yield (leaf_heading, parent_heading, info, reference) for every Information item."""
    if not record:
        return
    rec = record.get("Record", record)
    refs = {r.get("ReferenceNumber"): r for r in rec.get("Reference", [])}

    def walk(sections, parent=""):
        for sec in sections or []:
            heading = sec.get("TOCHeading", "")
            for info in sec.get("Information", []) or []:
                yield heading, parent, info, refs.get(info.get("ReferenceNumber"), {})
            yield from walk(sec.get("Section"), heading)

    yield from walk(rec.get("Section"))


def experimental_properties(record: dict[str, Any] | None, per_heading: int = 4) -> list[PropertyValue]:
    by_heading: dict[str, list[PropertyValue]] = {}
    for heading, _parent, info, ref in iter_information(record):
        text = _value_text(info.get("Value", {}))
        if not text:
            continue
        lst = by_heading.setdefault(heading, [])
        if len(lst) >= per_heading or any(p.value == text for p in lst):
            continue
        lst.append(PropertyValue(
            name=heading, value=text[:400], source=ref.get("SourceName", ""),
            reference=ref.get("URL", "") or info.get("URL", ""),
        ))
    order = {h: i for i, h in enumerate(EXPERIMENTAL_HEADINGS)}
    out: list[PropertyValue] = []
    for heading in sorted(by_heading, key=lambda h: order.get(h, 999)):
        out.extend(by_heading[heading])
    return out


_CRYSTAL = re.compile(r"crystal|needle|prism|plate|monoclinic|orthorhombic|rhombic|leaflet", re.I)
_AMORPH = re.compile(r"amorphous|glass|powder|resin|gum|wax", re.I)


def solid_state(props: list[PropertyValue]) -> str:
    """Heuristic crystalline/amorphous call from 'Physical Description'/'Color/Form' text."""
    texts = " ".join(p.value for p in props if p.name in ("Physical Description", "Color/Form"))
    if not texts:
        return ""
    c, a = bool(_CRYSTAL.search(texts)), bool(_AMORPH.search(texts))
    if c and not a:
        return "crystalline (reported)"
    if a and not c:
        return "amorphous / powder (reported)"
    if c and a:
        return "crystalline and powder forms reported"
    return "not stated"


SPECTRUM_TYPES = [
    ("1H NMR", "1H NMR"), ("13C NMR", "13C NMR"), ("2D NMR", "2D NMR"), ("NMR", "NMR"),
    ("MS-MS", "MS/MS"), ("LC-MS", "LC-MS"), ("GC-MS", "GC-MS"), ("Mass", "MS"),
    ("IR", "IR"), ("Raman", "Raman"), ("UV", "UV-Vis"),
]


def _spectrum_type(heading: str) -> str:
    for needle, label in SPECTRUM_TYPES:
        if needle.lower() in heading.lower():
            return label
    return heading


def spectral_annotations(record: dict[str, Any] | None) -> list[dict[str, Any]]:
    """Summarise PubChem's 'Spectral Information' section by spectrum type."""
    groups: dict[str, dict[str, Any]] = {}
    for heading, parent, info, ref in iter_information(record):
        g = groups.setdefault(heading, {
            "type": _spectrum_type(heading), "heading": heading, "parent": parent,
            "records": 0, "sources": set(), "urls": [], "text": [],
        })
        g["records"] += 1
        if ref.get("SourceName"):
            g["sources"].add(ref["SourceName"])
        url = info.get("URL") or ref.get("URL")
        if url and url not in g["urls"] and len(g["urls"]) < 3:
            g["urls"].append(url)
        text = _value_text(info.get("Value", {}))
        name = info.get("Name", "")
        if text and len(g["text"]) < 3 and len(text) < 200 and heading.lower().startswith(("uv", "spectral prop", "ir")):
            g["text"].append(f"{name}: {text}" if name else text)
    out = []
    for g in groups.values():
        g["sources"] = sorted(g["sources"])
        out.append(g)
    # Group headings repeat per sub-record; collapse identical ones
    return sorted(out, key=lambda g: (g["type"], -g["records"]))


USE_HEADINGS = ["Drug and Medication Information", "Use and Manufacturing", "Pharmacology and Biochemistry",
                "Agrochemical Information"]

_THERAPEUTIC = ("therapeutic uses", "drug indication", "drug classes", "clinical", "fda", "livertox")
_PHARM_CLASS = ("mesh pharmacological classification", "atc code", "pharmacological classes", "mechanism of action")
_SKIP = ("methods of manufacturing", "formulations/preparations", "general manufacturing", "sampling procedures",
         "u.s. production", "u.s. imports", "u.s. exports", "absorption, distribution and excretion",
         "metabolism/metabolites", "biological half-life", "pharmacokinetics", "human metabolite information",
         "interactions")
_AGRI = re.compile(r"pesticid|herbicid|insecticid|fungicid|acaricid|nematicid|veterinar|animal feed|livestock|"
                   r"agricultur|crop|plant growth|repellent|antifeedant", re.I)
_INDUSTRY = re.compile(r"food|flavo|cosmetic|dye|colou?r|fragranc|perfum|industr|manufactur|additive|"
                       r"antioxidant in|preservative|beverage|supplement|polymer|reagent|laborator|solvent", re.I)


def categorise_use(heading: str, text: str) -> str:
    h = heading.lower()
    if any(k in h for k in _THERAPEUTIC):
        return "Therapeutic / clinical"
    if any(k in h for k in _PHARM_CLASS):
        return "Pharmacological class"
    if _AGRI.search(text) or "agrochemical" in h:
        return "Agricultural / veterinary"
    if _INDUSTRY.search(text) or "industry" in h or "consumer" in h:
        return "Industrial / commercial"
    if "pharmacodynamics" in h or "pharmacology" in h:
        return "Bioactivity / pharmacology"
    return "Other reported uses"


def uses_from_records(records: list[dict[str, Any] | None], per_heading: int = 6) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    counts: dict[str, int] = {}
    seen: set[str] = set()
    for record in records:
        for heading, _parent, info, ref in iter_information(record):
            h = heading.lower()
            if any(s in h for s in _SKIP):
                continue
            text = _value_text(info.get("Value", {}))
            if heading.lower() == "mesh pharmacological classification" and info.get("Name"):
                text = f"{info['Name']}: {text}" if text else info["Name"]
            if not text or text in seen or counts.get(heading, 0) >= per_heading:
                continue
            seen.add(text)
            counts[heading] = counts.get(heading, 0) + 1
            out.append({
                "category": categorise_use(heading, text),
                "heading": heading,
                "text": text[:500],
                "source": ref.get("SourceName", "PubChem"),
                "url": ref.get("URL", ""),
            })
    return out
