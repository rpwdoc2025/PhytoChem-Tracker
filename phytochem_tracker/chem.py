"""Local cheminformatics with RDKit: descriptors, depictions, files, masses.

These calculations run on the server, need no network, and are deterministic,
so they double as a cross-check for the values PubChem returns.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Any

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, Draw, Lipinski, QED, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

RDLogger.DisableLog("rdApp.*")

PROTON = 1.007276467  # mass of a proton (Da)

# Common ESI adducts: (label, charge, multiplier of M, mass delta in Da)
ADDUCTS_POS = [
    ("[M+H]+", 1, 1, PROTON),
    ("[M+Na]+", 1, 1, 22.989218),
    ("[M+K]+", 1, 1, 38.963158),
    ("[M+NH4]+", 1, 1, 18.033823),
    ("[M+H-H2O]+", 1, 1, PROTON - 18.010565),
    ("[2M+H]+", 1, 2, PROTON),
]
ADDUCTS_NEG = [
    ("[M-H]-", -1, 1, -PROTON),
    ("[M+Cl]-", -1, 1, 34.969402),
    ("[M+FA-H]-", -1, 1, 44.998201),
    ("[M-H2O-H]-", -1, 1, -PROTON - 18.010565),
    ("[2M-H]-", -1, 2, -PROTON),
]


@lru_cache(maxsize=2048)
def mol_from_smiles(smiles: str) -> Chem.Mol | None:
    if not smiles:
        return None
    return Chem.MolFromSmiles(smiles)


def mol_from_any(smiles: str = "", inchi: str = "") -> Chem.Mol | None:
    mol = mol_from_smiles(smiles) if smiles else None
    if mol is None and inchi:
        mol = Chem.MolFromInchi(inchi)
    return mol


def monoisotopic_mass(smiles: str) -> float | None:
    mol = mol_from_smiles(smiles)
    return round(Descriptors.ExactMolWt(mol), 6) if mol else None


def descriptors(smiles: str = "", inchi: str = "") -> dict[str, Any]:
    """Drug-likeness / isolation-relevant descriptors computed locally."""
    mol = mol_from_any(smiles, inchi)
    if mol is None:
        return {}
    hbd, hba = Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol)
    mw, logp = Descriptors.MolWt(mol), Crippen.MolLogP(mol)
    violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    try:
        qed = round(QED.qed(mol), 3)
    except Exception:  # noqa: BLE001
        qed = None
    return {
        "Formula (RDKit)": rdMolDescriptors.CalcMolFormula(mol),
        "Average MW (RDKit, g/mol)": round(mw, 3),
        "Monoisotopic mass (Da)": round(Descriptors.ExactMolWt(mol), 5),
        "Crippen cLogP": round(logp, 2),
        "TPSA (Å²)": round(rdMolDescriptors.CalcTPSA(mol), 1),
        "H-bond donors": hbd,
        "H-bond acceptors": hba,
        "Rotatable bonds": Lipinski.NumRotatableBonds(mol),
        "Aromatic rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "Fraction sp3 carbons": round(rdMolDescriptors.CalcFractionCSP3(mol), 2),
        "Stereocentres": len(Chem.FindMolChiralCenters(mol, includeUnassigned=True)),
        "Lipinski violations": violations,
        "QED drug-likeness": qed,
    }


def svg_2d(smiles: str, width: int = 420, height: int = 300) -> str:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return ""
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    opts.clearBackground = False
    opts.addStereoAnnotation = True
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def png_2d(smiles: str, size: tuple[int, int] = (420, 300)) -> bytes | None:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return None
    img = Draw.MolToImage(mol, size=size)
    import io

    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


def molblock_2d(smiles: str, name: str = "") -> str:
    mol = mol_from_smiles(smiles)
    if mol is None:
        return ""
    m = Chem.Mol(mol)
    AllChem.Compute2DCoords(m)
    if name:
        m.SetProp("_Name", name)
    return Chem.MolToMolBlock(m)


def molblock_3d(smiles: str, name: str = "", seed: int = 42) -> str:
    """Embed a 3D conformer (ETKDGv3 + MMFF) when PubChem has no 3D record."""
    mol = mol_from_smiles(smiles)
    if mol is None:
        return ""
    m = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(m, params) != 0:
        return ""
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    except Exception:  # noqa: BLE001
        pass
    if name:
        m.SetProp("_Name", name)
    return Chem.MolToMolBlock(m)


def sdf_from_molblock(molblock: str, props: dict[str, Any] | None = None) -> str:
    if not molblock:
        return ""
    lines = [molblock.rstrip("\n")]
    for k, v in (props or {}).items():
        if v not in (None, ""):
            lines.append(f">  <{k}>\n{v}\n")
    lines.append("$$$$\n")
    return "\n".join(lines)


def inchikey_from_smiles(smiles: str) -> str:
    mol = mol_from_smiles(smiles)
    return Chem.MolToInchiKey(mol) if mol else ""


def adduct_mz(neutral_mass: float, adduct: tuple[str, int, int, float]) -> float:
    _, charge, mult, delta = adduct
    return (mult * neutral_mass + delta) / abs(charge)
