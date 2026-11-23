"""Router agent: decide what kind of input the user gave.

The decision and its reason are logged, so a user can always see why the
orchestrator chose a given tool chain (slide 14: "no black box").
"""

from __future__ import annotations

import re
from dataclasses import dataclass

INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
# Binomial (optionally with hybrid sign and an infraspecific rank),
# e.g. "Camellia sinensis", "Mentha x piperita", "Camellia sinensis var. assamica".
BINOMIAL_RE = re.compile(
    r"^[A-Z][a-z]+(?:\s+(?:x|×))?\s+[a-z][a-z-]+"
    r"(?:\s+(?:subsp\.|ssp\.|var\.|f\.)\s+[a-z][a-z-]+)?$"
)
SPECTRUM_EXTENSIONS = (".mgf", ".mzml", ".msp")
# Characters that occur in SMILES but essentially never in compound names.
_SMILES_HINT = re.compile(r"[=#\[\]@\\/]|\(.*\)|^[BCNOPSFIcnosp][A-Za-z0-9]*\d")


@dataclass
class Route:
    input_type: str   # species | compound | spectrum
    id_type: str      # binomial | name | smiles | inchi | inchikey | mgf | mzml | msp
    value: str
    reason: str


def _looks_like_smiles(text: str) -> bool:
    if " " in text or len(text) < 2:
        return False
    if not _SMILES_HINT.search(text):
        return False
    try:
        from rdkit import Chem, RDLogger

        RDLogger.DisableLog("rdApp.*")
        return Chem.MolFromSmiles(text) is not None
    except ImportError:  # pragma: no cover - RDKit is a hard dependency in requirements
        return True


def route(value: str = "", filename: str | None = None, mode: str = "auto") -> Route:
    """Classify the user input.

    ``mode`` lets the user override auto-detection ("species", "compound",
    "spectrum"). For compounds we still detect the identifier type.
    """
    text = (value or "").strip()

    if filename or mode == "spectrum":
        ext = (filename or "").lower().rsplit(".", 1)[-1] if filename else ""
        id_type = ext if f".{ext}" in SPECTRUM_EXTENSIONS else "mgf"
        return Route("spectrum", id_type, filename or text, f"spectrum file uploaded ({filename or 'n/a'})")

    if not text:
        raise ValueError("Please enter a plant species, a compound, or upload a spectrum file.")

    if mode == "species":
        return Route("species", "binomial" if BINOMIAL_RE.match(text) else "common_name", text,
                     "user selected species mode")

    if INCHIKEY_RE.match(text):
        return Route("compound", "inchikey", text, "matches the 27-character InChIKey pattern")
    if text.startswith("InChI="):
        return Route("compound", "inchi", text, "starts with 'InChI='")
    if _looks_like_smiles(text):
        return Route("compound", "smiles", text, "parses as a valid SMILES string (RDKit)")
    if mode == "compound":
        return Route("compound", "name", text, "user selected compound mode; treated as a name")
    if BINOMIAL_RE.match(text):
        return Route("species", "binomial", text, "matches the 'Genus species' binomial pattern")
    return Route("compound", "name", text,
                 "not a binomial, SMILES, InChI or InChIKey; treated as a compound name "
                 "(switch to 'Plant species' mode for common plant names like 'peppermint')")
