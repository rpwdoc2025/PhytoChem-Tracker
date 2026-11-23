"""Data model, following the MVP specification (section 3).

Plain dataclasses keep the objects easy to serialise to JSON for the audit
log and exports (``dataclasses.asdict``).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass
class Provenance:
    """Where a piece of data came from. Attached to every tile."""

    source: str            # e.g. "PubChem PUG-View"
    query: str             # the URL or query string that produced the data
    retrieved_at: str      # ISO timestamp
    n: int | None = None   # count returned, when meaningful
    from_cache: bool = False
    note: str = ""


@dataclass
class PlantTaxon:
    name: str                       # accepted binomial, e.g. "Camellia sinensis"
    qid: str = ""                   # Wikidata item
    ncbi_taxid: str = ""
    common_names: list[str] = field(default_factory=list)
    kegg_org: str = ""              # KEGG organism code, e.g. "ath"
    resolved_from: str = ""         # what the user typed, if different


@dataclass
class Compound:
    name: str
    inchikey: str = ""
    cid: int | None = None
    smiles: str = ""
    inchi: str = ""
    formula: str = ""
    mw: float | None = None
    monoisotopic_mass: float | None = None
    synonyms: list[str] = field(default_factory=list)
    qid: str = ""                   # Wikidata item
    xrefs: dict[str, str] = field(default_factory=dict)  # kegg, chebi, knapsack, ...

    @property
    def key(self) -> str:
        """Stable identifier used to key per-compound results."""
        return self.inchikey or (f"CID{self.cid}" if self.cid else self.name)


@dataclass
class RankedCompound:
    rank: int
    compound: Compound
    paper_count: int
    lotus_references: int = 0
    pubmed_query: str = ""
    sources: list[str] = field(default_factory=list)


@dataclass
class PropertyValue:
    name: str
    value: str
    unit: str = ""
    source: str = ""
    reference: str = ""


@dataclass
class SpectrumRef:
    type: str            # MS, MS/MS, GC-MS, 1H NMR, 13C NMR, IR, UV-Vis, ...
    source: str
    url: str = ""
    accession: str = ""
    description: str = ""
    peaks: list[tuple[float, float]] = field(default_factory=list)


@dataclass
class UseRecord:
    category: str        # therapeutic / agricultural-veterinary / industrial-commercial / pharmacological class / other
    text: str
    source: str
    url: str = ""


@dataclass
class PathwayRef:
    name: str
    db: str
    pathway_id: str = ""
    url: str = ""
    enzymes: list[str] = field(default_factory=list)  # EC numbers


@dataclass
class GeneMap:
    ec: str
    enzyme_name: str
    species: str
    gene_id: str
    db: str
    url: str = ""
    protein: str = ""
    evidence: str = ""


@dataclass
class CompoundDetails:
    """Everything tiles 2-7 show for one compound."""

    compound: Compound
    computed: dict[str, Any] = field(default_factory=dict)
    experimental: list[PropertyValue] = field(default_factory=list)
    solid_state: str = ""
    spectra: list[SpectrumRef] = field(default_factory=list)
    uses: list[UseRecord] = field(default_factory=list)
    structure: dict[str, Any] = field(default_factory=dict)
    pathways: list[PathwayRef] = field(default_factory=list)
    enzymes_ec: list[str] = field(default_factory=list)
    genes: list[GeneMap] = field(default_factory=list)
    genes_status: str = "not run"
    links: dict[str, str] = field(default_factory=dict)
    provenance: dict[str, list[Provenance]] = field(default_factory=dict)
    errors: dict[str, str] = field(default_factory=dict)

    def add_prov(self, tile: str, prov: Provenance) -> None:
        self.provenance.setdefault(tile, []).append(prov)


@dataclass
class Report:
    run_id: str
    input_type: str
    input_value: str
    taxon: PlantTaxon | None = None
    top: list[RankedCompound] = field(default_factory=list)
    occurrences: list[dict[str, Any]] = field(default_factory=list)  # compound mode: taxa it occurs in
    spectrum_matches: list[dict[str, Any]] = field(default_factory=list)
    spectrum_summary: dict[str, Any] = field(default_factory=dict)
    details: dict[str, CompoundDetails] = field(default_factory=dict)
    summary: str = ""
    summary_model: str = ""
    provenance: dict[str, list[Provenance]] = field(default_factory=dict)
    errors: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    elapsed_s: float = 0.0

    def add_prov(self, tile: str, prov: Provenance) -> None:
        self.provenance.setdefault(tile, []).append(prov)

    def compounds(self) -> list[Compound]:
        return [d.compound for d in self.details.values()]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)
