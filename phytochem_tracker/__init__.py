"""PhytoChem-Tracker: an agentic dashboard that unifies phytochemical data.

Given a plant species, a compound (name / SMILES / InChI / InChIKey) or an
MS/MS spectrum file, the orchestrator dispatches specialised agents that query
open scientific databases and assemble a 7-tile report:

1. Top-10 phytochemicals (LOTUS occurrences ranked by PubMed paper counts)
2. Physico-chemical properties (PubChem computed + experimental, RDKit)
3. Spectra (PubChem spectral annotations, MoNA reference spectra, links)
4. Uses & applications (PubChem annotations, MeSH classes, links)
5. Chemical structures (2D/3D depictions, SDF/MOL downloads)
6. Biosynthetic pathways (KEGG, with Plant Reactome / PlantCyc links)
7. Gene IDs & enzymes (KEGG EC -> species genes, UniProt)

Every number shown comes from a database response, never from an LLM.
"""

__version__ = "3.0.0"
