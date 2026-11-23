# Architecture and scientific methods

This document explains *how* each tile is produced, so that results can be judged, reproduced
and, where needed, challenged. It maps directly onto the MVP specification
(scope → architecture → data model → ranking logic → API plan → non-functional requirements).

## 1. Request flow

1. **Router agent** (`router.py`) classifies the input with explicit rules and records its
   reason: InChIKey pattern → `InChI=` prefix → valid SMILES (RDKit parse) → binomial pattern
   → otherwise a compound name. The user can override the mode.
2. **Entity resolution** turns what the user typed into stable identifiers:
   * species → Wikidata taxon (`P225` scientific name, or `P1843` English common name),
     NCBI taxonomy ID (`P685`) and KEGG organism code;
   * compound → PubChem CID → InChIKey, then Wikidata cross-references (KEGG `P665`,
     ChEBI `P683`, KNApSAcK `P2064`).
   InChIKey is the join key across all sources.
3. **Per-compound agents** run in a thread pool. Each agent writes its own tile and
   provenance records, and stores any failure in `details.errors[tile]`.
4. **Gene mapping** runs eagerly for the top *N* compounds (default 3) and on demand for the rest,
   because it is the most call-intensive step.

## 2. Tile 1: ranking ("Top-10")

1. Collect every compound with a LOTUS *found in taxon* statement for the species, including
   varieties/subspecies (`wdt:P171*`), with the number of distinct supporting references.
2. Pre-select the `n_candidates` (default 30) best-documented compounds, i.e. those with the most
   LOTUS references. Entries with neither a name nor a PubChem CID are skipped (they cannot be
   searched in PubMed by name).
3. For each candidate, build a PubMed query from its name plus up to two literature synonyms from
   PubChem (registry numbers, database codes and long systematic names are filtered out):
   `("quercetin"[tiab] OR "sophoretin"[tiab]) AND "Camellia sinensis"[tiab]`.
4. Rank by paper count; ties are broken by the number of LOTUS references, then by name.
   Every query string and count is shown in the tile and stored in the audit log.

**Limitations.** The paper count is a proxy for research attention, not for abundance in the
plant. Title/abstract searching misses papers that mention the compound only in the full text,
and very common names can over-count. The candidate pre-selection means a compound with few
LOTUS references but many papers can be missed: raise *LOTUS candidates* in Settings for a more
exhaustive ranking.

## 3. Tiles 2–7

| Tile | Method |
|---|---|
| Physico-chemical | PubChem computed properties (batch PUG-REST call); experimental values from the PUG-View *Experimental Properties* section, deduplicated, up to 4 per property, each with its original source (HSDB, NTP, …). "Solid state" is a keyword call (crystal/needle/prism vs amorphous/powder) on *Physical Description* and *Color/Form*, labelled "(reported)". RDKit descriptors are computed locally as a cross-check. |
| Spectra | PubChem *Spectral Information* grouped by spectrum type with record counts and sources; up to 5 MoNA reference spectra by InChIKey, MS/MS first. |
| Uses | PubChem sections *Drug and Medication Information*, *Use and Manufacturing*, *Pharmacology and Biochemistry*, *Agrochemical Information*. Category = section heading (therapeutic, MeSH pharmacological class) or keyword rules (pesticide/veterinary → agricultural; food/cosmetic/dye → industrial). The rules are in `pubchem.categorise_use`. |
| Structures | RDKit 2D depiction; PubChem 3D conformer when available, otherwise RDKit ETKDGv3 embedding + MMFF optimisation (labelled as computed). |
| Pathways | KEGG `link/pathway` for the compound; global overview maps are listed last. If the species has a KEGG genome, the organism-specific map is linked. Enzymes = KEGG `link/enzyme` for the compound, intersected with each of the first three specific pathways. |
| Genes | For up to 8 EC numbers (pathway enzymes first): KEGG GENES for the organism (where KEGG has the genome) and UniProtKB entries with that EC in the species' NCBI taxon (Swiss-Prot first). Phytozome has no anonymous API, so gene IDs are given for searching there. |

## 4. Spectrum dereplication

* Candidate masses: RDKit monoisotopic masses from each LOTUS compound's SMILES.
* Adducts: [M+H]+, [M+Na]+, [M+K]+, [M+NH4]+, [M+H−H2O]+, [2M+H]+ (positive);
  [M−H]−, [M+Cl]−, [M+FA−H]−, [M−H2O−H]−, [2M−H]− (negative); polarity is taken from the file.
* Tolerance: ±10 ppm by default.
* Similarity: greedy cosine on square-root-scaled intensities, peaks ≥1 % of the base peak,
  precursor region removed, fragment tolerance 0.02 Da (the GNPS/matchms *CosineGreedy*
  convention).

A precursor match is only a candidate: isomers (e.g. quercetin vs. morin) have identical masses.
The cosine score and orthogonal evidence (retention time, standards) are needed for confident
annotation (MSI level 2 or 1).

## 5. Non-functional design

| Requirement (spec) | Implementation |
|---|---|
| Latency, caching | `requests-cache` SQLite backend, 72 h TTL, 404s cached, stale-if-error |
| Rate limits | per-host minimum interval (`config.RATE_LIMITS`); NCBI key raises PubMed to 9/s |
| Retries | urllib3 `Retry`, 3 attempts, exponential backoff, honours `Retry-After` |
| Auditability | `audit.AuditLog`: JSONL per run, API keys redacted, exported with the JSON report |
| Security | read-only calls only; keys come from secrets/env and are never logged |
| Hallucination guard | LLM never supplies values; its summary is labelled and generated only from the retrieved JSON |

## 6. Extending

* **New source:** add a module in `connectors/`, call it from an agent in `agents.py`, add a
  fixture route in `demo.py` and a test.
* **Endpoint moved:** change the URL in `config.py` (or set `PCT_<NAME>_BASE`).
* **LangGraph / CrewAI:** the agents are pure functions of `(client, details, …)`, so they can be
  wrapped as graph nodes without changes.
