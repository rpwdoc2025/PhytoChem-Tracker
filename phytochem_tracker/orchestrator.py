"""Orchestrator: routes the input and coordinates the agents.

    user input ──► Router ──┬─ species  ─► taxon resolution ─► Literature agent (Top-10)
                            ├─ compound ─► PubChem resolution ─► occurrences (LOTUS)
                            └─ spectrum ─► parse ─► precursor dereplication ─► MoNA cosine
                                              │
                    per compound, in parallel ▼
          PhysChem · Spectral · Uses · Structure · Pathway ─► Gene agent ─► (optional) LLM summary

The control flow is deterministic and fully logged; an LLM is only used, if
configured, to (a) map a common plant name to a binomial when Wikidata cannot,
and (b) write a summary that is restricted to the retrieved data.
"""

from __future__ import annotations

import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Callable

from . import agents, chem, config, llm
from .audit import AuditLog
from .config import Settings
from .connectors import kegg, links, pubchem, wikidata
from .http import Client
from .models import Compound, CompoundDetails, PlantTaxon, RankedCompound, Report
from .router import Route
from .spectra import Spectrum, cosine_similarity, match_precursor, parse_spectra

ProgressCB = Callable[[str, str], None]


def _noop(agent: str, msg: str) -> None:  # pragma: no cover
    pass


class Orchestrator:
    def __init__(self, settings: Settings | None = None, audit: AuditLog | None = None,
                 progress: ProgressCB | None = None, max_workers: int = 6):
        self.settings = settings or Settings()
        self.audit = audit or AuditLog()
        limits = dict(config.RATE_LIMITS)
        if self.settings.ncbi_api_key:
            limits["eutils.ncbi.nlm.nih.gov"] = 9.0
        session = None
        if self.settings.offline_demo:
            from .demo import DemoSession

            session, limits = DemoSession(), {"demo.invalid": 1000.0}
        self.client = Client(self.audit, use_cache=self.settings.use_cache,
                             cache_hours=self.settings.cache_hours, rate_limits=limits, session=session)
        self._progress = progress or _noop
        self.max_workers = max_workers

    # ------------------------------------------------------------------
    def progress(self, agent: str, msg: str) -> None:
        self.audit.log("agent", agent=agent, message=msg)
        self._progress(agent, msg)

    # ------------------------------------------------------------------
    def run(self, route: Route, file_bytes: bytes | None = None, context_species: str = "",
            spectrum_index: int | None = None) -> Report:
        t0 = time.perf_counter()
        report = Report(run_id=self.audit.run_id, input_type=route.input_type, input_value=route.value)
        self.audit.log("route", input_type=route.input_type, id_type=route.id_type, value=route.value,
                       reason=route.reason)
        self.progress("Router agent", f"Input classified as {route.input_type} ({route.id_type}): {route.reason}.")
        try:
            if route.input_type == "species":
                self._run_species(report, route.value)
            elif route.input_type == "compound":
                self._run_compound(report, route, context_species)
            else:
                self._run_spectrum(report, route.value, file_bytes or b"", context_species, spectrum_index)
        except Exception as exc:  # noqa: BLE001 - surface any failure in the report, not as a crash
            report.errors.append(f"{type(exc).__name__}: {exc}")
            self.audit.log("error", error=str(exc))
        if report.details and self.settings.llm_provider not in ("", "none") and self.settings.llm_api_key:
            self._summarise(report)
        report.elapsed_s = round(time.perf_counter() - t0, 2)
        self.audit.log("done", elapsed_s=report.elapsed_s, **self.audit.stats())
        return report

    # ------------------------------------------------------------------
    # Taxon resolution
    # ------------------------------------------------------------------
    def resolve_taxon(self, name: str) -> PlantTaxon | None:
        self.progress("Taxonomy agent", f"Resolving '{name}' on Wikidata…")
        t = wikidata.resolve_taxon(self.client, name)
        resolved_from = ""
        if not t and self.settings.llm_provider not in ("", "none") and self.settings.llm_api_key:
            guess = llm.binomial_for_common_name(self.settings, name, self.audit)
            if guess:
                self.progress("Taxonomy agent", f"LLM suggested '{guess}' for '{name}'; verifying on Wikidata…")
                t = wikidata.resolve_taxon(self.client, guess)
                resolved_from = name
        if not t:
            return None
        taxon = PlantTaxon(name=t["name"], qid=t["qid"], ncbi_taxid=t["ncbi_taxid"],
                           common_names=t["common_names"], resolved_from=resolved_from or (name if name != t["name"] else ""))
        taxon.kegg_org = kegg.organism_code(self.client, taxon.name)
        self.progress("Taxonomy agent", f"{taxon.name} → Wikidata {taxon.qid}, NCBI:txid{taxon.ncbi_taxid or '?'}, "
                                        f"KEGG organism: {taxon.kegg_org or 'none'}.")
        return taxon

    # ------------------------------------------------------------------
    # Mode 1: plant species
    # ------------------------------------------------------------------
    def _run_species(self, report: Report, name: str) -> None:
        taxon = self.resolve_taxon(name)
        if not taxon:
            report.errors.append(f"Could not find a taxon named '{name}' on Wikidata. Check the spelling, use the "
                                 "binomial (e.g. 'Mentha × piperita'), or enable an LLM to resolve common names.")
            return
        report.taxon = taxon
        lit = agents.LiteratureAgent(self.client, self.settings, self.progress)
        rows = lit.candidates(taxon)
        report.add_prov("top10", self.client.provenance(wikidata.SOURCE, f"compounds with P703 (found in taxon) {taxon.qid}",
                                                        n=len(rows)))
        if not rows:
            report.notes.append(f"LOTUS has no structure–organism pairs for {taxon.name} yet.")
            return
        top, props = lit.rank(taxon, rows)
        report.top = top
        for rc in top:
            report.add_prov("top10", self.client.provenance("PubMed E-utilities", rc.pubmed_query, n=rc.paper_count))
        report.notes.append(
            f"Ranked {lit.n_ranked} of {len(rows)} LOTUS compounds "
            f"(pre-selected by number of LOTUS references) by PubMed title/abstract co-mentions with '{taxon.name}'.")
        compounds = [rc.compound for rc in top]
        self._enrich(report, compounds[: self.settings.n_top], props, taxon)

    # ------------------------------------------------------------------
    # Mode 2: compound
    # ------------------------------------------------------------------
    def resolve_compound(self, value: str, id_type: str) -> tuple[Compound | None, dict[int, dict]]:
        self.progress("Resolver agent", f"Resolving {id_type} '{value[:60]}' on PubChem…")
        cids = pubchem.cids_for(self.client, value, id_type)
        if not cids and id_type == "name" and chem.mol_from_smiles(value) is not None:
            cids = pubchem.cids_for(self.client, value, "smiles")
        if not cids:
            if id_type == "smiles" and chem.mol_from_smiles(value) is not None:
                # Not in PubChem: still analyse the structure locally.
                return Compound(name=value, smiles=value, inchikey=chem.inchikey_from_smiles(value)), {}
            return None, {}
        props = pubchem.properties(self.client, cids[:1])
        rec = props.get(cids[0], {})
        comp = pubchem.compound_from_props(rec, fallback_name=value)
        comp.synonyms = pubchem.literature_synonyms(pubchem.synonyms(self.client, cids[0]), 5)
        return comp, props

    def _run_compound(self, report: Report, route: Route, context_species: str) -> None:
        comp, props = self.resolve_compound(route.value, route.id_type)
        if not comp:
            report.errors.append(f"PubChem has no compound matching '{route.value}'.")
            return
        taxon = self.resolve_taxon(context_species) if context_species else None
        report.taxon = taxon
        if comp.inchikey:
            self.progress("Literature agent", "Looking up organisms reported to contain this compound (LOTUS)…")
            try:
                report.occurrences = wikidata.taxa_for_compound(self.client, comp.inchikey)
                report.add_prov("top10", self.client.provenance(wikidata.SOURCE, f"P703 for InChIKey {comp.inchikey}",
                                                                n=len(report.occurrences)))
            except Exception as exc:  # noqa: BLE001
                report.errors.append(f"LOTUS occurrence lookup failed: {exc}")
        report.top = [RankedCompound(rank=1, compound=comp, paper_count=-1, sources=["PubChem"])]
        if not taxon:
            report.notes.append("No species context given: pathways use KEGG reference maps and genes are not mapped. "
                                "Enter a plant species in 'Species context' to map genes.")
        self._enrich(report, [comp], props, taxon)

    # ------------------------------------------------------------------
    # Mode 3: MS/MS spectrum
    # ------------------------------------------------------------------
    def _run_spectrum(self, report: Report, filename: str, data: bytes, context_species: str,
                      spectrum_index: int | None) -> None:
        spectra = parse_spectra(filename, data)
        self.progress("Spectral agent", f"Parsed {len(spectra)} MS/MS spectra from {filename}.")
        report.spectrum_summary = {
            "file": filename, "n_spectra": len(spectra),
            "spectra": [{"index": i, "title": s.title, "precursor_mz": s.precursor_mz, "polarity": s.polarity,
                         "n_peaks": len(s.peaks)} for i, s in enumerate(spectra[:500])],
        }
        if not spectra:
            report.errors.append("No spectra could be read from the file (supported: MGF, MSP, mzML).")
            return
        if not context_species:
            report.errors.append("Dereplication needs a plant species: enter one in 'Species context' so the precursor "
                                 "masses can be matched against compounds known from that plant (LOTUS).")
            return
        taxon = self.resolve_taxon(context_species)
        if not taxon:
            report.errors.append(f"Could not resolve species '{context_species}'.")
            return
        report.taxon = taxon
        rows = agents.LiteratureAgent(self.client, self.settings, self.progress).candidates(taxon)
        cands = []
        for r in rows:
            m = chem.monoisotopic_mass(r["smiles"]) if r["smiles"] else None
            if m:
                cands.append({"name": r["name"], "inchikey": r["inchikey"], "smiles": r["smiles"], "cid": r["cid"],
                              "qid": r["qid"], "monoisotopic_mass": m, "references": r["references"]})
        self.progress("Spectral agent", f"Matching precursors against {len(cands)} {taxon.name} compounds "
                                        f"(±{self.settings.ppm_tolerance} ppm)…")
        selected = range(len(spectra)) if spectrum_index is None else [spectrum_index]
        matches: list[dict[str, Any]] = []
        for i in selected:
            s = spectra[i]
            if s.precursor_mz is None:
                continue
            for h in match_precursor(s.precursor_mz, cands, s.polarity, self.settings.ppm_tolerance)[:10]:
                matches.append({"spectrum_index": i, "spectrum": s.title, "precursor_mz": s.precursor_mz, **h})
        report.add_prov("spectrum", self.client.provenance(wikidata.SOURCE, f"LOTUS candidates for {taxon.qid}",
                                                           n=len(cands), note="precursor m/z matching with RDKit masses"))
        # MS/MS similarity against MoNA references for the best candidates of the chosen spectrum
        focus = spectra[spectrum_index] if spectrum_index is not None else spectra[0]
        self._score_against_mona(focus, [m for m in matches if m["spectrum"] == focus.title][:5])
        matches.sort(key=lambda m: (m["spectrum_index"], -(m.get("cosine") or 0), abs(m["ppm_error"])))
        report.spectrum_matches = matches
        # Enrich the distinct matched compounds so tiles 2-7 work
        seen, comps = set(), []
        for m in matches:
            if m["inchikey"] in seen:
                continue
            seen.add(m["inchikey"])
            comps.append(Compound(name=m["name"], inchikey=m["inchikey"], smiles=m["smiles"], cid=m["cid"], qid=m["qid"]))
            if len(comps) >= self.settings.n_top:
                break
        if not comps:
            report.notes.append("No precursor matched a known compound of this species within the tolerance.")
            return
        cids = [c.cid for c in comps if c.cid]
        props = pubchem.properties(self.client, cids) if cids else {}
        for i, c in enumerate(comps):
            if c.cid in props:
                full = pubchem.compound_from_props(props[c.cid], fallback_name=c.name)
                full.qid = c.qid
                comps[i] = full
        report.top = [RankedCompound(rank=i + 1, compound=c, paper_count=-1, sources=["precursor match"])
                      for i, c in enumerate(comps)]
        self._enrich(report, comps, props, taxon)

    def _score_against_mona(self, spectrum: Spectrum, matches: list[dict[str, Any]]) -> None:
        from .connectors import mona

        for m in matches:
            try:
                refs = mona.spectra_for_inchikey(self.client, m["inchikey"], size=20)
            except Exception as exc:  # noqa: BLE001
                m["cosine_note"] = f"MoNA unavailable: {exc}"[:120]
                continue
            best = (0.0, 0, "")
            for r in refs:
                if r["ms_level"] and r["ms_level"] != "MS2":
                    continue
                try:
                    prec = float(r["precursor_mz"]) if r["precursor_mz"] else None
                except ValueError:
                    prec = None
                score, n = cosine_similarity(spectrum.peaks, r["peaks"], self.settings.fragment_tolerance,
                                             spectrum.precursor_mz, prec)
                if score > best[0]:
                    best = (score, n, r["id"])
            m["cosine"], m["matched_peaks"], m["mona_id"] = best
            m["mona_references"] = len(refs)

    # ------------------------------------------------------------------
    # Per-compound enrichment (tiles 2-7), in parallel
    # ------------------------------------------------------------------
    def _enrich(self, report: Report, compounds: list[Compound], props: dict[int, dict], taxon: PlantTaxon | None) -> None:
        if not compounds:
            return
        species = taxon.name if taxon else ""
        try:
            xrefs = wikidata.compound_xrefs(self.client, [c.inchikey for c in compounds])
        except Exception as exc:  # noqa: BLE001
            xrefs = {}
            report.errors.append(f"Wikidata cross-reference lookup failed: {exc}")
        for c in compounds:
            c.xrefs.update(xrefs.get(c.inchikey, {}))
            if not c.cid and c.xrefs.get("pubchem_cid", "").isdigit():
                c.cid = int(c.xrefs["pubchem_cid"])
            d = CompoundDetails(compound=c)
            d.links.update({f"compound::{k}": v for k, v in links.compound_links(c, species).items()})
            report.details[c.key] = d

        n_detail = len(compounds)
        self.progress("Orchestrator", f"Dispatching PhysChem, Spectral, Uses, Structure and Pathway agents "
                                      f"for {n_detail} compound(s) in parallel…")

        def work(c: Compound) -> str:
            d = report.details[c.key]
            agents._safe(d, "physchem", lambda: agents.physchem_agent(self.client, d, props.get(c.cid)))
            agents._safe(d, "structure", lambda: agents.structure_agent(self.client, d))
            agents._safe(d, "spectra", lambda: agents.spectral_agent(self.client, d))
            agents._safe(d, "uses", lambda: agents.uses_agent(self.client, d))
            agents._safe(d, "pathways", lambda: agents.pathway_agent(self.client, d, taxon))
            return c.name

        with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
            futures = {pool.submit(work, c): c for c in compounds}
            for k, fut in enumerate(as_completed(futures), 1):
                self.progress("Orchestrator", f"[{k}/{n_detail}] finished {fut.result()}.")

        if taxon:
            gene_targets = compounds[: max(1, self.settings.gene_compounds)]
            self.progress("Gene agent", f"Mapping EC numbers to {taxon.name} genes for "
                                        f"{len(gene_targets)} compound(s) (others on demand)…")
            with ThreadPoolExecutor(max_workers=3) as pool:
                list(pool.map(lambda c: self.map_genes(report.details[c.key], taxon), gene_targets))
            for c in compounds[len(gene_targets):]:
                report.details[c.key].genes_status = "not run yet — use 'Map genes' in tile 7"

    def map_genes(self, d: CompoundDetails, taxon: PlantTaxon) -> None:
        d.genes = []
        agents._safe(d, "genes", lambda: agents.gene_agent(self.client, d, taxon, self.settings))

    # ------------------------------------------------------------------
    def _summarise(self, report: Report) -> None:
        self.progress("Synthesis agent", f"Writing a grounded summary with {self.settings.llm_provider}…")
        try:
            report.summary, report.summary_model = llm.grounded_summary(self.settings, report, self.audit)
        except Exception as exc:  # noqa: BLE001
            report.errors.append(f"LLM summary failed: {exc}")
