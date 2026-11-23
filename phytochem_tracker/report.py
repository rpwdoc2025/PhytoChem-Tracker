"""Exports: JSON (full, with audit), CSV tables (zipped) and a PDF summary."""

from __future__ import annotations

import io
import json
import unicodedata
import zipfile
from datetime import datetime
from typing import Any

import pandas as pd

from .audit import AuditLog
from .models import Report


# --------------------------------------------------------------------------
# Tables (shared by the UI and the exports)
# --------------------------------------------------------------------------
def top_table(report: Report) -> pd.DataFrame:
    rows = []
    for rc in report.top:
        c = rc.compound
        rows.append({
            "Rank": rc.rank, "Compound": c.name,
            "PubMed papers": rc.paper_count if rc.paper_count >= 0 else None,
            "LOTUS refs": rc.lotus_references or None, "Formula": c.formula,
            "MW (g/mol)": round(c.mw, 2) if c.mw else None, "PubChem CID": c.cid, "InChIKey": c.inchikey,
            "PubMed query": rc.pubmed_query,
        })
    return pd.DataFrame(rows)


def properties_table(report: Report) -> pd.DataFrame:
    rows = []
    for d in report.details.values():
        base = {"Compound": d.compound.name, "InChIKey": d.compound.inchikey}
        for group, props in d.computed.items():
            for k, v in props.items():
                rows.append({**base, "Property": k, "Value": v, "Source": group})
        for p in d.experimental:
            rows.append({**base, "Property": p.name, "Value": p.value, "Source": p.source or "PubChem"})
    return pd.DataFrame(rows)


def spectra_table(report: Report) -> pd.DataFrame:
    return pd.DataFrame([{"Compound": d.compound.name, "Type": s.type, "Source": s.source, "Accession": s.accession,
                          "Description": s.description, "URL": s.url, "Peaks": len(s.peaks)}
                         for d in report.details.values() for s in d.spectra])


def uses_table(report: Report) -> pd.DataFrame:
    return pd.DataFrame([{"Compound": d.compound.name, "Category": u.category, "Use": u.text, "Source": u.source,
                          "URL": u.url} for d in report.details.values() for u in d.uses])


def pathways_table(report: Report) -> pd.DataFrame:
    return pd.DataFrame([{"Compound": d.compound.name, "Pathway": p.name, "Database": p.db, "ID": p.pathway_id,
                          "EC numbers": ", ".join(p.enzymes), "URL": p.url}
                         for d in report.details.values() for p in d.pathways])


def genes_table(report: Report) -> pd.DataFrame:
    return pd.DataFrame([{"Compound": d.compound.name, "EC": g.ec, "Enzyme": g.enzyme_name, "Species": g.species,
                          "Gene / ID": g.gene_id, "Database": g.db, "Protein": g.protein, "Evidence": g.evidence,
                          "URL": g.url} for d in report.details.values() for g in d.genes])


def spectrum_matches_table(report: Report) -> pd.DataFrame:
    cols = ["spectrum", "precursor_mz", "name", "adduct", "theoretical_mz", "ppm_error", "cosine", "matched_peaks",
            "mona_id", "references", "inchikey"]
    return pd.DataFrame([{k: m.get(k) for k in cols} for m in report.spectrum_matches])


# --------------------------------------------------------------------------
# Exports
# --------------------------------------------------------------------------
def to_json(report: Report, audit: AuditLog | None = None) -> bytes:
    d = report.to_dict()
    for det in d.get("details", {}).values():  # drop bulky structure files from JSON
        s = det.get("structure", {})
        for k in ("svg", "sdf3d", "mol2d", "sdf2d"):
            s.pop(k, None)
    if audit:
        d["audit"] = audit.events
    return json.dumps(d, indent=2, default=str).encode()


def to_csv_zip(report: Report) -> bytes:
    buf = io.BytesIO()
    tables = {"top_compounds.csv": top_table(report), "properties.csv": properties_table(report),
              "spectra.csv": spectra_table(report), "uses.csv": uses_table(report),
              "pathways.csv": pathways_table(report), "genes.csv": genes_table(report)}
    if report.spectrum_matches:
        tables["spectrum_matches.csv"] = spectrum_matches_table(report)
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, df in tables.items():
            zf.writestr(name, df.to_csv(index=False))
        sdfs = "".join(det.structure.get("sdf2d", "") for det in report.details.values())
        if sdfs:
            zf.writestr("structures_2d.sdf", sdfs)
    return buf.getvalue()


_REPL = {"α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta", "°": " deg", "Å": "A", "²": "2", "³": "3",
         "×": "x", "–": "-", "—": "-", "±": "+/-", "µ": "u", "μ": "u", "’": "'", "“": '"', "”": '"', "…": "..."}


def _latin1(text: Any) -> str:
    s = str(text if text is not None else "")
    for k, v in _REPL.items():
        s = s.replace(k, v)
    s = unicodedata.normalize("NFKD", s)
    return s.encode("latin-1", "replace").decode("latin-1")


def to_pdf(report: Report) -> bytes:
    from fpdf import FPDF

    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()
    pdf.set_fill_color(22, 101, 52)
    pdf.rect(0, 0, 210, 18, "F")
    pdf.set_text_color(255, 255, 255)
    pdf.set_font("Helvetica", "B", 14)
    pdf.set_xy(10, 5)
    pdf.cell(0, 8, "PhytoChem-Tracker report")
    pdf.set_text_color(0, 0, 0)
    pdf.set_xy(10, 24)
    pdf.set_font("Helvetica", "", 10)
    title = report.taxon.name if (report.taxon and report.input_type == "species") else report.input_value
    pdf.multi_cell(0, 5, _latin1(f"Query: {title} ({report.input_type})\nRun ID: {report.run_id}\n"
                                 f"Generated: {datetime.now():%Y-%m-%d %H:%M}"))
    if report.summary:
        pdf.ln(2)
        pdf.set_font("Helvetica", "B", 11)
        pdf.cell(0, 7, _latin1(f"AI summary of retrieved data ({report.summary_model})"), new_x="LMARGIN", new_y="NEXT")
        pdf.set_font("Helvetica", "", 9)
        pdf.multi_cell(0, 4.5, _latin1(report.summary))

    def section(title: str) -> None:
        pdf.ln(3)
        pdf.set_font("Helvetica", "B", 11)
        pdf.cell(0, 7, _latin1(title), new_x="LMARGIN", new_y="NEXT")
        pdf.set_font("Helvetica", "", 8.5)

    section("1. Top phytochemicals" if report.input_type == "species" else "1. Compound(s)")
    widths = [10, 62, 22, 30, 22, 44]
    header = ["#", "Compound", "Papers", "Formula", "MW", "InChIKey"]
    pdf.set_font("Helvetica", "B", 8.5)
    for w, h in zip(widths, header):
        pdf.cell(w, 6, h, border=1)
    pdf.ln()
    pdf.set_font("Helvetica", "", 8)
    for rc in report.top:
        c = rc.compound
        vals = [rc.rank, c.name[:40], rc.paper_count if rc.paper_count >= 0 else "-", c.formula,
                f"{c.mw:.2f}" if c.mw else "", c.inchikey[:27]]
        for w, v in zip(widths, vals):
            pdf.cell(w, 5.5, _latin1(v), border=1)
        pdf.ln()

    for d in report.details.values():
        c = d.compound
        section(f"{c.name}  (CID {c.cid or '-'})")
        lines = []
        pub = d.computed.get("PubChem (computed)", {})
        for k in ("Molecular formula", "Molecular weight (g/mol)", "XLogP3 (lipophilicity)", "TPSA (Å²)"):
            if k in pub:
                lines.append(f"{k}: {pub[k]}")
        for p in d.experimental:
            if p.name in ("Melting Point", "Boiling Point", "Solubility", "Physical Description"):
                lines.append(f"{p.name}: {p.value[:150]} [{p.source}]")
        if d.solid_state:
            lines.append(f"Solid state: {d.solid_state}")
        cats = sorted({u.category for u in d.uses})
        if cats:
            lines.append("Use categories: " + "; ".join(cats))
        if d.spectra:
            lines.append("Spectra available: " + ", ".join(sorted({s.type for s in d.spectra})))
        if d.pathways:
            lines.append("Pathways (KEGG): " + "; ".join(p.name for p in d.pathways[:6]))
        if d.genes:
            lines.append(f"Genes/proteins mapped: {len(d.genes)} (" +
                         ", ".join(sorted({g.gene_id for g in d.genes})[:8]) + ")")
        pdf.multi_cell(0, 4.3, _latin1("\n".join(lines) or "No data retrieved."))

    section("Sources and reproducibility")
    pdf.multi_cell(0, 4.3, _latin1(
        "Data were retrieved live from LOTUS/Wikidata, PubChem (PUG-REST, PUG-View), PubMed E-utilities, KEGG, "
        "UniProt and MoNA. Every query, with timestamp and cache status, is recorded in the JSON export "
        f"(audit section) for run {report.run_id}. " + " ".join(report.notes)))
    return bytes(pdf.output())
