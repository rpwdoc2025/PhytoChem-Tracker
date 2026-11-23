"""PhytoChem-Tracker v3 — Streamlit dashboard.

Run locally:   streamlit run app.py
"""

from __future__ import annotations

import html
import json
import os
from typing import Any

import pandas as pd
import plotly.graph_objects as go
import streamlit as st
import streamlit.components.v1 as components

from phytochem_tracker import __version__, report as rep
from phytochem_tracker.audit import AuditLog
from phytochem_tracker.config import Settings
from phytochem_tracker.demo import DEMO_MGF, DEMO_SPECIES
from phytochem_tracker.llm import KEY_ENV
from phytochem_tracker.models import CompoundDetails, Report
from phytochem_tracker.orchestrator import Orchestrator
from phytochem_tracker.router import route
from phytochem_tracker.spectra import parse_spectra

st.set_page_config(page_title="PhytoChem-Tracker", page_icon="🌿", layout="wide")

st.markdown("""
<style>
.block-container {padding-top: 1.6rem;}
.pct-hero {background: linear-gradient(120deg,#14532d 0%,#166534 45%,#4d7c0f 100%); color:#fff;
           padding: 1.1rem 1.4rem; border-radius: 14px; margin-bottom: 1rem;}
.pct-hero h1 {color:#fff; margin:0; font-size:1.9rem;}
.pct-hero p {margin:.25rem 0 0 0; opacity:.9;}
.pct-card {border:1px solid rgba(22,101,52,.25); border-radius:12px; padding:.8rem 1rem; height:100%;
           background: rgba(22,101,52,.04);}
.pct-card h4 {margin:.1rem 0 .3rem 0;}
.pct-badge {display:inline-block; padding:.08rem .5rem; border-radius:999px; font-size:.78rem;
            background:#dcfce7; color:#14532d; margin-right:.3rem; border:1px solid #86efac;}
.pct-demo {background:#fef9c3; border:1px solid #facc15; color:#713f12; padding:.5rem .8rem; border-radius:10px;}
</style>
""", unsafe_allow_html=True)


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def secret(name: str) -> str:
    try:
        return str(st.secrets.get(name, "")) or os.getenv(name, "")
    except Exception:  # noqa: BLE001 - no secrets file
        return os.getenv(name, "")


def provenance_box(items, label: str = "Provenance") -> None:
    if not items:
        return
    with st.expander(f"🔎 {label} ({len(items)})"):
        st.dataframe(pd.DataFrame([{
            "Source": p.source, "Query": p.query, "Results": p.n, "Cached": p.from_cache,
            "Retrieved (UTC)": p.retrieved_at, "Note": p.note} for p in items]),
            hide_index=True, width="stretch")


def link_list(d: CompoundDetails, prefix: str) -> None:
    items = [(k.split("::", 1)[1], v) for k, v in d.links.items() if k.startswith(prefix + "::")]
    if items:
        st.markdown(" · ".join(f"[{html.escape(k)}]({v.replace(' ', '%20').replace('(', '%28').replace(')', '%29')})"
                               for k, v in items))


def stem_plot(peaks, title: str = "", mirror=None, mirror_label: str = "reference") -> go.Figure:
    fig = go.Figure()

    def add(pk, sign, name, color):
        if not pk:
            return
        mx = max(i for _, i in pk) or 1
        xs, ys = [], []
        for mz, it in pk:
            xs += [mz, mz, None]
            ys += [0, sign * 100 * it / mx, None]
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", name=name, line=dict(color=color, width=1.6),
                                 hoverinfo="skip"))
        fig.add_trace(go.Scatter(x=[p[0] for p in pk], y=[sign * 100 * p[1] / mx for p in pk], mode="markers",
                                 marker=dict(size=4, color=color), name=name, showlegend=False,
                                 hovertemplate="m/z %{x:.4f}<br>rel. int. %{y:.0f}<extra></extra>"))

    add(peaks, 1, "query" if mirror else "spectrum", "#166534")
    if mirror:
        add(mirror, -1, mirror_label, "#b45309")
    fig.update_layout(title=title, height=320, margin=dict(l=10, r=10, t=40, b=10),
                      xaxis_title="m/z", yaxis_title="relative intensity (%)", template="plotly_white")
    return fig


def report_settings_ppm() -> float:
    s = st.session_state.get("settings")
    return s.ppm_tolerance if s else 10.0


def mol_viewer_3d(sdf: str, height: int = 360) -> None:
    payload = json.dumps(sdf)
    components.html(f"""
<div id="v" style="width:100%;height:{height}px;position:relative;border-radius:10px;overflow:hidden"></div>
<script src="https://cdn.jsdelivr.net/npm/3dmol@2.4.2/build/3Dmol-min.js"></script>
<script>
  const v = $3Dmol.createViewer(document.getElementById('v'), {{backgroundColor: 'white'}});
  v.addModel({payload}, 'sdf');
  v.setStyle({{}}, {{stick: {{radius: 0.16}}, sphere: {{scale: 0.23}}}});
  v.zoomTo(); v.render(); v.spin('y', 0.6);
</script>""", height=height + 10)


# --------------------------------------------------------------------------
# Shareable links: ?q=Camellia+sinensis&mode=species&context=...&demo=1&run=1
# --------------------------------------------------------------------------
_qp = st.query_params
_MODES = {"auto": "Auto-detect", "species": "Plant species", "compound": "Compound", "spectrum": "MS/MS spectrum"}
if "qp_loaded" not in st.session_state:
    st.session_state["qp_loaded"] = True
    if _qp.get("q"):
        st.session_state["query"] = _qp["q"]
    if _qp.get("context"):
        st.session_state["context"] = _qp["context"]
    if _qp.get("mode") in _MODES:
        st.session_state["mode"] = _MODES[_qp["mode"]]
    st.session_state["autorun"] = _qp.get("run") == "1"

# --------------------------------------------------------------------------
# Sidebar: inputs and settings
# --------------------------------------------------------------------------
with st.sidebar:
    st.markdown("### 🌿 Query")
    mode_label = st.radio("Input type", ["Auto-detect", "Plant species", "Compound", "MS/MS spectrum"],
                          horizontal=False, key="mode")
    mode = {"Auto-detect": "auto", "Plant species": "species", "Compound": "compound",
            "MS/MS spectrum": "spectrum"}[mode_label]

    placeholder = {"auto": "e.g. Camellia sinensis, quercetin, or an InChIKey",
                   "species": "e.g. Camellia sinensis or peppermint",
                   "compound": "name, SMILES, InChI or InChIKey", "spectrum": ""}[mode]
    uploaded, use_demo_mgf = None, False
    if mode != "spectrum":
        query = st.text_input("Search", key="query", placeholder=placeholder)

        def _use_example() -> None:
            if st.session_state.get("example"):
                st.session_state["query"] = st.session_state["example"]
                st.session_state["example"] = None

        st.pills("Examples", ["Camellia sinensis", "Curcuma longa", "Mentha × piperita", "quercetin", "caffeine"],
                 key="example", on_change=_use_example, label_visibility="collapsed")
    else:
        query = ""
        uploaded = st.file_uploader("MS/MS file", type=["mgf", "msp", "mzml"],
                                    help="MGF, MSP or mzML with MS2 scans and precursor m/z.")
        use_demo_mgf = st.checkbox("Use the bundled example MGF (tea extract)", value=uploaded is None)

    context_species = ""
    if mode in ("compound", "spectrum", "auto"):
        context_species = st.text_input(
            "Species context (optional for compounds, required for spectra)", key="context",
            placeholder="e.g. Camellia sinensis",
            help="Used to map enzymes to this plant's genes and, for spectra, as the candidate list "
                 "for dereplication.")

    with st.expander("⚙️ Settings"):
        n_top = st.slider("Compounds in ranked list", 3, 15, 10)
        n_cand = st.slider("LOTUS candidates sent to PubMed ranking", 10, 80, 30,
                           help="Higher = more thorough ranking, more API calls on the first run.")
        n_syn = st.slider("Extra synonyms per compound in PubMed query", 0, 4, 2)
        gene_n = st.slider("Compounds with eager gene mapping", 1, 10, 3)
        ppm = st.number_input("Precursor tolerance (ppm)", 1.0, 50.0, 10.0, 1.0)
        infra = st.checkbox("Include varieties/subspecies of the species", True)
        use_cache = st.checkbox("Use local cache (72 h)", True)
        offline = st.checkbox("Offline demo (fixture data, no network)",
                              value=secret("PCT_OFFLINE_DEMO") == "1" or _qp.get("demo") == "1",
                              help="Runs the whole pipeline on built-in sample responses. Values are illustrative.")
        ncbi_key = st.text_input("NCBI API key (optional)", value=secret("NCBI_API_KEY"), type="password",
                                 help="Raises the PubMed limit from 3 to 10 requests/s.")
        ncbi_email = st.text_input("Contact e-mail for NCBI (optional)", value=secret("NCBI_EMAIL"))

    with st.expander("🤖 AI assistant (optional)"):
        st.caption("Used only to map common plant names to binomials and to write a summary restricted to the "
                   "retrieved data. All values in the tiles come from databases.")
        provider = st.selectbox("Provider", ["none", "anthropic", "openai", "gemini"],
                                index=["none", "anthropic", "openai", "gemini"].index(
                                    secret("PCT_LLM_PROVIDER") or "none"),
                                format_func=lambda p: {"none": "Off", "anthropic": "Claude (Anthropic)",
                                                       "openai": "ChatGPT (OpenAI)", "gemini": "Gemini (Google)"}[p])
        llm_key = st.text_input("API key", value=secret(KEY_ENV.get(provider, "")) if provider != "none" else "",
                                type="password", disabled=provider == "none")
        llm_model = st.text_input("Model (blank = default)", value="", disabled=provider == "none")

    run = st.button("▶ Run agents", type="primary", width="stretch")
    st.caption(f"v{__version__} · [GitHub](https://github.com/rpwdoc2025/PhytoChem-Tracker)")

settings = Settings(
    n_top=n_top, n_candidates=n_cand, synonyms_per_compound=n_syn, gene_compounds=gene_n, ppm_tolerance=ppm,
    include_infraspecific=infra, ncbi_api_key=ncbi_key, ncbi_email=ncbi_email, llm_provider=provider,
    llm_api_key=llm_key, llm_model=llm_model, use_cache=use_cache, offline_demo=offline,
)

# --------------------------------------------------------------------------
# Header
# --------------------------------------------------------------------------
st.markdown("""<div class="pct-hero"><h1>🌿 PhytoChem-Tracker</h1>
<p>One query → seven tiles of phytochemical intelligence: literature-ranked compounds, properties, spectra, uses,
structures, biosynthetic pathways and genes, pulled live from open databases with a full audit trail.</p></div>""",
            unsafe_allow_html=True)

# --------------------------------------------------------------------------
# Run
# --------------------------------------------------------------------------
if st.session_state.pop("autorun", False):
    run = True

if run:
    st.query_params.from_dict({k: v for k, v in {
        "q": query, "mode": mode, "context": context_species, "demo": "1" if offline else ""}.items() if v})
    try:
        if mode == "spectrum":
            if uploaded is not None:
                fname, data = uploaded.name, uploaded.getvalue()
            elif use_demo_mgf:
                fname, data = "demo_tea_extract.mgf", DEMO_MGF.encode()
                context_species = context_species or DEMO_SPECIES
            else:
                st.error("Upload an MS/MS file or tick 'Use the bundled example MGF'.")
                st.stop()
            r = route(filename=fname, mode="spectrum")
        else:
            fname, data = None, None
            r = route(query, mode=mode)
    except ValueError as exc:
        st.error(str(exc))
        st.stop()

    audit = AuditLog()
    with st.status("Agents working…", expanded=True) as status:
        def progress(agent: str, msg: str) -> None:
            status.write(f"**{agent}** — {msg}")

        orch = Orchestrator(settings, audit, progress)
        result = orch.run(r, file_bytes=data, context_species=context_species.strip())
        status.update(label=f"Done in {result.elapsed_s:.1f} s", state="complete" if not result.errors else "error",
                      expanded=False)
    st.session_state.update(report=result, audit=audit, settings=settings, spectrum_file=(fname, data),
                            demo=offline)

report: Report | None = st.session_state.get("report")
audit: AuditLog | None = st.session_state.get("audit")

# --------------------------------------------------------------------------
# Landing page
# --------------------------------------------------------------------------
if report is None:
    c1, c2, c3 = st.columns(3)
    c1.markdown("""<div class="pct-card"><h4>🌱 Plant species</h4>Top-10 phytochemicals from <b>LOTUS</b>
occurrence data, ranked by <b>PubMed</b> papers that mention the compound together with the plant.</div>""",
                unsafe_allow_html=True)
    c2.markdown("""<div class="pct-card"><h4>⚗️ Compound</h4>Name, SMILES, InChI or InChIKey → properties, spectra,
uses, 2D/3D structure, KEGG pathways and genes, plus the plants it occurs in.</div>""", unsafe_allow_html=True)
    c3.markdown("""<div class="pct-card"><h4>📈 MS/MS spectrum</h4>Upload MGF/MSP/mzML: precursor m/z matched to the
plant's known compounds (adducts, ppm), then cosine-scored against <b>MoNA</b> reference spectra.</div>""",
                unsafe_allow_html=True)
    st.markdown("#### Data sources")
    st.markdown(" ".join(f'<span class="pct-badge">{s}</span>' for s in [
        "LOTUS / Wikidata", "PubChem PUG-REST & PUG-View", "PubMed E-utilities", "KEGG", "UniProt", "MoNA",
        "RDKit", "links: KNApSAcK · IMPPAT · Dr. Duke's · NP-MRD · GNPS · MassBank · PlantCyc · Plant Reactome · Phytozome"]),
        unsafe_allow_html=True)
    st.info("Pick an input type in the sidebar, enter a query (or click an example) and press **Run agents**. "
            "First runs take 30–90 s while the agents respect each database's rate limits; repeat queries are "
            "served from the local cache in seconds. Tick **Offline demo** to try the app without network access.")
    st.stop()

# --------------------------------------------------------------------------
# Results header
# --------------------------------------------------------------------------
if st.session_state.get("demo"):
    st.markdown('<div class="pct-demo">⚠️ <b>Offline demo:</b> results come from built-in sample responses. '
                'Identifiers and structures are real; paper counts, gene records and the reference spectrum are '
                'illustrative. Untick “Offline demo” for live data.</div>', unsafe_allow_html=True)

stats = audit.stats() if audit else {}
title = report.taxon.name if (report.taxon and report.input_type == "species") else (
    report.spectrum_summary.get("file") if report.input_type == "spectrum" else report.input_value)
st.subheader(f"Results for *{title}*")
m1, m2, m3, m4 = st.columns(4)
m1.metric("Compounds", len(report.details))
m2.metric("API calls", stats.get("http_calls", 0))
m3.metric("Served from cache", stats.get("from_cache", 0))
m4.metric("Run time", f"{report.elapsed_s:.1f} s")
if report.taxon:
    t = report.taxon
    extra = f" (resolved from '{t.resolved_from}')" if t.resolved_from else ""
    st.caption(f"Taxon: **{t.name}**{extra} · Wikidata [{t.qid}](https://www.wikidata.org/wiki/{t.qid}) · "
               f"NCBI txid {t.ncbi_taxid or '—'} · KEGG organism {t.kegg_org or 'not in KEGG'}"
               + (f" · common names: {', '.join(t.common_names[:4])}" if t.common_names else ""))
for e in report.errors:
    st.warning(e)
for n in report.notes:
    st.caption(f"ℹ️ {n}")
if report.summary:
    with st.container(border=True):
        st.markdown(f"**AI summary of the retrieved data** · `{report.summary_model}`")
        st.markdown(report.summary)

# --------------------------------------------------------------------------
# Spectrum results
# --------------------------------------------------------------------------
if report.input_type == "spectrum" and report.spectrum_summary:
    st.markdown("### 📈 Spectrum dereplication")
    fname, data = st.session_state.get("spectrum_file", (None, None))
    spectra = parse_spectra(fname, data) if fname and data else []
    sp_df = pd.DataFrame(report.spectrum_summary.get("spectra", []))
    if not sp_df.empty:
        idx = st.selectbox("Spectrum", sp_df["index"].tolist(),
                           format_func=lambda i: f"#{i + 1}: {sp_df.loc[i, 'title']} (m/z {sp_df.loc[i, 'precursor_mz']})")
        matches = [m for m in report.spectrum_matches if m["spectrum_index"] == idx]
        left, right = st.columns([3, 2])
        with right:
            if matches:
                cols = {"name": "Candidate", "adduct": "Adduct", "theoretical_mz": "Theor. m/z", "ppm_error": "Δ ppm",
                        "cosine": "Cosine", "matched_peaks": "Matched peaks", "mona_id": "MoNA reference"}
                mdf = pd.DataFrame(matches)
                mdf = mdf[[k for k in cols if k in mdf.columns]].rename(columns=cols)
                st.dataframe(mdf, hide_index=True, width="stretch")
                st.caption(f"{len(matches)} candidate(s) within ±{report_settings_ppm()} ppm.")
            else:
                st.info("No candidate within tolerance for this precursor.")
        with left:
            if idx < len(spectra):
                ref_peaks, ref_label = None, ""
                best = next((m for m in matches if m.get("mona_id")), None)
                if best:
                    for d in report.details.values():
                        for s in d.spectra:
                            if s.accession == best["mona_id"]:
                                ref_peaks, ref_label = s.peaks, f"MoNA {s.accession}"
                st.plotly_chart(stem_plot(spectra[idx].peaks, spectra[idx].title, ref_peaks, ref_label),
                                width="stretch")
        st.caption("Cosine = greedy fragment matching on square-root intensities (±0.02 Da), as in GNPS. "
                   "For network-level annotation, submit the file to [GNPS2](https://gnps2.org/).")

if not report.details:
    st.stop()

# --------------------------------------------------------------------------
# Compound selector + seven tiles
# --------------------------------------------------------------------------
keys = list(report.details)
label = {k: report.details[k].compound.name for k in keys}
sel = st.selectbox("Compound shown in tiles 2–7", keys, format_func=lambda k: label[k])
d = report.details[sel]
c = d.compound

tabs = st.tabs(["1 · Top phytochemicals", "2 · Physico-chemical", "3 · Spectra", "4 · Uses",
                "5 · Structures", "6 · Pathways", "7 · Genes & enzymes"])

# ---- Tile 1
with tabs[0]:
    if report.input_type == "species":
        df = rep.top_table(report)
        if not df.empty:
            a, b = st.columns([3, 2])
            with a:
                st.dataframe(df.drop(columns=["PubMed query"]), hide_index=True, width="stretch")
            with b:
                fig = go.Figure(go.Bar(x=df["PubMed papers"][::-1], y=df["Compound"][::-1], orientation="h",
                                       marker_color="#166534"))
                fig.update_layout(height=320, margin=dict(l=10, r=10, t=10, b=10), template="plotly_white",
                                  xaxis_title="PubMed papers")
                st.markdown("**PubMed papers co-mentioning the plant**")
                st.plotly_chart(fig, width="stretch")
            with st.expander("PubMed queries used for ranking"):
                st.dataframe(df[["Rank", "Compound", "PubMed papers", "PubMed query"]], hide_index=True,
                             width="stretch")
    elif report.input_type == "compound":
        st.markdown(f"**{c.name}** — occurs in the following organisms (LOTUS / Wikidata), most-referenced first:")
        if report.occurrences:
            st.dataframe(pd.DataFrame(report.occurrences), hide_index=True, width="stretch",
                         column_config={"taxon": "Organism", "wikidata": "Wikidata", "references": "LOTUS refs"})
        else:
            st.info("No occurrence records found in LOTUS for this compound.")
    else:
        if report.top:
            st.markdown(f"Compounds of *{report.taxon.name if report.taxon else 'the plant'}* matched to the "
                        "spectrum file (best match first):")
            st.dataframe(rep.top_table(report)[["Rank", "Compound", "Formula", "MW (g/mol)", "InChIKey"]],
                         hide_index=True, width="stretch")
    link_list(d, "compound")
    provenance_box(report.provenance.get("top10", []) + report.provenance.get("spectrum", []))

# ---- Tile 2
with tabs[1]:
    a, b = st.columns(2)
    for col, group in zip((a, b), ("PubChem (computed)", "RDKit (local)")):
        with col:
            st.markdown(f"**{group}**")
            props = d.computed.get(group, {})
            if props:
                st.dataframe(pd.DataFrame({"Property": list(props), "Value": [str(v) for v in props.values()]}),
                             hide_index=True, width="stretch")
            else:
                st.caption("Not available.")
    st.markdown("**Experimental properties (PubChem annotations, with original sources)**")
    if d.solid_state:
        st.markdown(f'<span class="pct-badge">Solid state: {html.escape(d.solid_state)}</span>',
                    unsafe_allow_html=True)
    if d.experimental:
        st.dataframe(pd.DataFrame([{"Property": p.name, "Value": p.value, "Source": p.source, "Reference": p.reference}
                                   for p in d.experimental]), hide_index=True, width="stretch",
                     column_config={"Reference": st.column_config.LinkColumn()})
    else:
        st.caption("PubChem has no experimental property annotations for this compound.")
    if "physchem" in d.errors:
        st.warning(d.errors["physchem"])
    provenance_box(d.provenance.get("physchem", []))

# ---- Tile 3
with tabs[2]:
    if d.spectra:
        st.dataframe(pd.DataFrame([{"Type": s.type, "Source": s.source, "Details": s.description, "Link": s.url}
                                   for s in d.spectra]), hide_index=True, width="stretch",
                     column_config={"Link": st.column_config.LinkColumn()})
        refs = [s for s in d.spectra if s.peaks]
        if refs:
            pick = st.selectbox("Reference spectrum", range(len(refs)),
                                format_func=lambda i: f"{refs[i].accession} · {refs[i].description}")
            st.plotly_chart(stem_plot(refs[pick].peaks, f"{c.name} — {refs[pick].accession}"),
                            width="stretch")
    else:
        st.caption("No spectral annotations retrieved.")
    st.markdown("**Search spectral databases:**")
    link_list(d, "spectra")
    for k in ("spectra", "spectra_mona"):
        if k in d.errors:
            st.caption(f"⚠️ {d.errors[k]}")
    provenance_box(d.provenance.get("spectra", []))

# ---- Tile 4
with tabs[3]:
    if d.uses:
        cats = sorted({u.category for u in d.uses})
        st.markdown(" ".join(f'<span class="pct-badge">{html.escape(x)}</span>' for x in cats),
                    unsafe_allow_html=True)
        for cat in cats:
            with st.expander(cat, expanded=True):
                for u in [u for u in d.uses if u.category == cat]:
                    src = f"[{u.source}]({u.url})" if u.url else u.source
                    st.markdown(f"- {u.text}  \n  <small>— {src}</small>", unsafe_allow_html=True)
        st.caption("Categories are assigned from PubChem section headings and keywords (see `pubchem.categorise_use`).")
    else:
        st.caption("PubChem has no use annotations for this compound.")
    st.markdown("**Ethnobotanical and medicinal-plant resources:**")
    link_list(d, "uses")
    provenance_box(d.provenance.get("uses", []))

# ---- Tile 5
with tabs[4]:
    s = d.structure
    a, b = st.columns(2)
    with a:
        st.markdown("**2D structure (RDKit)**")
        if s.get("svg"):
            st.markdown(f'<div style="background:white;border-radius:10px;padding:6px">{s["svg"]}</div>',
                        unsafe_allow_html=True)
        else:
            st.caption("No structure available.")
    with b:
        st.markdown(f"**3D structure** · {s.get('sdf3d_source', 'not available')}")
        if s.get("sdf3d"):
            mol_viewer_3d(s["sdf3d"])
    ids = {k: v for k, v in s.get("identifiers", {}).items() if v}
    if ids:
        st.dataframe(pd.DataFrame({"Identifier": list(ids), "Value": [str(v) for v in ids.values()]}),
                     hide_index=True, width="stretch")
    safe = "".join(ch if ch.isalnum() else "_" for ch in c.name)[:40] or "compound"
    dl = st.columns(4)
    if s.get("sdf2d"):
        dl[0].download_button("⬇ SDF (2D)", s["sdf2d"], f"{safe}_2d.sdf", "chemical/x-mdl-sdfile")
        dl[1].download_button("⬇ MOL (2D)", s.get("mol2d", ""), f"{safe}.mol", "chemical/x-mdl-molfile")
    if s.get("sdf3d"):
        dl[2].download_button("⬇ SDF (3D)", s["sdf3d"], f"{safe}_3d.sdf", "chemical/x-mdl-sdfile")
    if c.smiles:
        dl[3].download_button("⬇ SMILES", c.smiles, f"{safe}.smi", "text/plain")
    provenance_box(d.provenance.get("structure", []))

# ---- Tile 6
with tabs[5]:
    if d.pathways:
        st.dataframe(pd.DataFrame([{"Pathway": p.name, "Database": p.db, "ID": p.pathway_id,
                                    "Enzymes (EC) acting on this compound": ", ".join(p.enzymes), "Map": p.url}
                                   for p in d.pathways]), hide_index=True, width="stretch",
                     column_config={"Map": st.column_config.LinkColumn(display_text="open map")})
        if d.enzymes_ec:
            st.caption(f"KEGG lists {len(d.enzymes_ec)} enzyme(s) that use or produce this compound: "
                       + ", ".join(d.enzymes_ec[:20]))
    else:
        st.caption("No KEGG pathway linked to this compound.")
    st.markdown("**Plant pathway databases:**")
    link_list(d, "pathways")
    if "pathways" in d.errors:
        st.warning(d.errors["pathways"])
    provenance_box(d.provenance.get("pathways", []))

# ---- Tile 7
with tabs[6]:
    st.caption(f"Status: {d.genes_status}")
    if d.genes:
        st.dataframe(pd.DataFrame([{"EC": g.ec, "Enzyme": g.enzyme_name, "Gene / ID": g.gene_id, "Database": g.db,
                                    "Protein": g.protein, "Species": g.species, "Evidence": g.evidence, "Link": g.url}
                                   for g in d.genes]), hide_index=True, width="stretch",
                     column_config={"Link": st.column_config.LinkColumn()})
    gene_species = report.taxon.name if report.taxon else ""
    if not report.taxon:
        gene_species = st.text_input("Plant species for gene mapping", value="Arabidopsis thaliana",
                                     help="No species context was given; choose a plant (default: model plant).")
    if st.button(f"🧬 Map genes for {c.name}" + (f" in {gene_species}" if gene_species else "")):
        with st.spinner("Gene agent: EC numbers → KEGG GENES / UniProt…"):
            orch = Orchestrator(st.session_state.get("settings", settings), audit)
            taxon = report.taxon or orch.resolve_taxon(gene_species)
            if taxon:
                orch.map_genes(d, taxon)
                st.rerun()
            else:
                st.error(f"Could not resolve '{gene_species}'.")
    st.caption("Phytozome requires a JGI login and has no open API; search the gene IDs above at "
               "[Phytozome](https://phytozome-next.jgi.doe.gov/). UniProt entries include TrEMBL predictions "
               "from sequenced plant genomes.")
    if "genes" in d.errors:
        st.warning(d.errors["genes"])
    provenance_box(d.provenance.get("genes", []))

# --------------------------------------------------------------------------
# Exports and audit
# --------------------------------------------------------------------------
st.divider()
e1, e2, e3 = st.columns(3)
e1.download_button("⬇ Full report (JSON + audit)", rep.to_json(report, audit), f"phytochem_{report.run_id}.json",
                   "application/json", width="stretch")
e2.download_button("⬇ Tables (CSV, zipped)", rep.to_csv_zip(report), f"phytochem_{report.run_id}.zip",
                   "application/zip", width="stretch")
try:
    e3.download_button("⬇ PDF summary", rep.to_pdf(report), f"phytochem_{report.run_id}.pdf", "application/pdf",
                       width="stretch")
except Exception as exc:  # noqa: BLE001
    e3.caption(f"PDF unavailable: {exc}")

with st.expander("🧾 Audit log: every routing decision and API call in this run"):
    if audit:
        ev = pd.DataFrame(audit.events)
        st.caption(f"Run {audit.run_id} · {stats.get('http_calls', 0)} HTTP calls · "
                   f"{stats.get('from_cache', 0)} from cache · {stats.get('errors', 0)} errors · "
                   f"sources: {', '.join(stats.get('sources', []))}")
        cols = [x for x in ["t", "elapsed_s", "kind", "agent", "message", "source", "status", "from_cache", "ms", "url",
                            "reason"] if x in ev.columns]
        st.dataframe(ev[cols].astype(str).replace("nan", ""), hide_index=True, width="stretch")
