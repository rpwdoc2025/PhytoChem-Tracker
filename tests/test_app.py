"""UI smoke tests: drive the Streamlit app headlessly in offline-demo mode."""

from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

APP = str(Path(__file__).resolve().parents[1] / "app.py")


def run_app(mode, query=None, context=None):
    at = AppTest.from_file(APP, default_timeout=180)
    at.run()
    at.sidebar.radio(key="mode").set_value(mode).run()
    for cb in at.sidebar.checkbox:
        if cb.label.startswith("Offline demo"):
            cb.check()
    if query is not None:
        at.sidebar.text_input(key="query").set_value(query)
    if context is not None:
        at.sidebar.text_input(key="context").set_value(context)
    next(b for b in at.sidebar.button if "Run agents" in b.label).click().run()
    return at


def test_landing_page():
    at = AppTest.from_file(APP, default_timeout=60).run()
    assert not at.exception
    assert any("Data sources" in m.value for m in at.markdown)


@pytest.mark.parametrize("mode,query,context,n_compounds", [
    ("Plant species", "Camellia sinensis", None, "4"),
    ("Compound", "quercetin", "", "1"),
    ("MS/MS spectrum", None, None, "2"),
])
def test_modes_render_seven_tiles(mode, query, context, n_compounds):
    at = run_app(mode, query, context)
    assert not at.exception, [e.value for e in at.exception]
    assert len(at.tabs) == 7
    assert at.metric[0].value == n_compounds


def test_switch_compound_and_map_genes():
    at = run_app("Plant species", "Camellia sinensis")
    next(s for s in at.selectbox if s.label.startswith("Compound shown")).select_index(3).run()
    next(b for b in at.button if b.label.startswith("🧬")).click().run()
    assert not at.exception
    assert any(c.value.startswith("Status: ") and "gene/protein records" in c.value for c in at.caption)
