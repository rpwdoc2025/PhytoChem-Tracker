"""PubMed E-utilities: paper counts used to rank compounds (the Literature agent)."""

from __future__ import annotations

from .. import config
from ..http import Client

SOURCE = "PubMed E-utilities"


def _quote_term(term: str) -> str:
    return '"' + term.replace('"', "") + '"'


def build_query(compound_names: list[str], species: str) -> str:
    """("quercetin"[tiab] OR "sophoretin"[tiab]) AND "Camellia sinensis"[tiab]

    Searching title/abstract (tiab) keeps counts specific: a compound counts
    only when it is mentioned together with the plant in the same paper.
    """
    names = [n for n in compound_names if n]
    comp = " OR ".join(f"{_quote_term(n)}[tiab]" for n in names)
    if len(names) > 1:
        comp = f"({comp})"
    return f"{comp} AND {_quote_term(species)}[tiab]"


def count(client: Client, term: str, api_key: str = "", email: str = "") -> int:
    params = {"db": "pubmed", "term": term, "retmode": "json", "retmax": 0, "tool": "PhytoChem-Tracker"}
    if api_key:
        params["api_key"] = api_key
    if email:
        params["email"] = email
    data = client.get_json(f"{config.EUTILS_BASE}/esearch.fcgi", source=SOURCE, params=params)
    res = data.get("esearchresult", {})
    if "ERROR" in res:
        raise ValueError(f"PubMed error: {res['ERROR']}")
    return int(res.get("count", 0))
