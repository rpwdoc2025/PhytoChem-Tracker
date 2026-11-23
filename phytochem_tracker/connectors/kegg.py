"""KEGG REST: compound -> pathways and enzymes (EC); EC -> genes in an organism.

KEGG's REST API returns tab-separated text. Academic use of the REST API is
free; see https://www.kegg.jp/kegg/legal.html.
"""

from __future__ import annotations

from urllib.parse import quote

from .. import config
from ..http import Client, HTTPError

SOURCE = "KEGG REST"

# Global/overview maps are less informative than specific biosynthetic maps.
GLOBAL_MAPS = {"map01100", "map01110", "map01120", "map01060", "map01061", "map01062", "map01063",
               "map01064", "map01065", "map01066", "map01070", "map01200", "map01210", "map01212",
               "map01230", "map01232", "map01240", "map01250", "map01310", "map01320"}


def _lines(text: str) -> list[list[str]]:
    return [ln.split("\t") for ln in text.strip().splitlines() if ln.strip()]


def _get(client: Client, path: str) -> str:
    resp = client.request("GET", f"{config.KEGG_BASE}/{path}", source=SOURCE, ok_statuses=(200, 404))
    return resp.text if resp.status_code == 200 else ""


def find_compound(client: Client, names: list[str], formula: str = "") -> str:
    """Return a KEGG compound ID (e.g. C00389) matching a name, checked against the formula."""
    candidates: dict[str, str] = {}
    if formula:
        for row in _lines(_get(client, f"find/compound/{quote(formula)}/formula")):
            if len(row) >= 2:
                candidates[row[0].replace("cpd:", "")] = row[1]
    wanted = {n.lower() for n in names if n}
    if candidates:
        ids = list(candidates)[:10]
        for row in _lines(_get(client, "list/" + "+".join(f"cpd:{i}" for i in ids))):
            if len(row) >= 2:
                cid = row[0].replace("cpd:", "")
                kegg_names = {n.strip().lower() for n in row[1].split(";")}
                if kegg_names & wanted:
                    return cid
    for name in names[:2]:
        for row in _lines(_get(client, f"find/compound/{quote(name)}")):
            if len(row) >= 2 and name.lower() in {n.strip().lower() for n in row[1].split(";")}:
                cid = row[0].replace("cpd:", "")
                if not candidates or cid in candidates:
                    return cid
    return ""


def pathway_names(client: Client, org: str = "") -> dict[str, str]:
    """{'map00941': 'Flavonoid biosynthesis', ...} (or org-specific IDs when ``org`` is given)."""
    text = _get(client, f"list/pathway/{org}" if org else "list/pathway")
    return {r[0].replace("path:", ""): r[1] for r in _lines(text) if len(r) >= 2}


def pathways_for_compound(client: Client, kegg_id: str) -> list[str]:
    rows = _lines(_get(client, f"link/pathway/cpd:{kegg_id}"))
    maps = [r[1].replace("path:", "") for r in rows if len(r) >= 2]
    specific = [m for m in maps if m not in GLOBAL_MAPS]
    return specific + [m for m in maps if m in GLOBAL_MAPS]


def enzymes_for_compound(client: Client, kegg_id: str) -> list[str]:
    rows = _lines(_get(client, f"link/enzyme/cpd:{kegg_id}"))
    return sorted({r[1].replace("ec:", "") for r in rows if len(r) >= 2})


def enzymes_for_pathway(client: Client, map_id: str) -> list[str]:
    rows = _lines(_get(client, f"link/enzyme/path:{map_id}"))
    return sorted({r[1].replace("ec:", "") for r in rows if len(r) >= 2})


def enzyme_names(client: Client, ecs: list[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for i in range(0, len(ecs), 10):
        chunk = ecs[i:i + 10]
        for row in _lines(_get(client, "list/" + "+".join(f"ec:{e}" for e in chunk))):
            if len(row) >= 2:
                out[row[0].replace("ec:", "")] = row[1].split(";")[0].strip()
    return out


def organism_code(client: Client, species: str) -> str:
    """KEGG three/four-letter organism code for a binomial (e.g. 'Arabidopsis thaliana' -> 'ath')."""
    target = species.lower().replace(" × ", " x ")
    try:
        text = _get(client, "list/organism")
    except HTTPError:
        return ""
    for row in _lines(text):
        if len(row) >= 3:
            name = row[2].lower()
            if name == target or name.startswith(target + " "):
                return row[1]
    return ""


def genes_for_ec(client: Client, org: str, ec: str) -> list[str]:
    rows = _lines(_get(client, f"link/{org}/ec:{ec}"))
    return [r[1] for r in rows if len(r) >= 2]
