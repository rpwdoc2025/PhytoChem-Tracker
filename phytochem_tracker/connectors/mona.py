"""MassBank of North America (MoNA): reference MS/MS spectra by InChIKey.

MoNA moved to the Spring Filter query language in 2023; older deployments
used RSQL. We try the current syntax first and fall back to the legacy one,
so the connector keeps working on mirrors. Failures are non-fatal: the tile
still shows deep links to MoNA, MassBank and GNPS.
"""

from __future__ import annotations

from typing import Any

from .. import config
from ..http import Client, HTTPError
from ..spectra import parse_peak_string

SOURCE = "MoNA REST"


def _queries(inchikey: str) -> list[str]:
    return [
        f"exists(compound.metaData.name:'InChIKey' and compound.metaData.value:'{inchikey}')",
        f"compound.metaData=q='name==\"InChIKey\" and value==\"{inchikey}\"'",
    ]


def _meta(items: list[dict[str, Any]] | None) -> dict[str, str]:
    out = {}
    for m in items or []:
        name, value = m.get("name"), m.get("value")
        if name and value is not None and name not in out:
            out[name] = str(value)
    return out


def parse_records(data: Any) -> list[dict[str, Any]]:
    records = data.get("content", data) if isinstance(data, dict) else data
    out = []
    for rec in records or []:
        meta = _meta(rec.get("metaData"))
        comp = (rec.get("compound") or [{}])[0]
        names = comp.get("names") or [{}]
        out.append({
            "id": rec.get("id", ""),
            "name": names[0].get("name", ""),
            "ms_level": meta.get("ms level", ""),
            "precursor_type": meta.get("precursor type", ""),
            "precursor_mz": meta.get("precursor m/z", ""),
            "ionization": meta.get("ionization mode", ""),
            "instrument": meta.get("instrument type", meta.get("instrument", "")),
            "collision_energy": meta.get("collision energy", ""),
            "peaks": parse_peak_string(rec.get("spectrum", "")),
            "url": f"https://mona.fiehnlab.ucdavis.edu/spectra/display/{rec.get('id', '')}",
        })
    return out


def spectra_for_inchikey(client: Client, inchikey: str, size: int = 10) -> list[dict[str, Any]]:
    last: Exception | None = None
    for q in _queries(inchikey):
        try:
            data = client.get_json(f"{config.MONA_BASE}/spectra/search", source=SOURCE,
                                   params={"query": q, "size": size, "page": 0}, timeout=40)
            recs = parse_records(data)
            if recs:
                return recs
        except (HTTPError, ValueError) as exc:
            last = exc
    if last:
        raise HTTPError(f"MoNA search failed: {last}")
    return []
