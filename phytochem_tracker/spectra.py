"""Mass-spectrometry utilities: file parsing, precursor matching, similarity.

Workflow for an uploaded MS/MS file:
1. ``parse_spectra`` reads MGF / MSP (built-in parsers) or mzML (pyteomics).
2. ``match_precursor`` compares the precursor m/z with the monoisotopic
   masses of candidate compounds (e.g. every LOTUS compound reported in the
   plant) for common adducts, within a ppm tolerance: classic dereplication.
3. ``cosine_similarity`` scores the uploaded MS/MS against reference spectra
   (MoNA) with the standard greedy peak-matching cosine.
"""

from __future__ import annotations

import io
import math
from dataclasses import dataclass, field
from typing import Any, Iterable

from .chem import ADDUCTS_NEG, ADDUCTS_POS, adduct_mz


@dataclass
class Spectrum:
    title: str
    precursor_mz: float | None
    charge: int | None = None
    polarity: str = ""          # "positive" / "negative" / ""
    rt: float | None = None
    peaks: list[tuple[float, float]] = field(default_factory=list)
    meta: dict[str, str] = field(default_factory=dict)


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------
def _charge(text: str) -> tuple[int | None, str]:
    t = text.strip()
    if not t:
        return None, ""
    sign = "-" if t.endswith("-") or t.startswith("-") else "+"
    digits = "".join(ch for ch in t if ch.isdigit()) or "1"
    return int(digits) * (-1 if sign == "-" else 1), ("negative" if sign == "-" else "positive")


def parse_mgf(text: str) -> list[Spectrum]:
    spectra: list[Spectrum] = []
    cur: dict[str, Any] | None = None
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith(("#", ";", "!")):
            continue
        upper = line.upper()
        if upper == "BEGIN IONS":
            cur = {"meta": {}, "peaks": []}
            continue
        if upper == "END IONS":
            if cur is not None:
                meta = cur["meta"]
                pep = meta.get("PEPMASS", "").split()
                charge, pol = _charge(meta.get("CHARGE", ""))
                ion = meta.get("IONMODE", meta.get("POLARITY", "")).lower()
                if ion.startswith("pos"):
                    pol = "positive"
                elif ion.startswith("neg"):
                    pol = "negative"
                rt = meta.get("RTINSECONDS") or meta.get("RTINMINUTES")
                spectra.append(Spectrum(
                    title=meta.get("TITLE", meta.get("NAME", f"spectrum {len(spectra) + 1}")),
                    precursor_mz=float(pep[0]) if pep else None,
                    charge=charge, polarity=pol,
                    rt=float(rt) if rt else None,
                    peaks=sorted(cur["peaks"]), meta=meta,
                ))
            cur = None
            continue
        if cur is None:
            continue
        if "=" in line and not line[0].isdigit():
            k, v = line.split("=", 1)
            cur["meta"][k.strip().upper()] = v.strip()
        else:
            parts = line.replace(",", " ").split()
            try:
                cur["peaks"].append((float(parts[0]), float(parts[1]) if len(parts) > 1 else 1.0))
            except (ValueError, IndexError):
                pass
    return spectra


def parse_msp(text: str) -> list[Spectrum]:
    spectra: list[Spectrum] = []
    block: list[str] = []
    for line in text.splitlines() + [""]:
        if line.strip():
            block.append(line.strip())
            continue
        if not block:
            continue
        meta: dict[str, str] = {}
        peaks: list[tuple[float, float]] = []
        in_peaks = False
        for ln in block:
            if in_peaks or (ln[0].isdigit() and ":" not in ln):
                for chunk in ln.replace(";", "\n").splitlines():
                    p = chunk.replace(",", " ").split()
                    if len(p) >= 2:
                        try:
                            peaks.append((float(p[0]), float(p[1])))
                        except ValueError:
                            pass
                continue
            k, _, v = ln.partition(":")
            meta[k.strip().upper()] = v.strip()
            if k.strip().upper() == "NUM PEAKS":
                in_peaks = True
        prec = meta.get("PRECURSORMZ") or meta.get("PRECURSOR_MZ")
        mode = meta.get("ION_MODE", meta.get("IONMODE", "")).lower()
        spectra.append(Spectrum(
            title=meta.get("NAME", f"spectrum {len(spectra) + 1}"),
            precursor_mz=float(prec) if prec else None,
            polarity="positive" if mode.startswith("p") else "negative" if mode.startswith("n") else "",
            peaks=sorted(peaks), meta=meta,
        ))
        block = []
    return spectra


def parse_mzml(data: bytes, max_spectra: int = 500) -> list[Spectrum]:
    from pyteomics import mzml  # imported lazily: only needed for mzML

    out: list[Spectrum] = []
    with mzml.read(io.BytesIO(data)) as reader:
        for spec in reader:
            if spec.get("ms level", 1) < 2:
                continue
            try:
                ion = spec["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]
                prec = float(ion.get("selected ion m/z"))
            except (KeyError, IndexError, TypeError):
                prec = None
            pol = "positive" if "positive scan" in spec else "negative" if "negative scan" in spec else ""
            peaks = list(zip(map(float, spec["m/z array"]), map(float, spec["intensity array"])))
            out.append(Spectrum(title=spec.get("id", f"scan {len(out) + 1}"), precursor_mz=prec,
                                polarity=pol, peaks=peaks))
            if len(out) >= max_spectra:
                break
    return out


def parse_spectra(filename: str, data: bytes) -> list[Spectrum]:
    name = filename.lower()
    if name.endswith(".mzml"):
        return parse_mzml(data)
    text = data.decode("utf-8", errors="replace")
    if name.endswith(".msp"):
        return parse_msp(text)
    return parse_mgf(text)


def parse_peak_string(s: str) -> list[tuple[float, float]]:
    """MoNA/MassBank 'mz:int mz:int' format."""
    peaks = []
    for tok in (s or "").split():
        if ":" in tok:
            mz, it = tok.split(":", 1)
            try:
                peaks.append((float(mz), float(it)))
            except ValueError:
                pass
    return sorted(peaks)


# --------------------------------------------------------------------------
# Dereplication and similarity
# --------------------------------------------------------------------------
def ppm_error(observed: float, theoretical: float) -> float:
    return (observed - theoretical) / theoretical * 1e6


def match_precursor(
    precursor_mz: float,
    candidates: Iterable[dict[str, Any]],
    polarity: str = "",
    ppm: float = 10.0,
) -> list[dict[str, Any]]:
    """Match a precursor m/z against candidate neutral monoisotopic masses.

    ``candidates`` are dicts with at least ``name`` and ``monoisotopic_mass``.
    Returns matches sorted by absolute ppm error.
    """
    adducts = ADDUCTS_POS if polarity == "positive" else ADDUCTS_NEG if polarity == "negative" else ADDUCTS_POS + ADDUCTS_NEG
    hits = []
    for cand in candidates:
        m = cand.get("monoisotopic_mass")
        if not m:
            continue
        for add in adducts:
            theo = adduct_mz(m, add)
            err = ppm_error(precursor_mz, theo)
            if abs(err) <= ppm:
                hits.append({**cand, "adduct": add[0], "theoretical_mz": round(theo, 5), "ppm_error": round(err, 2)})
    hits.sort(key=lambda h: abs(h["ppm_error"]))
    return hits


def _normalise(peaks: list[tuple[float, float]], precursor_mz: float | None = None,
               min_rel: float = 0.01, top_n: int = 100) -> list[tuple[float, float]]:
    """Square-root intensities, drop the precursor region and noise, keep top N."""
    if not peaks:
        return []
    pk = [(mz, it) for mz, it in peaks if it > 0 and (precursor_mz is None or mz < precursor_mz - 1.5)]
    if not pk:
        return []
    mx = max(it for _, it in pk)
    pk = [(mz, it) for mz, it in pk if it / mx >= min_rel]
    pk = sorted(pk, key=lambda p: -p[1])[:top_n]
    return sorted((mz, math.sqrt(it)) for mz, it in pk)


def cosine_similarity(a: list[tuple[float, float]], b: list[tuple[float, float]],
                      tolerance: float = 0.02, precursor_a: float | None = None,
                      precursor_b: float | None = None) -> tuple[float, int]:
    """Greedy cosine score (as in GNPS/matchms CosineGreedy). Returns (score, matched peaks)."""
    pa, pb = _normalise(a, precursor_a), _normalise(b, precursor_b)
    if not pa or not pb:
        return 0.0, 0
    pairs = []
    j0 = 0
    for i, (mza, ia) in enumerate(pa):
        while j0 < len(pb) and pb[j0][0] < mza - tolerance:
            j0 += 1
        j = j0
        while j < len(pb) and pb[j][0] <= mza + tolerance:
            pairs.append((ia * pb[j][1], i, j))
            j += 1
    pairs.sort(reverse=True)
    used_a, used_b, score = set(), set(), 0.0
    for prod, i, j in pairs:
        if i in used_a or j in used_b:
            continue
        used_a.add(i)
        used_b.add(j)
        score += prod
    norm = math.sqrt(sum(i * i for _, i in pa)) * math.sqrt(sum(i * i for _, i in pb))
    return (round(score / norm, 4) if norm else 0.0), len(used_a)
