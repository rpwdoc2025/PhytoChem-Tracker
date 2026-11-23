"""Optional LLM layer (Claude, ChatGPT or Gemini).

Hallucination guard (slide 16): the LLM never supplies data. It is used for
two narrow tasks only:

1. ``binomial_for_common_name`` – suggest a binomial for a vernacular name
   ("peppermint" -> "Mentha × piperita"); the suggestion is then verified on
   Wikidata before it is used.
2. ``grounded_summary`` – write a short narrative from the retrieved JSON,
   with instructions to use no facts or numbers that are not in that JSON.

Everything works without an API key; these features are simply skipped.
"""

from __future__ import annotations

import json
import os
import re
import time
from typing import Any

import requests

from .audit import AuditLog
from .config import Settings
from .models import Report

DEFAULT_MODELS = {
    "anthropic": os.getenv("PCT_ANTHROPIC_MODEL", "claude-sonnet-5-5"),
    "openai": os.getenv("PCT_OPENAI_MODEL", "gpt-4.1-mini"),
    "gemini": os.getenv("PCT_GEMINI_MODEL", "gemini-2.5-flash"),
}
KEY_ENV = {"anthropic": "ANTHROPIC_API_KEY", "openai": "OPENAI_API_KEY", "gemini": "GEMINI_API_KEY"}


def api_key_from_env(provider: str) -> str:
    return os.getenv(KEY_ENV.get(provider, ""), "")


def complete(settings: Settings, system: str, user: str, audit: AuditLog | None = None, max_tokens: int = 900) -> tuple[str, str]:
    provider = settings.llm_provider
    key = settings.llm_api_key or api_key_from_env(provider)
    model = settings.llm_model or DEFAULT_MODELS.get(provider, "")
    if not key:
        raise RuntimeError(f"No API key for {provider}")
    t0 = time.perf_counter()
    if provider == "anthropic":
        r = requests.post("https://api.anthropic.com/v1/messages", timeout=90, headers={
            "x-api-key": key, "anthropic-version": "2023-06-01", "content-type": "application/json"},
            json={"model": model, "max_tokens": max_tokens, "system": system,
                  "messages": [{"role": "user", "content": user}]})
        r.raise_for_status()
        text = "".join(b.get("text", "") for b in r.json().get("content", []))
    elif provider == "openai":
        r = requests.post("https://api.openai.com/v1/chat/completions", timeout=90,
                          headers={"Authorization": f"Bearer {key}"},
                          json={"model": model, "max_tokens": max_tokens,
                                "messages": [{"role": "system", "content": system}, {"role": "user", "content": user}]})
        r.raise_for_status()
        text = r.json()["choices"][0]["message"]["content"]
    elif provider == "gemini":
        r = requests.post(f"https://generativelanguage.googleapis.com/v1beta/models/{model}:generateContent",
                          timeout=90, headers={"x-goog-api-key": key},
                          json={"systemInstruction": {"parts": [{"text": system}]},
                                "contents": [{"role": "user", "parts": [{"text": user}]}],
                                "generationConfig": {"maxOutputTokens": max_tokens}})
        r.raise_for_status()
        parts = r.json()["candidates"][0]["content"]["parts"]
        text = "".join(p.get("text", "") for p in parts)
    else:
        raise ValueError(f"Unknown LLM provider '{provider}'")
    if audit:
        audit.log("llm", provider=provider, model=model, prompt_chars=len(system) + len(user),
                  response_chars=len(text), ms=round((time.perf_counter() - t0) * 1000))
    return text.strip(), f"{provider}:{model}"


def binomial_for_common_name(settings: Settings, name: str, audit: AuditLog | None = None) -> str:
    system = ("You map plant common names to accepted botanical binomials. Reply with the binomial only "
              "(e.g. 'Mentha × piperita'), or NONE if unsure.")
    try:
        text, _ = complete(settings, system, f"Common name: {name}", audit, max_tokens=30)
    except Exception:  # noqa: BLE001
        return ""
    text = text.strip().strip(".").strip("'\"")
    return "" if not re.match(r"^[A-Z][a-z]+ (×\s|x\s)?[a-z-]+", text) else text


def _compact(report: Report) -> dict[str, Any]:
    """Only the retrieved facts the summary may use."""
    out: dict[str, Any] = {"input": report.input_value, "input_type": report.input_type}
    if report.taxon:
        out["species"] = report.taxon.name
    out["top_compounds"] = [{"rank": rc.rank, "name": rc.compound.name, "pubmed_papers": rc.paper_count}
                            for rc in report.top]
    comps = []
    for d in list(report.details.values())[:10]:
        comps.append({
            "name": d.compound.name,
            "formula": d.compound.formula,
            "mw": d.compound.mw,
            "melting_point": [p.value for p in d.experimental if p.name == "Melting Point"][:2],
            "solubility": [p.value for p in d.experimental if p.name == "Solubility"][:2],
            "solid_state": d.solid_state,
            "use_categories": sorted({u.category for u in d.uses}),
            "pharmacological_classes": [u.text.split(":")[0] for u in d.uses if u.category == "Pharmacological class"][:5],
            "pathways": [p.name for p in d.pathways][:5],
            "spectra_available": sorted({s.type for s in d.spectra}),
            "n_gene_records": len(d.genes),
        })
    out["compounds"] = comps
    if report.spectrum_matches:
        out["spectrum_matches"] = [{k: m.get(k) for k in ("spectrum", "name", "adduct", "ppm_error", "cosine")}
                                   for m in report.spectrum_matches[:10]]
    return out


def grounded_summary(settings: Settings, report: Report, audit: AuditLog | None = None) -> tuple[str, str]:
    system = (
        "You are a natural-products chemist writing a briefing for a lab meeting. Use ONLY the facts in the JSON "
        "provided. Do not add any compound, number, use, or claim that is not in the JSON; if something is missing, "
        "say it is not available. Write 2 short paragraphs (max 180 words) followed by 3 bullet-point suggestions "
        "for next experimental steps that follow from the data."
    )
    data = json.dumps(_compact(report), default=str)[:12000]
    return complete(settings, system, f"Retrieved data (JSON):\n{data}", audit)
