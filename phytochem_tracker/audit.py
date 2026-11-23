"""Audit trail: one JSONL file per run with every decision and HTTP call.

Faculty asked for reproducibility and post-docs for "no black box": the audit
log records what the router decided and why, every query sent to every
database, whether it was served from cache, and how long it took.
"""

from __future__ import annotations

import json
import re
import threading
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from . import config

_SECRET_PARAMS = re.compile(r"(api_key|apikey|key|token)=([^&]+)", re.I)


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def redact(url: str) -> str:
    """Remove API keys from URLs before they are logged or displayed."""
    return _SECRET_PARAMS.sub(r"\1=REDACTED", url)


class AuditLog:
    def __init__(self, run_id: str | None = None, write_file: bool = True):
        self.run_id = run_id or datetime.now().strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:6]
        self.events: list[dict[str, Any]] = []
        self._lock = threading.Lock()
        self._t0 = time.perf_counter()
        self.path: Path | None = (config.log_dir() / f"{self.run_id}.jsonl") if write_file else None

    def log(self, kind: str, **data: Any) -> None:
        event = {"t": now_iso(), "elapsed_s": round(time.perf_counter() - self._t0, 3), "kind": kind, **data}
        if "url" in event:
            event["url"] = redact(str(event["url"]))
        with self._lock:
            self.events.append(event)
            if self.path:
                try:
                    with self.path.open("a", encoding="utf-8") as fh:
                        fh.write(json.dumps(event, default=str) + "\n")
                except OSError:
                    self.path = None  # read-only filesystem: keep in memory only

    def http_calls(self) -> list[dict[str, Any]]:
        return [e for e in self.events if e["kind"] == "http"]

    def stats(self) -> dict[str, Any]:
        calls = self.http_calls()
        return {
            "http_calls": len(calls),
            "from_cache": sum(1 for c in calls if c.get("from_cache")),
            "errors": sum(1 for c in calls if not c.get("ok", True)),
            "sources": sorted({c.get("source", "?") for c in calls}),
        }
