"""HTTP client shared by all connectors.

Features required by the risk-mitigation slide:
* SQLite response cache with a TTL (requests-cache) -> survives API downtime
  and makes repeated queries instant.
* Automatic retries with exponential backoff on 429/5xx, honouring
  ``Retry-After``.
* Per-host rate limiting so we stay inside each provider's request budget.
* Every call is written to the audit log (URL, status, cache hit, timing).
"""

from __future__ import annotations

import threading
import time
from typing import Any
from urllib.parse import urlparse

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from . import config
from .audit import AuditLog, now_iso, redact
from .models import Provenance

try:  # requests-cache is optional; the app still works without it.
    import requests_cache
except ImportError:  # pragma: no cover
    requests_cache = None


class HTTPError(RuntimeError):
    def __init__(self, message: str, status: int | None = None, url: str = ""):
        super().__init__(message)
        self.status = status
        self.url = url


class RateLimiter:
    """Minimal thread-safe per-host limiter (minimum interval between calls)."""

    def __init__(self, limits: dict[str, float]):
        self.limits = dict(limits)
        self._next: dict[str, float] = {}
        self._lock = threading.Lock()

    def wait(self, host: str) -> None:
        rate = self.limits.get(host)
        if not rate:
            return
        interval = 1.0 / rate
        with self._lock:
            now = time.monotonic()
            slot = max(now, self._next.get(host, 0.0))
            self._next[host] = slot + interval
        delay = slot - time.monotonic()
        if delay > 0:
            time.sleep(delay)


_SESSION_LOCK = threading.Lock()
_SESSIONS: dict[tuple[bool, int], requests.Session] = {}


def _build_session(use_cache: bool, cache_hours: int) -> requests.Session:
    if use_cache and requests_cache is not None:
        session: requests.Session = requests_cache.CachedSession(
            str(config.cache_dir() / "http_cache"),
            backend="sqlite",
            expire_after=cache_hours * 3600,
            allowable_codes=(200, 404),  # cache "not found" too: it is an answer
            allowable_methods=("GET", "POST"),
            stale_if_error=True,         # serve stale data if the API is down
        )
    else:
        session = requests.Session()
    retry = Retry(
        total=3,
        backoff_factor=1.0,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET", "POST"),
        respect_retry_after_header=True,
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_maxsize=16)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers["User-Agent"] = config.USER_AGENT
    return session


def get_session(use_cache: bool = True, cache_hours: int = 72) -> requests.Session:
    key = (use_cache, cache_hours)
    with _SESSION_LOCK:
        if key not in _SESSIONS:
            _SESSIONS[key] = _build_session(use_cache, cache_hours)
        return _SESSIONS[key]


class Client:
    """Thin wrapper used by connectors. One instance per run."""

    def __init__(
        self,
        audit: AuditLog | None = None,
        use_cache: bool = True,
        cache_hours: int = 72,
        rate_limits: dict[str, float] | None = None,
        timeout: float = 30.0,
        session: requests.Session | None = None,
    ):
        self.audit = audit or AuditLog(write_file=False)
        self.session = session if session is not None else get_session(use_cache, cache_hours)
        self.limiter = RateLimiter(rate_limits or config.RATE_LIMITS)
        self.timeout = timeout

    def request(
        self,
        method: str,
        url: str,
        *,
        source: str,
        params: dict[str, Any] | None = None,
        data: dict[str, Any] | None = None,
        headers: dict[str, str] | None = None,
        ok_statuses: tuple[int, ...] = (200,),
        timeout: float | None = None,
    ) -> requests.Response:
        host = urlparse(url).netloc
        from_cache = False
        t0 = time.perf_counter()
        status: int | None = None
        try:
            # Only throttle real network calls: check the cache first.
            cached = self._is_cached(method, url, params, data)
            if not cached:
                self.limiter.wait(host)
            resp = self.session.request(
                method, url, params=params, data=data, headers=headers,
                timeout=timeout or self.timeout,
            )
            status = resp.status_code
            from_cache = bool(getattr(resp, "from_cache", False))
            ok = status in ok_statuses
            self.audit.log(
                "http", source=source, method=method, url=resp.url, status=status,
                from_cache=from_cache, ok=ok, ms=round((time.perf_counter() - t0) * 1000),
            )
            if not ok:
                raise HTTPError(f"{source}: HTTP {status}", status=status, url=redact(resp.url))
            return resp
        except requests.RequestException as exc:
            self.audit.log(
                "http", source=source, method=method, url=url, status=status,
                from_cache=False, ok=False, error=str(exc)[:300],
                ms=round((time.perf_counter() - t0) * 1000),
            )
            raise HTTPError(f"{source}: {type(exc).__name__}: {exc}", status=status, url=redact(url)) from exc

    def _is_cached(self, method, url, params, data) -> bool:
        cache = getattr(self.session, "cache", None)
        if cache is None:
            return False
        try:
            req = requests.Request(method, url, params=params, data=data).prepare()
            return bool(cache.contains(request=req))
        except Exception:  # noqa: BLE001 - cache introspection is best-effort
            return False

    def get_json(self, url: str, *, source: str, **kw: Any) -> Any:
        return self.request("GET", url, source=source, **kw).json()

    def get_text(self, url: str, *, source: str, **kw: Any) -> str:
        return self.request("GET", url, source=source, **kw).text

    def provenance(self, source: str, query: str, n: int | None = None, note: str = "") -> Provenance:
        """Build a Provenance record, reusing the cache flag of the last matching call."""
        q = redact(query)
        from_cache = False
        for ev in reversed(self.audit.events):
            if ev.get("kind") == "http" and ev.get("source") == source:
                from_cache = bool(ev.get("from_cache"))
                break
        return Provenance(source=source, query=q, retrieved_at=now_iso(), n=n, from_cache=from_cache, note=note)
