import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def pytest_configure(config):
    # Keep tests hermetic: cache and audit files go to a temp folder.
    os.environ.setdefault("PCT_CACHE_DIR", str(Path(os.getenv("TMPDIR", "/tmp")) / "pct_test_cache"))
