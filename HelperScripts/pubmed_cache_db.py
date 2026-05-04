"""SQLite-backed PubMed metadata and full-text method cache."""

from __future__ import annotations

import os
import sqlite3
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterable

from HelperScripts.flybase_data import resolve_pubmed_cache_dir

try:
    import fcntl
except ImportError:  # pragma: no cover - fcntl is available on macOS/Linux HPC.
    fcntl = None


PROJECT_ROOT = Path(__file__).resolve().parents[1]
LOCAL_CACHE_DIR = PROJECT_ROOT / ".flai_cache"
DEFAULT_DB_NAME = "pubmed_cache.sqlite"
SCHEMA_VERSION = 1

PUBMED_METADATA_COLUMNS = [
    "pmid",
    "pmcid",
    "title",
    "abstract",
    "year",
    "journal",
    "authors",
    "doi",
    "source",
    "updated_at",
]


def _clean(value: Any) -> str:
    return str(value or "").strip()


def normalize_pmid(value: Any) -> str:
    text = _clean(value)
    return text if text.isdigit() else ""


def authors_to_text(value: Any) -> str:
    if isinstance(value, (list, tuple, set)):
        return "; ".join(_clean(item) for item in value if _clean(item))
    return _clean(value)


def _is_writable_location(path: Path) -> bool:
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        test_path = path.parent / f".{path.name}.{os.getpid()}.write-test"
        test_path.write_text("ok", encoding="utf-8")
        test_path.unlink(missing_ok=True)
        return True
    except OSError:
        return False


def resolve_pubmed_cache_db_path(explicit: str | os.PathLike[str] | None = None) -> tuple[Path, str]:
    """Resolve the PubMed cache DB path, preferring explicit/env, then Turbo, then local."""
    if explicit:
        return Path(explicit).expanduser(), "explicit"

    env_path = os.environ.get("FLAI_PUBMED_CACHE_DB")
    if env_path:
        return Path(env_path).expanduser(), "env"

    turbo_path = resolve_pubmed_cache_dir() / DEFAULT_DB_NAME
    if _is_writable_location(turbo_path):
        return turbo_path, "turbo"

    local_path = LOCAL_CACHE_DIR / DEFAULT_DB_NAME
    local_path.parent.mkdir(parents=True, exist_ok=True)
    return local_path, "local-fallback"


class PubMedCacheDB:
    """Small SQLite repository for PubMed metadata and full-text method caches."""

    def __init__(
        self,
        db_path: str | os.PathLike[str] | None = None,
        *,
        busy_timeout_ms: int = 60_000,
    ):
        self.path, self.source = resolve_pubmed_cache_db_path(db_path)
        self.busy_timeout_ms = int(busy_timeout_ms)
        self._initialized = False

    @contextmanager
    def _write_lock(self):
        lock_path = self.path.with_suffix(self.path.suffix + ".lock")
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        with open(lock_path, "a", encoding="utf-8") as handle:
            if fcntl is not None:
                deadline = time.monotonic() + max(self.busy_timeout_ms / 1000, 1)
                while True:
                    try:
                        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
                        break
                    except BlockingIOError:
                        if time.monotonic() >= deadline:
                            raise
                        time.sleep(0.1)
            try:
                yield
            finally:
                if fcntl is not None:
                    fcntl.flock(handle.fileno(), fcntl.LOCK_UN)

    def connect(self) -> sqlite3.Connection:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        conn = sqlite3.connect(str(self.path), timeout=max(self.busy_timeout_ms / 1000, 1))
        conn.row_factory = sqlite3.Row
        conn.execute(f"PRAGMA busy_timeout={self.busy_timeout_ms}")
        conn.execute("PRAGMA synchronous=NORMAL")
        return conn

    def initialize(self) -> None:
        if self._initialized:
            return
        with self._write_lock():
            if self._initialized:
                return
            with self.connect() as conn:
                conn.executescript(
                    """
                    CREATE TABLE IF NOT EXISTS pubmed_metadata (
                        pmid TEXT PRIMARY KEY,
                        pmcid TEXT NOT NULL DEFAULT '',
                        title TEXT NOT NULL DEFAULT '',
                        abstract TEXT NOT NULL DEFAULT '',
                        year TEXT NOT NULL DEFAULT '',
                        journal TEXT NOT NULL DEFAULT '',
                        authors TEXT NOT NULL DEFAULT '',
                        doi TEXT NOT NULL DEFAULT '',
                        source TEXT NOT NULL DEFAULT '',
                        updated_at TEXT NOT NULL DEFAULT ''
                    );

                    CREATE TABLE IF NOT EXISTS fulltext_method (
                        pmid TEXT PRIMARY KEY,
                        method TEXT NOT NULL DEFAULT '',
                        updated_at TEXT NOT NULL DEFAULT ''
                    );

                    CREATE TABLE IF NOT EXISTS cache_meta (
                        key TEXT PRIMARY KEY,
                        value TEXT NOT NULL DEFAULT ''
                    );
                    """
                )
                conn.execute(
                    "INSERT INTO cache_meta(key, value) VALUES('schema_version', ?) "
                    "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
                    (str(SCHEMA_VERSION),),
                )
        self._initialized = True

    def get_metadata(self, pmid: str) -> dict[str, str] | None:
        pmid_clean = normalize_pmid(pmid)
        if not pmid_clean:
            return None
        self.initialize()
        with self.connect() as conn:
            row = conn.execute(
                "SELECT pmid, pmcid, title, abstract, year, journal, authors, doi, source, updated_at "
                "FROM pubmed_metadata WHERE pmid = ?",
                (pmid_clean,),
            ).fetchone()
        return dict(row) if row else None

    def upsert_metadata(self, entry: dict[str, Any]) -> bool:
        row = {column: _clean(entry.get(column, "")) for column in PUBMED_METADATA_COLUMNS}
        row["pmid"] = normalize_pmid(row["pmid"])
        row["authors"] = authors_to_text(entry.get("authors", row["authors"]))
        if not row["pmid"]:
            return False
        self.initialize()
        with self._write_lock():
            with self.connect() as conn:
                conn.execute(
                    """
                    INSERT INTO pubmed_metadata(
                        pmid, pmcid, title, abstract, year, journal, authors, doi, source, updated_at
                    )
                    VALUES(:pmid, :pmcid, :title, :abstract, :year, :journal, :authors, :doi, :source, :updated_at)
                    ON CONFLICT(pmid) DO UPDATE SET
                        pmcid = CASE WHEN excluded.pmcid != '' THEN excluded.pmcid ELSE pubmed_metadata.pmcid END,
                        title = CASE WHEN excluded.title != '' THEN excluded.title ELSE pubmed_metadata.title END,
                        abstract = CASE WHEN excluded.abstract != '' THEN excluded.abstract ELSE pubmed_metadata.abstract END,
                        year = CASE WHEN excluded.year != '' THEN excluded.year ELSE pubmed_metadata.year END,
                        journal = CASE WHEN excluded.journal != '' THEN excluded.journal ELSE pubmed_metadata.journal END,
                        authors = CASE WHEN excluded.authors != '' THEN excluded.authors ELSE pubmed_metadata.authors END,
                        doi = CASE WHEN excluded.doi != '' THEN excluded.doi ELSE pubmed_metadata.doi END,
                        source = CASE WHEN excluded.source != '' THEN excluded.source ELSE pubmed_metadata.source END,
                        updated_at = CASE WHEN excluded.updated_at != '' THEN excluded.updated_at ELSE pubmed_metadata.updated_at END
                    """,
                    row,
                )
        return True

    def upsert_metadata_many(self, entries: Iterable[dict[str, Any]]) -> int:
        count = 0
        for entry in entries:
            if self.upsert_metadata(entry):
                count += 1
        return count

    def get_fulltext_method(self, pmid: str) -> str | None:
        pmid_clean = normalize_pmid(pmid)
        if not pmid_clean:
            return None
        self.initialize()
        with self.connect() as conn:
            row = conn.execute(
                "SELECT method FROM fulltext_method WHERE pmid = ?",
                (pmid_clean,),
            ).fetchone()
        method = _clean(row["method"]) if row else ""
        return method or None

    def upsert_fulltext_method(self, pmid: str, method: str, updated_at: str = "") -> bool:
        pmid_clean = normalize_pmid(pmid)
        method_clean = _clean(method)
        if not pmid_clean or not method_clean:
            return False
        self.initialize()
        with self._write_lock():
            with self.connect() as conn:
                conn.execute(
                    """
                    INSERT INTO fulltext_method(pmid, method, updated_at)
                    VALUES(?, ?, ?)
                    ON CONFLICT(pmid) DO UPDATE SET
                        method = CASE WHEN excluded.method != '' THEN excluded.method ELSE fulltext_method.method END,
                        updated_at = CASE WHEN excluded.updated_at != '' THEN excluded.updated_at ELSE fulltext_method.updated_at END
                    """,
                    (pmid_clean, method_clean, _clean(updated_at)),
                )
        return True

    def integrity_check(self) -> str:
        self.initialize()
        with self.connect() as conn:
            row = conn.execute("PRAGMA integrity_check").fetchone()
        return str(row[0]) if row else "missing integrity_check result"

    def table_counts(self) -> dict[str, int]:
        self.initialize()
        with self.connect() as conn:
            return {
                "pubmed_metadata": int(conn.execute("SELECT COUNT(*) FROM pubmed_metadata").fetchone()[0]),
                "fulltext_method": int(conn.execute("SELECT COUNT(*) FROM fulltext_method").fetchone()[0]),
            }


_DEFAULT_DB: PubMedCacheDB | None = None


def get_default_pubmed_cache_db() -> PubMedCacheDB:
    global _DEFAULT_DB
    if _DEFAULT_DB is None:
        _DEFAULT_DB = PubMedCacheDB()
        _DEFAULT_DB.initialize()
    return _DEFAULT_DB
