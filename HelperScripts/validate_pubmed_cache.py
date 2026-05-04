#!/usr/bin/env python3
"""Validate legacy PubMed CSV caches and migrate valid rows into SQLite."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from HelperScripts.flybase_data import resolve_pubmed_cache_dir
from HelperScripts.metadata_resolver import normalize_authors
from HelperScripts.pubmed_cache_db import (
    PUBMED_METADATA_COLUMNS,
    PubMedCacheDB,
    authors_to_text,
    normalize_pmid,
    resolve_pubmed_cache_db_path,
)


def _clean(value: Any) -> str:
    return str(value or "").strip()


def _completeness_score(entry: dict[str, Any]) -> tuple[int, str]:
    fields = ["pmcid", "title", "abstract", "year", "journal", "authors", "doi", "source"]
    return sum(1 for field in fields if _clean(entry.get(field))), _clean(entry.get("updated_at"))


def _merge_entry(existing: dict[str, str] | None, incoming: dict[str, str]) -> dict[str, str]:
    if existing is None:
        return dict(incoming)
    winner, loser = (incoming, existing) if _completeness_score(incoming) >= _completeness_score(existing) else (existing, incoming)
    merged = dict(winner)
    for key, value in loser.items():
        if not _clean(merged.get(key)) and _clean(value):
            merged[key] = _clean(value)
    return merged


def _read_csv(path: Path, columns: list[str]) -> tuple[list[dict[str, str]], list[dict[str, Any]]]:
    if not path.exists():
        return [], [{"path": str(path), "reason": "missing_file"}]
    try:
        df = pd.read_csv(path, dtype=str, keep_default_na=False)
    except Exception as exc:
        return [], [{"path": str(path), "reason": "read_error", "error": str(exc)}]
    df = df.drop(columns=[c for c in df.columns if str(c).startswith("Unnamed")], errors="ignore")
    for column in columns:
        if column not in df.columns:
            df[column] = ""
    rows = []
    invalid = []
    for idx, raw in df[columns].iterrows():
        pmid = normalize_pmid(raw.get("pmid", ""))
        if not pmid:
            invalid.append({"path": str(path), "row_index": int(idx), "reason": "invalid_pmid"})
            continue
        entry = {column: _clean(raw.get(column, "")) for column in columns}
        entry["pmid"] = pmid
        rows.append(entry)
    return rows, invalid


def load_legacy_metadata(path: Path) -> tuple[dict[str, dict[str, str]], list[dict[str, Any]], int]:
    rows, invalid = _read_csv(path, PUBMED_METADATA_COLUMNS)
    deduped: dict[str, dict[str, str]] = {}
    duplicates = 0
    for row in rows:
        row["authors"] = authors_to_text(normalize_authors(row.get("authors", "")))
        pmid = row["pmid"]
        if pmid in deduped:
            duplicates += 1
        deduped[pmid] = _merge_entry(deduped.get(pmid), row)
    return deduped, invalid, duplicates


def load_legacy_fulltext(path: Path) -> tuple[dict[str, dict[str, str]], list[dict[str, Any]], int]:
    rows, invalid = _read_csv(path, ["pmid", "method", "updated_at"])
    deduped: dict[str, dict[str, str]] = {}
    duplicates = 0
    for row in rows:
        if not row.get("method"):
            invalid.append({"path": str(path), "pmid": row["pmid"], "reason": "blank_method"})
            continue
        pmid = row["pmid"]
        if pmid in deduped:
            duplicates += 1
        deduped[pmid] = _merge_entry(deduped.get(pmid), row)
    return deduped, invalid, duplicates


def _verify_metadata(db: PubMedCacheDB, entries: dict[str, dict[str, str]]) -> list[dict[str, Any]]:
    failures = []
    for pmid, expected in entries.items():
        found = db.get_metadata(pmid)
        if not found:
            failures.append({"pmid": pmid, "reason": "missing_from_sqlite"})
            continue
        for field, value in expected.items():
            if field == "pmid" or not _clean(value):
                continue
            if _clean(found.get(field)) != _clean(value):
                failures.append({"pmid": pmid, "field": field, "reason": "value_mismatch"})
                break
    return failures


def _verify_fulltext(db: PubMedCacheDB, entries: dict[str, dict[str, str]]) -> list[dict[str, Any]]:
    failures = []
    for pmid, expected in entries.items():
        found = db.get_fulltext_method(pmid)
        if _clean(found) != _clean(expected.get("method")):
            failures.append({"pmid": pmid, "reason": "fulltext_method_mismatch"})
    return failures


def run(args: argparse.Namespace) -> dict[str, Any]:
    default_cache_dir = resolve_pubmed_cache_dir()
    metadata_csv = Path(args.metadata_csv or default_cache_dir / "pmid_to_title_abstract.csv").expanduser()
    fulltext_csv = Path(args.fulltext_csv or default_cache_dir / "pmid_to_fulltext_method.csv").expanduser()
    db_path, db_source = resolve_pubmed_cache_db_path(args.db_path or None)
    db = PubMedCacheDB(db_path)

    metadata, invalid_metadata, metadata_duplicates = load_legacy_metadata(metadata_csv)
    fulltext, invalid_fulltext, fulltext_duplicates = load_legacy_fulltext(fulltext_csv)

    report: dict[str, Any] = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "apply": bool(args.apply),
        "metadata_csv": str(metadata_csv),
        "fulltext_csv": str(fulltext_csv),
        "sqlite_db": str(db.path),
        "sqlite_db_source": db_source,
        "metadata_valid_rows": len(metadata),
        "metadata_duplicate_rows": metadata_duplicates,
        "metadata_invalid_rows": invalid_metadata,
        "fulltext_valid_rows": len(fulltext),
        "fulltext_duplicate_rows": fulltext_duplicates,
        "fulltext_invalid_rows": invalid_fulltext,
        "sqlite_integrity_before": "not_checked",
        "sqlite_integrity_after": "not_checked",
        "verification_failures": [],
        "table_counts": {},
    }

    if db.path.exists() or args.apply:
        try:
            report["sqlite_integrity_before"] = db.integrity_check()
        except Exception as exc:
            report["sqlite_integrity_before"] = f"error: {exc}"

    if args.apply:
        db.initialize()
        for entry in metadata.values():
            db.upsert_metadata(entry)
        for entry in fulltext.values():
            db.upsert_fulltext_method(entry["pmid"], entry["method"], entry.get("updated_at", ""))
        report["verification_failures"].extend(_verify_metadata(db, metadata))
        report["verification_failures"].extend(_verify_fulltext(db, fulltext))
        report["sqlite_integrity_after"] = db.integrity_check()
        report["table_counts"] = db.table_counts()

    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Validate and migrate PubMed cache CSVs into SQLite.")
    parser.add_argument("--metadata-csv", default="", help="Legacy pmid_to_title_abstract.csv path.")
    parser.add_argument("--fulltext-csv", default="", help="Legacy pmid_to_fulltext_method.csv path.")
    parser.add_argument("--db-path", default="", help="SQLite DB path; defaults to FLAI/Turbo/local resolver.")
    parser.add_argument("--apply", action="store_true", help="Apply validated rows to SQLite.")
    parser.add_argument("--fail-on-invalid", action="store_true", help="Exit nonzero if invalid rows are found.")
    parser.add_argument("--report-path", default="", help="Optional JSON report path.")
    args = parser.parse_args(argv)

    report = run(args)
    payload = json.dumps(report, indent=2, sort_keys=True)
    if args.report_path:
        report_path = Path(args.report_path).expanduser()
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(payload + "\n", encoding="utf-8")
    print(payload)

    has_invalid = bool(report["metadata_invalid_rows"] or report["fulltext_invalid_rows"])
    has_failures = bool(report["verification_failures"])
    integrity_bad = any(
        value not in {"ok", "not_checked"}
        for value in (report["sqlite_integrity_before"], report["sqlite_integrity_after"])
    )
    if has_failures or integrity_bad or (args.fail_on_invalid and has_invalid):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
