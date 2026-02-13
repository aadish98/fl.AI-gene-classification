#!/usr/bin/env python3
"""Sanitize and validate PubMed cache CSV files."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

import pandas as pd


_TURBO_CANDIDATES = [Path("/Volumes/umms-rallada"), Path("/nfs/turbo/umms-rallada")]
TURBO_ROOT = next((p for p in _TURBO_CANDIDATES if p.exists()), _TURBO_CANDIDATES[0])
DEFAULT_CACHE_DIR = TURBO_ROOT / "UM Lab Users" / "Aadish" / "Data" / "PubMed Cache"
DEFAULT_PUBMED_CACHE = DEFAULT_CACHE_DIR / "pmid_to_title_abstract.csv"
DEFAULT_FULLTEXT_CACHE = DEFAULT_CACHE_DIR / "pmid_to_fulltext_method.csv"

PUBMED_COLS = ["pmid", "pmcid", "title", "abstract", "year", "journal", "authors", "doi", "source", "updated_at"]
FULLTEXT_COLS = ["pmid", "method", "updated_at"]


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        df = pd.read_csv(path, dtype=str, keep_default_na=False)
    except Exception:
        return pd.DataFrame()
    drop_cols = [c for c in df.columns if str(c).startswith("Unnamed")]
    if drop_cols:
        df = df.drop(columns=drop_cols)
    return df


def _final_nonempty(series: pd.Series) -> str:
    for value in reversed(list(series)):
        text = str(value or "").strip()
        if text:
            return text
    return ""


def sanitize_pubmed_cache(df: pd.DataFrame) -> tuple[pd.DataFrame, Dict]:
    stats = {"rows_in": int(len(df)), "rows_invalid_pmid": 0, "rows_out": 0, "deduped": 0}
    if df.empty:
        return pd.DataFrame(columns=PUBMED_COLS), stats
    for col in PUBMED_COLS:
        if col not in df.columns:
            df[col] = ""
    df = df[PUBMED_COLS].copy()
    df["pmid"] = df["pmid"].apply(lambda x: str(x).strip())
    invalid_mask = ~df["pmid"].str.isdigit()
    stats["rows_invalid_pmid"] = int(invalid_mask.sum())
    df = df[~invalid_mask].copy()

    grouped = []
    for pmid, g in df.groupby("pmid", sort=True):
        grouped.append({
            "pmid": pmid,
            "pmcid": _final_nonempty(g["pmcid"]),
            "title": _final_nonempty(g["title"]),
            "abstract": _final_nonempty(g["abstract"]),
            "year": _final_nonempty(g["year"]),
            "journal": _final_nonempty(g["journal"]),
            "authors": _final_nonempty(g["authors"]),
            "doi": _final_nonempty(g["doi"]),
            "source": _final_nonempty(g["source"]),
            "updated_at": _final_nonempty(g["updated_at"]),
        })
    out = pd.DataFrame(grouped, columns=PUBMED_COLS)
    stats["rows_out"] = int(len(out))
    stats["deduped"] = int(len(df) - len(out))
    return out, stats


def sanitize_fulltext_cache(df: pd.DataFrame) -> tuple[pd.DataFrame, Dict]:
    stats = {"rows_in": int(len(df)), "rows_invalid_pmid": 0, "rows_out": 0, "deduped": 0}
    if df.empty:
        return pd.DataFrame(columns=FULLTEXT_COLS), stats
    for col in FULLTEXT_COLS:
        if col not in df.columns:
            df[col] = ""
    df = df[FULLTEXT_COLS].copy()
    df["pmid"] = df["pmid"].apply(lambda x: str(x).strip())
    invalid_mask = ~df["pmid"].str.isdigit()
    stats["rows_invalid_pmid"] = int(invalid_mask.sum())
    df = df[~invalid_mask].copy()
    grouped = []
    for pmid, g in df.groupby("pmid", sort=True):
        grouped.append({
            "pmid": pmid,
            "method": _final_nonempty(g["method"]),
            "updated_at": _final_nonempty(g["updated_at"]),
        })
    out = pd.DataFrame(grouped, columns=FULLTEXT_COLS)
    stats["rows_out"] = int(len(out))
    stats["deduped"] = int(len(df) - len(out))
    return out, stats


def main():
    parser = argparse.ArgumentParser(description="Sanitize and validate PubMed cache CSV files.")
    parser.add_argument("--pubmed-cache-path", default=str(DEFAULT_PUBMED_CACHE))
    parser.add_argument("--fulltext-cache-path", default=str(DEFAULT_FULLTEXT_CACHE))
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument("--apply", action="store_true", help="Write sanitized files in-place.")
    mode_group.add_argument("--dry-run", action="store_true", help="Preview sanitation only (default).")
    parser.add_argument(
        "--report-path",
        default="Logs/sanitize_pubmed_caches_report.json",
        help="Path for JSON report output.",
    )
    args = parser.parse_args()

    pubmed_path = Path(args.pubmed_cache_path).expanduser()
    fulltext_path = Path(args.fulltext_cache_path).expanduser()
    report_path = Path(args.report_path).expanduser()

    pubmed_df = _read_csv(pubmed_path)
    fulltext_df = _read_csv(fulltext_path)
    pubmed_clean, pubmed_stats = sanitize_pubmed_cache(pubmed_df)
    fulltext_clean, fulltext_stats = sanitize_fulltext_cache(fulltext_df)

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "mode": "apply" if args.apply else "dry-run",
        "files": {
            "pubmed_cache": {
                "path": str(pubmed_path),
                **pubmed_stats,
            },
            "fulltext_cache": {
                "path": str(fulltext_path),
                **fulltext_stats,
            },
        },
    }

    if args.apply:
        if pubmed_path.parent:
            pubmed_path.parent.mkdir(parents=True, exist_ok=True)
        if fulltext_path.parent:
            fulltext_path.parent.mkdir(parents=True, exist_ok=True)
        pubmed_clean.to_csv(pubmed_path, index=False)
        fulltext_clean.to_csv(fulltext_path, index=False)

    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=True, indent=2)

    print(f"Mode: {report['mode']}")
    print(
        f"PubMed rows in/out: {pubmed_stats['rows_in']}/{pubmed_stats['rows_out']} "
        f"(invalid={pubmed_stats['rows_invalid_pmid']} deduped={pubmed_stats['deduped']})"
    )
    print(
        f"Fulltext rows in/out: {fulltext_stats['rows_in']}/{fulltext_stats['rows_out']} "
        f"(invalid={fulltext_stats['rows_invalid_pmid']} deduped={fulltext_stats['deduped']})"
    )
    print(f"Report: {report_path}")


if __name__ == "__main__":
    main()
