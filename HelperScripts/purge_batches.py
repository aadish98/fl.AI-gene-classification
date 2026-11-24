#!/usr/bin/env python3
"""Purge bad OpenAI checkpoint batches so non-force reruns reprocess them.

This utility is intentionally dry-run by default. With --apply it edits
.batch_state metadata and removes matching run-store records for genes whose
saved output contains quota/credit/API-failure markers.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]

DEFAULT_INPUT_DIRS = [
    PROJECT_ROOT / "Data" / "Tx-Omics-Conserved-Correlates-Large" / "Set1",
    PROJECT_ROOT / "Data" / "Tx-Omics-Conserved-Correlates-Large" / "Set2",
    PROJECT_ROOT / "Data" / "Tx-Omics-Conserved-CSW-Large" / "Set1",
    PROJECT_ROOT / "Data" / "Tx-Omics-Conserved-CSW-Large" / "Set2",
]

DEFAULT_MARKERS = [
    "insufficient_quota",
    "quota",
    "credit",
    "billing",
    "rate limit",
    "rate_limit",
    "API error",
    "Error:",
]


@dataclass
class BadBatch:
    input_dir: Path
    state_dir: Path
    state_path: Path
    batch_path: Path
    batch_index: int
    species: str
    bad_gene_ids: set[str] = field(default_factory=set)
    reasons: list[str] = field(default_factory=list)


def _utc_stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def _read_json(path: Path, default: Any) -> Any:
    try:
        with path.open("r", encoding="utf-8") as handle:
            return json.load(handle)
    except FileNotFoundError:
        return default
    except Exception as exc:
        print(f"[warn] Could not read JSON {path}: {exc}", file=sys.stderr)
        return default


def _atomic_json_dump(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(f".{path.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp")
    try:
        with tmp_path.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=True, indent=2, sort_keys=True)
        tmp_path.replace(path)
    finally:
        tmp_path.unlink(missing_ok=True)


def _backup(path: Path, timestamp: str) -> Path | None:
    if not path.exists():
        return None
    backup_path = path.with_name(f"{path.name}.bak_{timestamp}")
    shutil.copy2(path, backup_path)
    return backup_path


def _clean(value: Any) -> str:
    return str(value or "").strip()


def _gene_id(record: dict[str, Any]) -> str:
    payload = record.get("gene_record")
    if isinstance(payload, dict):
        return _clean(payload.get("gene_id") or payload.get("fbgn_id"))
    return _clean(record.get("gene_id") or record.get("fbgn_id"))


def _selected_record_text(record: dict[str, Any]) -> str:
    parts: list[str] = []
    summaries = record.get("summaries")
    if isinstance(summaries, list):
        for summary in summaries:
            if not isinstance(summary, dict):
                continue
            parts.append(_clean(summary.get("qc_justification")))
            parts.append(_clean(summary.get("summary")))
    hit_info = record.get("hit_info")
    if isinstance(hit_info, dict):
        parts.append(_clean(hit_info.get("rationale")))
    return "\n".join(part for part in parts if part)


def _record_bad_reasons(record: dict[str, Any], markers: list[str]) -> list[str]:
    text = _selected_record_text(record)
    if not text:
        return []
    text_lower = text.lower()
    out = []
    for marker in markers:
        marker_clean = marker.strip()
        if marker_clean and marker_clean.lower() in text_lower:
            out.append(marker_clean)
    return sorted(set(out), key=str.lower)


def _batch_index_from_path(batch_path: Path, payload: dict[str, Any]) -> int | None:
    value = payload.get("batch_index")
    try:
        return int(value)
    except (TypeError, ValueError):
        match = re.search(r"batch_(\d+)\.json$", batch_path.name)
        if match:
            return int(match.group(1))
    return None


def _find_bad_batches(input_dir: Path, species: str, markers: list[str]) -> list[BadBatch]:
    batch_root = input_dir / ".batch_state"
    if not batch_root.is_dir():
        return []

    bad_batches: list[BadBatch] = []
    for state_dir in sorted(batch_root.glob(f"{species}_*")):
        if not state_dir.is_dir():
            continue
        state_path = state_dir / "state.json"
        for batch_path in sorted(state_dir.glob("batch_*.json")):
            payload = _read_json(batch_path, default={})
            if not isinstance(payload, dict):
                continue
            batch_index = _batch_index_from_path(batch_path, payload)
            if batch_index is None:
                print(f"[warn] Could not resolve batch index for {batch_path}", file=sys.stderr)
                continue

            bad_gene_ids: set[str] = set()
            reasons: list[str] = []
            for record in payload.get("records", []) or []:
                if not isinstance(record, dict):
                    continue
                record_reasons = _record_bad_reasons(record, markers)
                if not record_reasons:
                    continue
                gene_id = _gene_id(record)
                if gene_id:
                    bad_gene_ids.add(gene_id)
                reasons.extend(record_reasons)

            if bad_gene_ids:
                bad_batches.append(
                    BadBatch(
                        input_dir=input_dir,
                        state_dir=state_dir,
                        state_path=state_path,
                        batch_path=batch_path,
                        batch_index=batch_index,
                        species=species,
                        bad_gene_ids=bad_gene_ids,
                        reasons=sorted(set(reasons), key=str.lower),
                    )
                )
    return bad_batches


def _run_store_path(input_dir: Path, species: str) -> Path:
    safe_species = re.sub(r"[^A-Za-z0-9_.-]+", "_", species.strip() or "fly")
    return input_dir / ".batch_state" / f"run_store_{safe_species}.json"


def _remove_state_batch(state_path: Path, batch_index: int) -> bool:
    state = _read_json(state_path, default={})
    if not isinstance(state, dict):
        print(f"[warn] Invalid state file format: {state_path}", file=sys.stderr)
        return False
    changed = False

    completed = state.get("completed_batches")
    if isinstance(completed, list):
        new_completed = []
        for value in completed:
            try:
                keep = int(value) != batch_index
            except (TypeError, ValueError):
                keep = True
            if keep:
                new_completed.append(value)
        if len(new_completed) != len(completed):
            state["completed_batches"] = new_completed
            changed = True

    batch_files = state.get("batch_files")
    if isinstance(batch_files, dict):
        key = str(batch_index)
        if key in batch_files:
            batch_files.pop(key, None)
            changed = True

    if changed:
        _atomic_json_dump(state_path, state)
    return changed


def _remove_run_store_genes(run_store_path: Path, gene_ids: set[str]) -> int:
    run_store = _read_json(run_store_path, default={})
    if not isinstance(run_store, dict):
        print(f"[warn] Invalid run-store format: {run_store_path}", file=sys.stderr)
        return 0
    genes = run_store.get("genes")
    if not isinstance(genes, dict):
        return 0

    keys_to_remove = []
    for key, record in genes.items():
        if not isinstance(record, dict):
            continue
        if _gene_id(record) in gene_ids:
            keys_to_remove.append(key)

    for key in keys_to_remove:
        genes.pop(key, None)
    if keys_to_remove:
        _atomic_json_dump(run_store_path, run_store)
    return len(keys_to_remove)


def _quarantine_batch(batch_path: Path, input_dir: Path, timestamp: str, *, delete: bool) -> Path | None:
    if not batch_path.exists():
        return None
    if delete:
        batch_path.unlink()
        return None
    batch_root = input_dir / ".batch_state"
    quarantine_dir = batch_root / ".purged_batches" / timestamp / batch_path.parent.name
    quarantine_dir.mkdir(parents=True, exist_ok=True)
    target_path = quarantine_dir / batch_path.name
    shutil.move(str(batch_path), str(target_path))
    return target_path


def _print_bad_batches(bad_batches: list[BadBatch]) -> None:
    if not bad_batches:
        print("No bad batches found.")
        return
    print(f"Found {len(bad_batches)} bad batch(es):")
    for item in bad_batches:
        print(f"- {item.batch_path}")
        print(f"  input_dir: {item.input_dir}")
        print(f"  batch_index: {item.batch_index}")
        print(f"  bad_gene_count: {len(item.bad_gene_ids)}")
        print(f"  markers: {', '.join(item.reasons)}")


def _apply_purge(bad_batches: list[BadBatch], *, delete_batches: bool) -> None:
    timestamp = _utc_stamp()

    batches_by_state: dict[Path, list[BadBatch]] = {}
    genes_by_input_species: dict[tuple[Path, str], set[str]] = {}
    for item in bad_batches:
        batches_by_state.setdefault(item.state_path, []).append(item)
        genes_by_input_species.setdefault((item.input_dir, item.species), set()).update(item.bad_gene_ids)

    for state_path in sorted(batches_by_state):
        backup_path = _backup(state_path, timestamp)
        if backup_path:
            print(f"Backed up state: {backup_path}")
        for item in batches_by_state[state_path]:
            state_changed = _remove_state_batch(state_path, item.batch_index)
            moved_to = _quarantine_batch(
                item.batch_path,
                item.input_dir,
                timestamp,
                delete=delete_batches,
            )
            action = "deleted" if delete_batches else f"moved to {moved_to}"
            print(
                f"Purged batch {item.batch_index} from {item.state_dir.name}: "
                f"state_changed={state_changed}, batch_file={action}"
            )

    for (input_dir, species), gene_ids in sorted(
        genes_by_input_species.items(), key=lambda item: (str(item[0][0]), item[0][1])
    ):
        run_store_path = _run_store_path(input_dir, species)
        backup_path = _backup(run_store_path, timestamp)
        if backup_path:
            print(f"Backed up run store: {backup_path}")
        removed = _remove_run_store_genes(run_store_path, gene_ids)
        print(f"Removed {removed} run-store record(s) from {run_store_path}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Find checkpoint batches contaminated by OpenAI quota/credit/API errors "
            "and purge them so a non-force rerun reprocesses affected genes."
        )
    )
    parser.add_argument(
        "--input-dir",
        action="append",
        default=[],
        help=(
            "Pipeline input directory to scan. Repeatable. Defaults to the four "
            "Tx-Omics-Conserved Large split folders."
        ),
    )
    parser.add_argument(
        "--species",
        default="human",
        help="Species key to scan in .batch_state (default: human).",
    )
    parser.add_argument(
        "--marker",
        action="append",
        default=[],
        help=(
            "Bad-output marker to search for. Repeatable. When omitted, built-in "
            "quota/credit/API-error markers are used."
        ),
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Apply changes. Without this flag, only prints a dry-run report.",
    )
    parser.add_argument(
        "--delete-batches",
        action="store_true",
        help=(
            "Delete bad batch files instead of moving them to "
            ".batch_state/.purged_batches/<timestamp>/."
        ),
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    input_dirs = [Path(path).expanduser().resolve() for path in args.input_dir]
    if not input_dirs:
        input_dirs = DEFAULT_INPUT_DIRS
    species = str(args.species or "human").strip() or "human"
    markers = [str(marker) for marker in args.marker if str(marker).strip()] or DEFAULT_MARKERS

    bad_batches: list[BadBatch] = []
    for input_dir in input_dirs:
        if not input_dir.is_dir():
            print(f"[warn] Input directory not found, skipping: {input_dir}", file=sys.stderr)
            continue
        bad_batches.extend(_find_bad_batches(input_dir, species, markers))

    print(f"Mode: {'APPLY' if args.apply else 'DRY RUN'}")
    print(f"Species: {species}")
    print(f"Markers: {', '.join(markers)}")
    print(f"Scanned input dirs: {len(input_dirs)}")
    _print_bad_batches(bad_batches)

    if not args.apply:
        print("Dry run only. Re-run with --apply to purge the listed batches.")
        return 0
    if not bad_batches:
        return 0

    _apply_purge(bad_batches, delete_batches=bool(args.delete_batches))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
