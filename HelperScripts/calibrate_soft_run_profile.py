"""Rebuild the soft-run token-average profile from saved run-store telemetry.

Scans Data/**/.batch_state/run_store*.json, pools the per-request OpenAI usage
recorded in each gene's `openai_usage.requests`, and writes the averages that
`flai_gene_classification._estimate_costs` reads into
.flai_system/soft_run_reference_profiles.json.

If telemetry exists for the currently configured (summary, classification) model
pair, a profile keyed on that pair is written from those events only. Otherwise
(or with --default) all telemetry is pooled into a model-independent default
profile. Species keyword pass rates in the file are left untouched.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import flai_gene_classification as flai  # noqa: E402

SUMMARY_REQUEST_TYPES = {"summary_chunk", "summary_rewrite"}
CLASSIFICATION_REQUEST_TYPE = "classification"
MIN_REFERENCES_PER_STORE = 50


def _iter_run_stores(data_dir: Path):
    for path in sorted(data_dir.glob("**/.batch_state/run_store*.json")):
        if path.suffix == ".json" and ".lock" not in path.name:
            yield path


def _collect(path: Path) -> dict[tuple[str, str], Counter]:
    """Return token counters keyed by (summary_model, classification_model) for one run store."""
    store = flai._json_load(path, default=None)
    if not isinstance(store, dict):
        return {}
    events: list[tuple[str, str, str, dict, str]] = []
    for gene_id, record in (store.get("genes") or {}).items():
        requests = ((record or {}).get("openai_usage") or {}).get("requests") or []
        for req in requests:
            events.append((
                str(req.get("request_type", "")),
                str(req.get("model", "")),
                str(gene_id),
                req.get("usage") or {},
                str(req.get("paper_id", "")),
            ))
    summary_models = {m for t, m, _, _, _ in events if t in SUMMARY_REQUEST_TYPES}
    classification_models = {m for t, m, _, _, _ in events if t == CLASSIFICATION_REQUEST_TYPE}
    if len(summary_models) != 1 or len(classification_models) != 1:
        pair = ("", "")
    else:
        pair = (next(iter(summary_models)), next(iter(classification_models)))

    counts: Counter = Counter()
    refs: set[tuple[str, str]] = set()
    classified_genes: set[str] = set()
    for req_type, _model, gene_id, usage, paper_id in events:
        tokens_in = flai._to_float(usage.get("input_tokens"))
        tokens_out = flai._to_float(usage.get("output_tokens"))
        if req_type in SUMMARY_REQUEST_TYPES:
            counts["summary_requests"] += 1
            counts["summary_input"] += tokens_in
            counts["summary_output"] += tokens_out
            refs.add((gene_id, paper_id))
        elif req_type == CLASSIFICATION_REQUEST_TYPE:
            counts["classification_requests"] += 1
            counts["classification_input"] += tokens_in
            counts["classification_output"] += tokens_out
            classified_genes.add(gene_id)
    counts["references"] = len(refs)
    counts["classified_genes"] = len(classified_genes)
    if counts["references"] < MIN_REFERENCES_PER_STORE:
        return {}
    return {pair: counts}


def _build_profile(counts: Counter, sources: list[str], label: str) -> dict[str, object]:
    refs = max(counts["references"], 1)
    genes = max(counts["classified_genes"], 1)
    return {
        "label": label,
        "summary_requests_per_reference": counts["summary_requests"] / refs,
        "summary_input_tokens_per_reference": counts["summary_input"] / refs,
        "summary_output_tokens_per_reference": counts["summary_output"] / refs,
        "classification_requests_per_gene": counts["classification_requests"] / genes,
        "classification_input_tokens_per_gene": counts["classification_input"] / genes,
        "classification_output_tokens_per_gene": counts["classification_output"] / genes,
        "observed_references": int(counts["references"]),
        "observed_gpt_classifications": int(counts["classified_genes"]),
        "observed_model_requests": int(counts["summary_requests"] + counts["classification_requests"]),
        "source_runs": sources,
        "source": "saved run stores under Data/**/.batch_state",
        "calibrated_at_utc": flai._utc_now_iso(),
        "method": "pooled per-request telemetry",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--data-dir", default=str(PROJECT_ROOT / "Data"), help="Directory to scan for run stores.")
    parser.add_argument("--default", action="store_true", help="Write the model-independent default profile from all telemetry.")
    parser.add_argument("--dry-run", action="store_true", help="Print the result without writing the profile file.")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    by_pair: dict[tuple[str, str], Counter] = defaultdict(Counter)
    sources_by_pair: dict[tuple[str, str], list[str]] = defaultdict(list)
    pooled: Counter = Counter()
    pooled_sources: list[str] = []
    for path in _iter_run_stores(data_dir):
        for pair, counts in _collect(path).items():
            run_dir = str(path.parent.parent.relative_to(data_dir)) if data_dir in path.parents else str(path.parent.parent)
            by_pair[pair].update(counts)
            sources_by_pair[pair].append(run_dir)
            pooled.update(counts)
            pooled_sources.append(run_dir)

    if not pooled:
        print(f"No run-store telemetry with at least {MIN_REFERENCES_PER_STORE} references found under {data_dir}.")
        return 1

    current_pair = (flai._get_summary_model(), flai._get_classification_model())
    if not args.default and current_pair in by_pair:
        profile = _build_profile(
            by_pair[current_pair],
            sources_by_pair[current_pair],
            f"Calibrated on {current_pair[0]} + {current_pair[1]} run telemetry",
        )
        profile["summary_model"], profile["classification_model"] = current_pair
    else:
        profile = _build_profile(
            pooled,
            sorted(set(pooled_sources)),
            "Model-independent default profile (pooled run telemetry)",
        )
        profile["default"] = True
        profile["calibrated_on"] = "previous default models (see Docs/model_history.md)"

    print(f"Run stores used: {len(pooled_sources)}")
    for key in (
        "observed_references", "observed_gpt_classifications",
        "summary_requests_per_reference", "summary_input_tokens_per_reference",
        "summary_output_tokens_per_reference", "classification_requests_per_gene",
        "classification_input_tokens_per_gene", "classification_output_tokens_per_gene",
    ):
        value = profile[key]
        print(f"  {key}: {value:,.2f}" if isinstance(value, float) else f"  {key}: {value:,}")

    if args.dry_run:
        print("Dry run: profile file not written.")
        return 0

    payload = flai._load_soft_run_system_profiles()
    # Keep a measured profile for the current pair unless this run replaces it; drop profiles
    # for any other (retired) pair and any previous default profile.
    kept = [
        p for p in payload.get("profiles", [])
        if isinstance(p, dict) and not p.get("default") and "default" in profile
        and (p.get("summary_model"), p.get("classification_model")) == current_pair
    ]
    payload["profiles"] = kept + [profile]
    flai._save_soft_run_system_profiles(payload)
    print(f"Wrote {flai.SOFT_RUN_REFERENCE_PROFILES_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
