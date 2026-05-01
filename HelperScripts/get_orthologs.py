#!/usr/bin/env python3
"""Create fly-to-human or fly-to-mouse ortholog CSVs with fly traceability."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import re
import shutil
from pathlib import Path
from typing import Any, Mapping

import pandas as pd

from HelperScripts.GetFBgnIDs import DEFAULT_FLYBASE_DATA, find_latest_tsv, load_flybase_tsv
from HelperScripts.diopt_client import DIOPT_TAXON_HUMAN, DIOPT_TAXON_MOUSE, DIOPTClient
from HelperScripts.species_data_utils import default_data_root, download_if_needed, normalize_columns


NCBI_GENE_INFO_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz"
DEFAULT_DIOPT_WORKERS = 8


TARGETS = {
    "human": {
        "taxon": DIOPT_TAXON_HUMAN,
        "csv_suffix": "_human.csv",
        "output_dir": "Human",
        "entrez_col": "human_entrez_gene_id",
        "symbol_col": "human_gene_symbol",
        "authority_col": "hgnc_id",
    },
    "mouse": {
        "taxon": DIOPT_TAXON_MOUSE,
        "csv_suffix": "_mouse.csv",
        "output_dir": "Mouse",
        "entrez_col": "mouse_entrez_gene_id",
        "symbol_col": "mouse_gene_symbol",
        "authority_col": "mgi_id",
    },
}


def _clean(value: Any) -> str:
    return str(value or "").strip()


def load_fbgn_to_primary_symbol(flybase_data_dir: Path) -> dict[str, str]:
    synonym_path = find_latest_tsv(flybase_data_dir / "Genes", "fb_synonym")
    df = load_flybase_tsv(synonym_path, keep_default_na=False)
    df = df[df.organism_abbreviation == "Dmel"]
    df = df[df["primary_FBid"].astype(str).str.startswith("FBgn", na=False)]
    fbgns = df["primary_FBid"].astype(str).str.strip()
    symbols = df["current_symbol"].astype(str).str.strip()
    non_empty = fbgns.ne("")
    return dict(zip(fbgns[non_empty], symbols[non_empty]))


def load_fbgn_to_entrez(cache_dir: Path | None = None, refresh: bool = False) -> dict[str, str]:
    """Load FlyBase FBgn -> NCBI Entrez Gene ID from NCBI gene_info."""
    cache_root = Path(cache_dir) if cache_dir else default_data_root() / "diopt"
    path = download_if_needed(NCBI_GENE_INFO_URL, cache_root / "gene_info.gz", refresh=refresh)
    mapping: dict[str, str] = {}
    wanted_cols = {"tax_id", "geneid", "dbxrefs"}

    def wanted_col(column: str) -> bool:
        normalized = re.sub(r"\s+", "_", str(column).strip().lstrip("#").lower())
        return normalized in wanted_cols

    for chunk in pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        chunksize=500_000,
        usecols=wanted_col,
    ):
        chunk = normalize_columns(chunk)
        if not {"tax_id", "geneid", "dbxrefs"}.issubset(set(chunk.columns)):
            continue
        fly = chunk.loc[chunk["tax_id"].eq("7227"), ["geneid", "dbxrefs"]]
        if fly.empty:
            continue
        fbgn_matches = fly["dbxrefs"].astype(str).str.extractall(r"(FBgn\d+)")[0]
        if fbgn_matches.empty:
            continue
        match_rows = fbgn_matches.index.get_level_values(0)
        entrez_ids = fly.loc[match_rows, "geneid"].astype(str).str.strip().to_numpy()
        fbgn_values = fbgn_matches.astype(str).to_numpy()
        mapping.update(
            {
                fbgn: entrez
                for fbgn, entrez in zip(fbgn_values, entrez_ids)
                if entrez
            }
        )
    return mapping


def _csv_files(input_directory: Path) -> list[Path]:
    out = []
    for path in sorted(input_directory.glob("*.csv")):
        name = path.name
        if name.endswith("_human.csv") or name.endswith("_mouse.csv"):
            continue
        if name.endswith("_classification.csv"):
            continue
        out.append(path)
    return out


def _output_paths(csv_path: Path, target: Mapping[str, Any]) -> tuple[Path, Path]:
    out_dir = csv_path.parent / str(target["output_dir"])
    out_path = out_dir / f"{csv_path.stem}{target['csv_suffix']}"
    legacy_out_path = csv_path.with_name(f"{csv_path.stem}{target['csv_suffix']}")
    return out_path, legacy_out_path


def _input_records(df: pd.DataFrame, input_gene_col: str) -> list[dict[str, str]]:
    fbgns = df["flybase_gene_id"].astype(str).str.strip()
    filtered = df.loc[fbgns.str.startswith("FBgn", na=False)].copy()
    if filtered.empty:
        return []
    filtered["_fbgn"] = filtered["flybase_gene_id"].astype(str).str.strip()
    if input_gene_col in filtered.columns:
        filtered["_input_symbol"] = filtered[input_gene_col].astype(str).str.strip()
    else:
        filtered["_input_symbol"] = ""
    filtered = filtered.drop_duplicates(subset="_fbgn", keep="first")
    return filtered[["_fbgn", "_input_symbol"]].to_dict("records")


def _coerce_worker_count(diopt_workers: int | None, lookup_count: int) -> int:
    if lookup_count <= 0:
        return 1
    if diopt_workers is None:
        diopt_workers = DEFAULT_DIOPT_WORKERS
    return max(1, min(int(diopt_workers), lookup_count))


def _fetch_orthologs_for_records(
    records: list[dict[str, str]],
    *,
    client: DIOPTClient,
    fbgn_to_entrez: Mapping[str, str],
    output_taxon: int,
    diopt_filter: str,
    diopt_workers: int | None,
) -> dict[str, list[dict[str, Any]]]:
    worker_count = _coerce_worker_count(diopt_workers, len(records))

    def fetch(record: dict[str, str]) -> tuple[str, list[dict[str, Any]]]:
        fbgn = record["_fbgn"]
        found = client.orthologs(
            fbgn,
            output_taxon=output_taxon,
            diopt_filter=diopt_filter,
            query_id=fbgn_to_entrez.get(fbgn, fbgn),
        )
        return fbgn, found

    if worker_count == 1:
        return dict(fetch(record) for record in records)

    with ThreadPoolExecutor(max_workers=worker_count) as pool:
        return dict(pool.map(fetch, records))


def convert_csv(
    csv_path: Path,
    *,
    orthologs: str,
    input_gene_col: str,
    diopt_filter: str,
    flybase_data_dir: Path,
    cache_dir: Path | None = None,
    reuse_existing: bool = True,
    fbgn_to_symbol: Mapping[str, str] | None = None,
    fbgn_to_entrez: Mapping[str, str] | None = None,
    diopt_client: DIOPTClient | None = None,
    diopt_workers: int | None = DEFAULT_DIOPT_WORKERS,
) -> Path:
    if orthologs not in TARGETS:
        raise ValueError(f"Unsupported ortholog target: {orthologs}")

    target = TARGETS[orthologs]
    out_path, legacy_out_path = _output_paths(csv_path, target)
    if reuse_existing and out_path.exists():
        print(f"Using existing {orthologs} ortholog CSV: {out_path}")
        return out_path
    if reuse_existing and legacy_out_path.exists():
        out_dir = out_path.parent
        out_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(legacy_out_path, out_path)
        print(f"Copied existing {orthologs} ortholog CSV to: {out_path}")
        return out_path

    if fbgn_to_symbol is None:
        fbgn_to_symbol = load_fbgn_to_primary_symbol(flybase_data_dir)
    if fbgn_to_entrez is None:
        fbgn_to_entrez = load_fbgn_to_entrez(cache_dir=cache_dir)
    if diopt_client is None:
        diopt_client = DIOPTClient(cache_dir=cache_dir)
    df = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
    if "flybase_gene_id" not in df.columns:
        raise ValueError(f"{csv_path} does not contain flybase_gene_id")

    rows: list[dict[str, str]] = []
    records = _input_records(df, input_gene_col)
    orthologs_by_fbgn = _fetch_orthologs_for_records(
        records,
        client=diopt_client,
        fbgn_to_entrez=fbgn_to_entrez,
        output_taxon=int(target["taxon"]),
        diopt_filter=diopt_filter,
        diopt_workers=diopt_workers,
    )
    for input_row in records:
        fbgn = _clean(input_row.get("_fbgn", ""))
        input_symbol = _clean(input_row.get("_input_symbol", ""))
        primary_symbol = fbgn_to_symbol.get(fbgn, input_symbol)
        trace = {
            "fly_gene_symbol_input": input_symbol,
            "fly_gene_symbol_primary": primary_symbol,
            "flybase_gene_id": fbgn,
        }

        found = orthologs_by_fbgn.get(fbgn, [])
        if not found:
            rows.append({
                **trace,
                target["symbol_col"]: "",
                target["authority_col"]: "",
                target["entrez_col"]: "",
                "diopt_score": "",
                "diopt_best_match": "",
                "diopt_supporting_algorithms": "",
                "status": "unmapped",
            })
            continue

        for ortholog in found:
            rows.append({
                **trace,
                target["symbol_col"]: _clean(ortholog.get("gene_symbol", "")),
                target["authority_col"]: _clean(ortholog.get("authority_id", "")),
                target["entrez_col"]: _clean(ortholog.get("entrez_gene_id", "")),
                "diopt_score": _clean(ortholog.get("diopt_score", "")),
                "diopt_best_match": _clean(ortholog.get("diopt_best_match", "")),
                "diopt_supporting_algorithms": _clean(ortholog.get("diopt_supporting_algorithms", "")),
                "status": "mapped",
            })

    columns = [
        "fly_gene_symbol_input",
        "fly_gene_symbol_primary",
        "flybase_gene_id",
        target["symbol_col"],
        target["authority_col"],
        target["entrez_col"],
        "diopt_score",
        "diopt_best_match",
        "diopt_supporting_algorithms",
        "status",
    ]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=columns).to_csv(out_path, index=False, encoding="utf-8-sig")
    return out_path


def convert_directory(
    input_directory: str | Path,
    *,
    orthologs: str,
    input_gene_col: str = "ext_gene",
    diopt_filter: str = "exclude_low_score_2",
    flybase_data_dir: str | Path = DEFAULT_FLYBASE_DATA,
    cache_dir: str | Path | None = None,
    reuse_existing: bool = True,
    diopt_workers: int | None = DEFAULT_DIOPT_WORKERS,
) -> list[Path]:
    input_dir = Path(input_directory)
    flybase_dir = Path(flybase_data_dir)
    cache_path = Path(cache_dir) if cache_dir else None
    target = TARGETS[orthologs]
    csv_paths = _csv_files(input_dir)
    needs_conversion = [
        csv_path
        for csv_path in csv_paths
        if not (
            reuse_existing
            and (
                _output_paths(csv_path, target)[0].exists()
                or _output_paths(csv_path, target)[1].exists()
            )
        )
    ]
    fbgn_to_symbol = load_fbgn_to_primary_symbol(flybase_dir) if needs_conversion else None
    fbgn_to_entrez = load_fbgn_to_entrez(cache_dir=cache_path) if needs_conversion else None
    diopt_client = DIOPTClient(cache_dir=cache_path) if needs_conversion else None
    outputs = []
    for csv_path in csv_paths:
        outputs.append(
            convert_csv(
                csv_path,
                orthologs=orthologs,
                input_gene_col=input_gene_col,
                diopt_filter=diopt_filter,
                flybase_data_dir=flybase_dir,
                cache_dir=cache_path,
                reuse_existing=reuse_existing,
                fbgn_to_symbol=fbgn_to_symbol,
                fbgn_to_entrez=fbgn_to_entrez,
                diopt_client=diopt_client,
                diopt_workers=diopt_workers,
            )
        )
    return outputs


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Create fly-to-human or fly-to-mouse ortholog CSVs.")
    parser.add_argument("input_directory")
    parser.add_argument("--orthologs", choices=sorted(TARGETS), required=True)
    parser.add_argument("--input-gene-col", default="ext_gene")
    parser.add_argument("--diopt-filter", default="exclude_low_score_2")
    parser.add_argument("--flybase-data-dir", default=str(DEFAULT_FLYBASE_DATA))
    parser.add_argument("--diopt-cache-dir", default="")
    parser.add_argument("--diopt-workers", type=int, default=DEFAULT_DIOPT_WORKERS)
    args = parser.parse_args(argv)

    outputs = convert_directory(
        args.input_directory,
        orthologs=args.orthologs,
        input_gene_col=args.input_gene_col,
        diopt_filter=args.diopt_filter,
        flybase_data_dir=args.flybase_data_dir,
        cache_dir=args.diopt_cache_dir or None,
        diopt_workers=args.diopt_workers,
    )
    for path in outputs:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
