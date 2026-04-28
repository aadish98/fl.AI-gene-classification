#!/usr/bin/env python3
"""Create fly-to-human or fly-to-mouse ortholog CSVs with fly traceability."""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path
from typing import Any

import pandas as pd

from HelperScripts.GetFBgnIDs import DEFAULT_FLYBASE_DATA, find_latest_tsv, load_flybase_tsv
from HelperScripts.diopt_client import DIOPT_TAXON_HUMAN, DIOPT_TAXON_MOUSE, fbgn_to_orthologs
from HelperScripts.species_data_utils import default_data_root, download_if_needed, normalize_columns


NCBI_GENE_INFO_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz"


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
    return {
        str(row["primary_FBid"]).strip(): str(row["current_symbol"]).strip()
        for _, row in df.iterrows()
        if str(row.get("primary_FBid", "")).strip()
    }


def load_fbgn_to_entrez(cache_dir: Path | None = None, refresh: bool = False) -> dict[str, str]:
    """Load FlyBase FBgn -> NCBI Entrez Gene ID from NCBI gene_info."""
    cache_root = Path(cache_dir) if cache_dir else default_data_root() / "diopt"
    path = download_if_needed(NCBI_GENE_INFO_URL, cache_root / "gene_info.gz", refresh=refresh)
    mapping: dict[str, str] = {}
    for chunk in pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, chunksize=500_000):
        chunk = normalize_columns(chunk)
        if not {"tax_id", "geneid", "dbxrefs"}.issubset(set(chunk.columns)):
            continue
        fly = chunk[chunk["tax_id"] == "7227"]
        for _, row in fly.iterrows():
            entrez = str(row.get("geneid", "") or "").strip()
            dbxrefs = str(row.get("dbxrefs", "") or "")
            for match in re.findall(r"FBgn\d+", dbxrefs):
                if entrez:
                    mapping[match] = entrez
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


def convert_csv(
    csv_path: Path,
    *,
    orthologs: str,
    input_gene_col: str,
    diopt_filter: str,
    flybase_data_dir: Path,
    cache_dir: Path | None = None,
    reuse_existing: bool = True,
) -> Path:
    if orthologs not in TARGETS:
        raise ValueError(f"Unsupported ortholog target: {orthologs}")

    target = TARGETS[orthologs]
    out_dir = csv_path.parent / str(target["output_dir"])
    out_path = out_dir / f"{csv_path.stem}{target['csv_suffix']}"
    if reuse_existing and out_path.exists():
        print(f"Using existing {orthologs} ortholog CSV: {out_path}")
        return out_path
    legacy_out_path = csv_path.with_name(f"{csv_path.stem}{target['csv_suffix']}")
    if reuse_existing and legacy_out_path.exists():
        out_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(legacy_out_path, out_path)
        print(f"Copied existing {orthologs} ortholog CSV to: {out_path}")
        return out_path

    fbgn_to_symbol = load_fbgn_to_primary_symbol(flybase_data_dir)
    fbgn_to_entrez = load_fbgn_to_entrez(cache_dir=cache_dir)
    df = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
    if "flybase_gene_id" not in df.columns:
        raise ValueError(f"{csv_path} does not contain flybase_gene_id")

    rows: list[dict[str, str]] = []
    seen_input_fbgns = set()
    for _, input_row in df.iterrows():
        fbgn = _clean(input_row.get("flybase_gene_id", ""))
        if not fbgn.startswith("FBgn"):
            continue
        input_symbol = _clean(input_row.get(input_gene_col, ""))
        primary_symbol = fbgn_to_symbol.get(fbgn, input_symbol)
        trace = {
            "fly_gene_symbol_input": input_symbol,
            "fly_gene_symbol_primary": primary_symbol,
            "flybase_gene_id": fbgn,
        }
        if fbgn in seen_input_fbgns:
            continue
        seen_input_fbgns.add(fbgn)

        found = fbgn_to_orthologs(
            fbgn,
            output_taxon=int(target["taxon"]),
            diopt_filter=diopt_filter,
            query_id=fbgn_to_entrez.get(fbgn, fbgn),
            cache_dir=cache_dir,
        )
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
    out_dir.mkdir(parents=True, exist_ok=True)
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
) -> list[Path]:
    input_dir = Path(input_directory)
    flybase_dir = Path(flybase_data_dir)
    cache_path = Path(cache_dir) if cache_dir else None
    outputs = []
    for csv_path in _csv_files(input_dir):
        outputs.append(
            convert_csv(
                csv_path,
                orthologs=orthologs,
                input_gene_col=input_gene_col,
                diopt_filter=diopt_filter,
                flybase_data_dir=flybase_dir,
                cache_dir=cache_path,
                reuse_existing=reuse_existing,
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
    args = parser.parse_args(argv)

    outputs = convert_directory(
        args.input_directory,
        orthologs=args.orthologs,
        input_gene_col=args.input_gene_col,
        diopt_filter=args.diopt_filter,
        flybase_data_dir=args.flybase_data_dir,
        cache_dir=args.diopt_cache_dir or None,
    )
    for path in outputs:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
