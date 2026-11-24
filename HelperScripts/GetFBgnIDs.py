#!/usr/bin/env python3
"""
Map input gene symbols to FlyBase FBgn IDs in CSV files.

Each symbol is looked up exactly as written, in this order:
  1. the current symbol in the FlyBase synonym table;
  2. a symbol synonym in that table;
  3. if several genes share the synonym, the gene whose symbol FlyBase publications
     use most (HelperScripts/flybase_synonym_evidence.py), only when the lead is clear.
Greek letters are spelled out (zeta, alpha, ...) in a lookup key only when nothing
matched, because the FlyBase table spells them out. The input column is never altered
and nothing here calls an external API or an LLM.

Usage:
    python HelperScripts/GetFBgnIDs.py <input_directory> [input_gene_col]
"""

import argparse
import gzip
import os
import sys
import uuid
from dataclasses import dataclass
from glob import glob
from pathlib import Path
from typing import Optional

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from HelperScripts.build_synonym_evidence import EVIDENCE_PATH
from HelperScripts.flybase_data import ensure_flybase_data_files, resolve_flybase_data_dir
from HelperScripts.flybase_synonym_evidence import (
    SynonymEvidence,
    describe_candidates,
    describe_counts,
    load_evidence,
    pick_winner,
)

symbol_to_name = {
    'α': 'alpha',
    'β': 'beta',
    'γ': 'gamma',
    'δ': 'delta',  # Lowercase delta
    'ε': 'epsilon',
    'ζ': 'zeta',
    'η': 'eta',
    'θ': 'theta',
    'ι': 'iota',
    'κ': 'kappa',
    'λ': 'lambda',
    'μ': 'mu',
    'ν': 'nu',
    'ξ': 'xi',
    'ο': 'omicron',
    'π': 'pi',
    'ρ': 'rho',
    'σ': 'sigma',
    'τ': 'tau',
    'υ': 'upsilon',
    'φ': 'phi',
    'χ': 'chi',
    'ψ': 'psi',
    'ω': 'omega',
    'Α': 'Alpha',  # Uppercase Alpha
    'Β': 'Beta',
    'Γ': 'Gamma',
    'Δ': 'Delta',  # Uppercase Delta
    'Ε': 'Epsilon',
    'Ζ': 'Zeta',
    'Η': 'Eta',
    'Θ': 'Theta',
    'Ι': 'Iota',
    'Κ': 'Kappa',
    'Λ': 'Lambda',
    'Μ': 'Mu',
    'Ν': 'Nu',
    'Ξ': 'Xi',
    'Ο': 'Omicron',
    'Π': 'Pi',
    'Ρ': 'Rho',
    'Σ': 'Sigma',
    'Τ': 'Tau',
    'Υ': 'Upsilon',
    'Φ': 'Phi',
    'Χ': 'Chi',
    'Ψ': 'Psi',
    'Ω': 'Omega'
}


DEFAULT_FLYBASE_DATA = resolve_flybase_data_dir()


def atomic_write_csv(df: pd.DataFrame, csv_file: str | Path, **to_csv_kwargs):
    """Write a CSV via a unique temp file, then atomically replace the target."""
    path = Path(csv_file)
    tmp_path = path.with_name(f".{path.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp")
    try:
        df.to_csv(tmp_path, **to_csv_kwargs)
        os.replace(tmp_path, path)
    finally:
        tmp_path.unlink(missing_ok=True)


def replace_symbol(gene_series):
    for symbol, name in symbol_to_name.items():
        gene_series = gene_series.str.replace(symbol, name, regex=False)
    return gene_series


@dataclass(frozen=True)
class FlyGene:
    fbgn: str
    current: str
    synonyms: frozenset


_MISSING_VALUES = {"", "nan", "none", "-"}
GREEK_NOTE = "Greek letters spelled out for lookup"


class FlyBaseSymbols:
    """Exact-match symbol lookup: current symbol, then synonyms, then publication evidence."""

    def __init__(
        self,
        genes: dict[str, FlyGene],
        by_current: dict[str, str],
        by_synonym: dict[str, set[str]],
        evidence: Optional[SynonymEvidence] = None,
    ):
        self.genes = genes
        self.by_current = by_current
        self.by_synonym = by_synonym
        self.evidence = evidence

    def lookup(self, symbol: str) -> tuple[str, str]:
        """Return (FBgn, note). FBgn is '' when unmatched or ambiguous; ambiguity is explained in the note."""
        symbol = str(symbol or "").strip()
        if symbol.lower() in _MISSING_VALUES:
            return "", ""
        fbgn, note = self._lookup_exact(symbol)
        if fbgn or note:
            return fbgn, note
        spelled = replace_symbol(pd.Series([symbol])).iloc[0]
        if spelled != symbol:
            fbgn, note = self._lookup_exact(spelled)
            if fbgn or note:
                return fbgn, f"{GREEK_NOTE}; {note}" if note else GREEK_NOTE
        return "", ""

    def _lookup_exact(self, symbol: str) -> tuple[str, str]:
        fbgn = self.by_current.get(symbol)
        if fbgn:
            return fbgn, ""
        hits = sorted(self.by_synonym.get(symbol, ()))
        if len(hits) == 1:
            return hits[0], f"matched via synonym of {self.genes[hits[0]].current}"
        if len(hits) > 1:
            return self._resolve_shared_synonym(symbol, hits)
        return "", ""

    def _resolve_shared_synonym(self, symbol: str, hits: list[str]) -> tuple[str, str]:
        tsv_listing = ", ".join(f"{self.genes[fbgn].current} {fbgn}" for fbgn in hits)
        if self.evidence is None:
            return "", (
                "ambiguous synonym (evidence index missing; run HelperScripts/build_synonym_evidence.py): "
                f"{tsv_listing}"
            )
        candidates = self.evidence.candidates(symbol)
        if not candidates:
            return "", f"ambiguous synonym (no FlyBase citation data): {tsv_listing}"
        winner = pick_winner(candidates)
        if winner is None:
            return "", f"ambiguous synonym: {describe_candidates(candidates)}"
        return winner.fbgn, (
            f"synonym shared by {len(candidates)} genes; FlyBase citations: "
            f"{describe_counts(candidates)} \u2192 {winner.current}"
        )


def map_gene_ids(df: pd.DataFrame, symbols: FlyBaseSymbols, gene_col: str) -> pd.DataFrame:
    """Add flybase_gene_id and mapping_note. The gene column itself is never modified."""
    inputs = df[gene_col].astype(str)
    results = {symbol: symbols.lookup(symbol) for symbol in inputs.unique()}
    df["flybase_gene_id"] = inputs.map(lambda value: results[value][0])
    df["mapping_note"] = inputs.map(lambda value: results[value][1])

    notes = df["mapping_note"].astype(str)
    unmatched = df["flybase_gene_id"].eq("")
    print(f"TOTAL GENES: {len(df)}")
    print(f"Mapped: {int((~unmatched).sum())}")
    print(f"  via synonym: {int(notes.str.startswith('matched via synonym').sum())}")
    print(f"  shared synonym resolved by FlyBase citations: {int(notes.str.startswith('synonym shared').sum())}")
    print(f"Ambiguous (left unmapped): {int(notes.str.contains('ambiguous').sum())}")
    print(f"Unmatched: {int((unmatched & notes.eq('')).sum())}")
    shown = df.loc[unmatched, [gene_col, "mapping_note"]].drop_duplicates()
    if len(shown):
        print(shown.to_string(index=False))
    return df


def find_latest_tsv(directory: Path, pattern: str) -> Path:
    """Find latest FlyBase TSV, preferring .tsv.gz over .tsv."""
    gz_files = sorted(glob(str(directory / f"{pattern}*.tsv.gz")), reverse=True)
    if gz_files:
        return Path(gz_files[0])

    tsv_files = sorted(glob(str(directory / f"{pattern}*.tsv")), reverse=True)
    if tsv_files:
        return Path(tsv_files[0])

    raise FileNotFoundError(f"No TSV file matching '{pattern}' found in {directory}")


def load_flybase_tsv(filepath, **kwargs) -> pd.DataFrame:
    """
    Load a FlyBase TSV file, handling metadata/comments and header variants.
    Supports both .tsv and .tsv.gz files.
    """
    filepath = Path(filepath)

    if filepath.suffix == ".gz" or str(filepath).endswith(".tsv.gz"):
        opener = lambda f: gzip.open(f, "rt", encoding="utf-8")
    else:
        opener = lambda f: open(f, "r", encoding="utf-8")

    skip_rows = 0
    header_line = None
    use_custom_header = False
    header_found = False

    with opener(filepath) as f:
        for i, line in enumerate(f):
            stripped = line.strip()

            if not stripped:
                skip_rows = i + 1
                continue

            if header_found:
                if stripped.startswith("#"):
                    skip_rows = i + 1
                    continue
                skip_rows = i
                break

            if stripped.startswith("#"):
                if "\t" not in stripped:
                    skip_rows = i + 1
                    continue

                # Some FlyBase files use "##col1\tcol2..." headers.
                header_content = stripped.lstrip("#")
                cols = header_content.split("\t")
                col_names = [c.strip() for c in cols if c.strip()]

                if len(col_names) >= 2:
                    header_line = col_names
                    use_custom_header = True
                    header_found = True
                    skip_rows = i + 1
                    continue
                else:
                    skip_rows = i + 1
                    continue
            else:
                cols = stripped.split("\t")
                first_field = cols[0].strip() if cols else ""
                if first_field.startswith("FB") and len(first_field) > 4:
                    break
                else:
                    break

    if use_custom_header and header_line:
        df = pd.read_csv(
            filepath,
            sep="\t",
            skiprows=skip_rows,
            names=header_line,
            low_memory=False,
            on_bad_lines="warn",
            **kwargs,
        )
    else:
        df = pd.read_csv(
            filepath,
            sep="\t",
            skiprows=skip_rows,
            low_memory=False,
            on_bad_lines="warn",
            **kwargs,
        )

    return df


def _split_pipe(value: str) -> list[str]:
    return [part.strip() for part in str(value or "").split("|") if part.strip() and part.strip() != "-"]


def load_symbol_index(flybase_data_dir: Path, evidence: Optional[SynonymEvidence] = None) -> FlyBaseSymbols:
    columns_as_strings = {
        "current_symbol": str,
        "symbol_synonym(s)": str,
        "primary_FBid": str,
    }
    synonym_path = find_latest_tsv(flybase_data_dir / "Genes", "fb_synonym")
    print(f"Loading FlyBase synonym table: {synonym_path}")
    mappings_df = load_flybase_tsv(synonym_path, dtype=columns_as_strings, keep_default_na=False)
    mappings_df = mappings_df[mappings_df.organism_abbreviation == "Dmel"]
    mappings_df = mappings_df[mappings_df["primary_FBid"].str.startswith("FBgn", na=False)]

    genes: dict[str, FlyGene] = {}
    by_current: dict[str, str] = {}
    by_synonym: dict[str, set[str]] = {}
    for fbgn, current, symbol_synonyms in mappings_df[["primary_FBid", "current_symbol", "symbol_synonym(s)"]].itertuples(index=False, name=None):
        fbgn = str(fbgn).strip()
        current = str(current).strip()
        synonyms = frozenset(_split_pipe(symbol_synonyms))
        genes[fbgn] = FlyGene(fbgn, current, synonyms)
        if current:
            by_current.setdefault(current, fbgn)
        for synonym in synonyms:
            by_synonym.setdefault(synonym, set()).add(fbgn)
    return FlyBaseSymbols(genes, by_current, by_synonym, evidence)


def process_csv_file(csv_file, gene_column, symbols: FlyBaseSymbols):
    print(f"Processing file: {csv_file}")

    header_row = 0
    try:
        temp_df = pd.read_csv(csv_file, nrows=10, header=None, encoding="utf-8-sig")
        for i, row in temp_df.iterrows():
            if gene_column in row.astype(str).values:
                header_row = i
                print(f"Detected header at row {i}")
                break
    except Exception as e:
        print(f"Header detection skipped/failed (using default 0): {e}")

    custom_na_values = [
        "",
        "#N/A",
        "#N/A N/A",
        "#NA",
        "-1.#IND",
        "-1.#QNAN",
        "-NaN",
        "-nan",
        "1.#IND",
        "1.#QNAN",
        "<NA>",
        "N/A",
        "NA",
        "NULL",
        "NaN",
        "n/a",
        "nan",
        "null",
    ]
    df = pd.read_csv(
        csv_file,
        encoding="utf-8-sig",
        header=header_row,
        na_values=custom_na_values,
        keep_default_na=False,
        low_memory=False,
    )

    if gene_column not in df.columns:
        if "GO" in csv_file and "Genes" in df.columns:
            df["Genes"] = df["Genes"].astype(str).str.split(", ")
            df = df.explode("Genes").reset_index(drop=True)
            df[gene_column] = df["Genes"]
        else:
            print(f"'{gene_column}' column not found in {csv_file}. Skipping file.")
            return

    for generated in ("flybase_gene_id", "mapping_note"):
        if generated in df.columns:
            df.drop(columns=[generated], inplace=True)

    df = map_gene_ids(df, symbols, gene_column)
    df["flybase_gene_id"] = df["flybase_gene_id"].replace("", "-").fillna("-")
    df["mapping_note"] = df["mapping_note"].fillna("")
    atomic_write_csv(df, csv_file, index=False, encoding="utf-8-sig")

    print(f"File saved to {csv_file} with {df.shape[0]} rows.")


def _is_generated_csv(path: str | Path) -> bool:
    """Skip derived pipeline outputs when mapping original input symbols."""
    name = Path(path).name
    return (
        name.endswith("_human.csv")
        or name.endswith("_mouse.csv")
        or name.endswith("_classification.csv")
    )


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Map input fly gene symbols to FlyBase FBgn IDs in CSV files."
    )
    parser.add_argument("input_directory", help="Directory containing input CSV files")
    parser.add_argument(
        "input_gene_col",
        nargs="?",
        default="ext_gene",
        help="Input gene symbol column (default: ext_gene)",
    )
    parser.add_argument(
        "--flybase-data-dir",
        default=os.environ.get("FLYBASE_DATA_DIR", str(DEFAULT_FLYBASE_DATA)),
        help="Path to FlyBase data directory (default: Data/FlyBase path)",
    )
    args = parser.parse_args(argv)

    input_dir = Path(args.input_directory)
    if not input_dir.is_dir():
        print(f"Error: '{input_dir}' is not a valid directory")
        return 1

    flybase_data_dir = Path(args.flybase_data_dir)
    missing_flybase = ensure_flybase_data_files(flybase_data_dir, download_missing=True)
    synonym_missing = [item for item in missing_flybase if item["pattern"] == "fb_synonym"]
    if synonym_missing:
        print(f"Error: FlyBase synonym data is missing under: {flybase_data_dir}")
        for requirement in synonym_missing:
            print(f"  - Could not retrieve {requirement['label']}: {requirement['url']}")
        return 1

    csv_files = [
        csv_file
        for csv_file in sorted(glob(str(input_dir / "*.csv")))
        if not _is_generated_csv(csv_file)
    ]
    if not csv_files:
        print(f"No CSV files found in {input_dir}")
        return 0

    evidence = load_evidence()
    if evidence is None:
        print(
            f"Note: FlyBase citation index not found ({EVIDENCE_PATH}); shared synonyms will stay unmapped. "
            "Build it with: python HelperScripts/build_synonym_evidence.py"
        )
    try:
        symbol_index = load_symbol_index(flybase_data_dir, evidence)
    except Exception as e:
        print(f"Error loading FlyBase synonym mappings: {e}")
        return 1

    for csv_file in csv_files:
        process_csv_file(csv_file, args.input_gene_col, symbol_index)

    print("All files processed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
