"""Shared data-loading utilities for species gene catalogs and references."""

from __future__ import annotations

import gzip
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import pandas as pd
import requests


NCBI_GENE2PUBMED_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz"
NCBI_GENERIF_URL = "https://ftp.ncbi.nlm.nih.gov/gene/GeneRIF/generifs_basic.gz"


def default_data_root() -> Path:
    root = os.environ.get("FLAI_CACHE_DIR")
    if root:
        return Path(root).expanduser()
    return Path.home() / ".cache" / "flai-gene-classification"


def download_if_needed(url: str, dest: Path, refresh: bool = False, timeout: int = 120) -> Path:
    dest = Path(dest)
    if dest.exists() and not refresh:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    with requests.get(url, stream=True, timeout=timeout) as response:
        response.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
    tmp.replace(dest)
    return dest


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [
        re.sub(r"\s+", "_", str(col).strip().lstrip("#").lower())
        for col in df.columns
    ]
    return df


def split_values(value: object) -> list[str]:
    text = str(value or "").strip()
    if not text or text in {"-", "nan", "None"}:
        return []
    parts = re.split(r"[|;,]", text)
    return [p.strip() for p in parts if p.strip()]


def _read_taxon_chunks(path: Path, chunksize: int = 500_000):
    for chunk in pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, chunksize=chunksize):
        yield normalize_columns(chunk)


def load_gene2pubmed_pmids(path: Path, taxon_id: int, entrez_ids: Iterable[str]) -> dict[str, set[str]]:
    target_ids = {str(x).strip() for x in entrez_ids if str(x).strip()}
    out: dict[str, set[str]] = defaultdict(set)
    if not target_ids:
        return out
    for chunk in _read_taxon_chunks(path):
        tax_col = "tax_id"
        gene_col = "geneid" if "geneid" in chunk.columns else "gene_id"
        pmid_col = "pubmed_id"
        if not {tax_col, gene_col, pmid_col}.issubset(set(chunk.columns)):
            continue
        subset = chunk[(chunk[tax_col] == str(taxon_id)) & (chunk[gene_col].isin(target_ids))]
        for _, row in subset.iterrows():
            gene_id = str(row[gene_col]).strip()
            pmid = str(row[pmid_col]).strip()
            if gene_id and pmid.isdigit():
                out[gene_id].add(pmid)
    return out


def load_generif_pmids_and_snippets(path: Path, taxon_id: int, entrez_ids: Iterable[str]) -> tuple[dict[str, set[str]], dict[str, list[tuple[str, str]]]]:
    target_ids = {str(x).strip() for x in entrez_ids if str(x).strip()}
    pmids_by_gene: dict[str, set[str]] = defaultdict(set)
    snippets_by_gene: dict[str, list[tuple[str, str]]] = defaultdict(list)
    if not target_ids:
        return pmids_by_gene, snippets_by_gene

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        columns = [re.sub(r"\s+", "_", h.strip().lstrip("#").lower()) for h in header]
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(columns):
                continue
            row = dict(zip(columns, parts))
            tax_id = row.get("tax_id") or row.get("taxid") or row.get("tax")
            gene_id = row.get("gene_id") or row.get("geneid")
            if tax_id != str(taxon_id) or gene_id not in target_ids:
                continue
            pmid_text = row.get("pubmed_id_(pmid)_list") or row.get("pubmed_id") or row.get("pmid") or ""
            rif_text = row.get("generif_text") or row.get("gene_rif_text") or row.get("text") or ""
            for pmid in split_values(pmid_text):
                if pmid.isdigit():
                    pmids_by_gene[gene_id].add(pmid)
                    if rif_text:
                        snippets_by_gene[gene_id].append((pmid, rif_text.strip()))
    return pmids_by_gene, snippets_by_gene


def download_uniprot_tsv(path: Path, organism_id: int, refresh: bool = False) -> Path:
    if path.exists() and not refresh:
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    query = f"organism_id:{int(organism_id)} AND reviewed:true"
    params = {
        "query": query,
        "fields": "accession,gene_names,xref_geneid,lit_pubmed_id",
        "format": "tsv",
        "size": "500",
    }
    url = "https://rest.uniprot.org/uniprotkb/search"
    lines: list[str] = []
    while url:
        response = requests.get(url, params=params if "?" not in url else None, timeout=120)
        response.raise_for_status()
        text_lines = response.text.splitlines()
        if not lines:
            lines.extend(text_lines)
        else:
            lines.extend(text_lines[1:])
        link = response.headers.get("Link", "")
        next_url = ""
        for part in link.split(","):
            if 'rel="next"' in part:
                match = re.search(r"<([^>]+)>", part)
                if match:
                    next_url = match.group(1)
                    break
        url = next_url
        params = None
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def load_uniprot_pmids(path: Path, entrez_ids: Iterable[str]) -> dict[str, set[str]]:
    target_ids = {str(x).strip() for x in entrez_ids if str(x).strip()}
    out: dict[str, set[str]] = defaultdict(set)
    if not path.exists() or not target_ids:
        return out
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    df = normalize_columns(df)
    gene_col = next((c for c in df.columns if c in {"geneid", "xref_geneid", "gene_id"}), "")
    pmid_col = next((c for c in df.columns if "pubmed" in c), "")
    if not gene_col or not pmid_col:
        return out
    for _, row in df.iterrows():
        gene_ids = [g for value in split_values(row.get(gene_col, "")) for g in split_values(value)]
        pmids = [p for p in split_values(row.get(pmid_col, "")) if p.isdigit()]
        for gene_id in gene_ids:
            if gene_id in target_ids:
                out[gene_id].update(pmids)
    return out
