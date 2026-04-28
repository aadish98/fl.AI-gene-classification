"""Human gene catalog and reference-source loaders."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pandas as pd

from HelperScripts.gene_models import GeneCatalog, GeneRecord
from HelperScripts.species_data_utils import (
    NCBI_GENERIF_URL,
    NCBI_GENE2PUBMED_URL,
    default_data_root,
    download_if_needed,
    download_uniprot_tsv,
    load_gene2pubmed_pmids,
    load_generif_pmids_and_snippets,
    load_uniprot_pmids,
    normalize_columns,
    split_values,
)


TAXON_ID = 9606
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"

_DATA_DIR = default_data_root() / "human-data"
_REFRESH = False
_GENERIF_SNIPPETS: dict[str, list[tuple[str, str]]] = {}


def configure(data_dir: str | Path | None = None, refresh: bool = False):
    global _DATA_DIR, _REFRESH
    if data_dir:
        _DATA_DIR = Path(data_dir).expanduser()
    _REFRESH = bool(refresh)
    load_gene_catalog.cache_clear()


def _path(name: str) -> Path:
    return _DATA_DIR / name


def ensure_files():
    download_if_needed(HGNC_URL, _path("hgnc_complete_set.txt"), refresh=_REFRESH)
    download_if_needed(NCBI_GENE2PUBMED_URL, _path("gene2pubmed.gz"), refresh=_REFRESH)
    download_if_needed(NCBI_GENERIF_URL, _path("generifs_basic.gz"), refresh=_REFRESH)
    download_uniprot_tsv(_path("uniprot_sprot_human.tsv"), TAXON_ID, refresh=_REFRESH)


@lru_cache(maxsize=1)
def load_gene_catalog() -> GeneCatalog:
    ensure_files()
    df = pd.read_csv(_path("hgnc_complete_set.txt"), sep="\t", dtype=str, keep_default_na=False)
    df = normalize_columns(df)
    catalog = GeneCatalog(species="human")
    for _, row in df.iterrows():
        entrez_id = str(row.get("entrez_id", "") or "").strip()
        symbol = str(row.get("symbol", "") or "").strip()
        hgnc_id = str(row.get("hgnc_id", "") or "").strip()
        if not entrez_id or not symbol:
            continue
        synonyms = {symbol}
        for field in ("name", "alias_symbol", "prev_symbol", "alias_name", "prev_name"):
            synonyms.update(split_values(row.get(field, "")))
        record = GeneRecord(
            species="human",
            gene_id=entrez_id,
            symbol=symbol,
            authority_id=hgnc_id,
            synonyms={s for s in synonyms if s},
        )
        catalog.genes[entrez_id] = record
        for name in record.synonyms:
            catalog.symbol_to_gene_id.setdefault(name, entrez_id)
    return catalog


def get_gene2pubmed_pmids(entrez_ids) -> dict[str, set[str]]:
    ensure_files()
    return load_gene2pubmed_pmids(_path("gene2pubmed.gz"), TAXON_ID, entrez_ids)


def get_generif_pmids(entrez_ids) -> dict[str, set[str]]:
    ensure_files()
    global _GENERIF_SNIPPETS
    pmids, snippets = load_generif_pmids_and_snippets(_path("generifs_basic.gz"), TAXON_ID, entrez_ids)
    _GENERIF_SNIPPETS.update(snippets)
    return pmids


def get_generif_snippets(entrez_id: str) -> list[tuple[str, str]]:
    return list(_GENERIF_SNIPPETS.get(str(entrez_id).strip(), []))


def get_uniprot_pmids(entrez_ids) -> dict[str, set[str]]:
    ensure_files()
    return load_uniprot_pmids(_path("uniprot_sprot_human.tsv"), entrez_ids)
