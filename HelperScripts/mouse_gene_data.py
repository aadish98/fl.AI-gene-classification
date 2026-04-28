"""Mouse gene catalog and reference-source loaders."""

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


TAXON_ID = 10090
MGI_MARKERS_URL = "https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt"
MGI_ENTREZ_URL = "https://www.informatics.jax.org/downloads/reports/MGI_EntrezGene.rpt"

_DATA_DIR = default_data_root() / "mouse-data"
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
    download_if_needed(MGI_MARKERS_URL, _path("MRK_List2.rpt"), refresh=_REFRESH)
    download_if_needed(MGI_ENTREZ_URL, _path("MGI_EntrezGene.rpt"), refresh=_REFRESH)
    download_if_needed(NCBI_GENE2PUBMED_URL, _path("gene2pubmed.gz"), refresh=_REFRESH)
    download_if_needed(NCBI_GENERIF_URL, _path("generifs_basic.gz"), refresh=_REFRESH)
    download_uniprot_tsv(_path("uniprot_sprot_mouse.tsv"), TAXON_ID, refresh=_REFRESH)


def _find_col(columns, *candidates: str) -> str:
    normalized = {str(c).lower().replace(" ", "_"): c for c in columns}
    for candidate in candidates:
        key = candidate.lower().replace(" ", "_")
        if key in normalized:
            return normalized[key]
    for col in columns:
        low = str(col).lower()
        if all(part in low for part in candidates[0].lower().split()):
            return col
    return ""


@lru_cache(maxsize=1)
def load_gene_catalog() -> GeneCatalog:
    ensure_files()
    markers = pd.read_csv(_path("MRK_List2.rpt"), sep="\t", dtype=str, keep_default_na=False)
    entrez = pd.read_csv(_path("MGI_EntrezGene.rpt"), sep="\t", dtype=str, keep_default_na=False)
    markers = normalize_columns(markers)
    entrez = normalize_columns(entrez)

    mgi_col = _find_col(markers.columns, "mgi_marker_accession_id", "mgi_accession_id", "mgi_id")
    symbol_col = _find_col(markers.columns, "marker_symbol", "symbol")
    name_col = _find_col(markers.columns, "marker_name", "name")
    synonyms_col = _find_col(markers.columns, "marker_synonyms", "synonyms")
    entrez_mgi_col = _find_col(entrez.columns, "mgi_marker_accession_id", "mgi_accession_id", "mgi_id")
    entrez_col = _find_col(entrez.columns, "entrez_gene_id", "geneid", "gene_id")

    mgi_to_entrez = {}
    if entrez_mgi_col and entrez_col:
        for _, row in entrez.iterrows():
            mgi_id = str(row.get(entrez_mgi_col, "") or "").strip()
            entrez_id = str(row.get(entrez_col, "") or "").strip()
            if mgi_id and entrez_id:
                mgi_to_entrez[mgi_id] = entrez_id

    catalog = GeneCatalog(species="mouse")
    if not mgi_col or not symbol_col:
        return catalog

    for _, row in markers.iterrows():
        mgi_id = str(row.get(mgi_col, "") or "").strip()
        symbol = str(row.get(symbol_col, "") or "").strip()
        entrez_id = mgi_to_entrez.get(mgi_id, "")
        if not mgi_id or not symbol or not entrez_id:
            continue
        synonyms = {symbol}
        if name_col:
            synonyms.update(split_values(row.get(name_col, "")))
        if synonyms_col:
            synonyms.update(split_values(row.get(synonyms_col, "")))
        record = GeneRecord(
            species="mouse",
            gene_id=entrez_id,
            symbol=symbol,
            authority_id=mgi_id,
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
    return load_uniprot_pmids(_path("uniprot_sprot_mouse.tsv"), entrez_ids)
