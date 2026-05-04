"""DIOPT ortholog lookup client with on-disk caching."""

from __future__ import annotations

import json
import os
import re
import threading
import time
import uuid
from pathlib import Path
from typing import Any

import requests


DIOPT_INPUT_TAXON_FLY = 7227
DIOPT_TAXON_HUMAN = 9606
DIOPT_TAXON_MOUSE = 10090
DEFAULT_ERROR_CACHE_TTL_SECONDS = 6 * 60 * 60

FILTER_ALIASES = {
    "exclude_low_score_2": "exclude_score_less_2",
    "exclude_low_score_1": "exclude_score_less_1",
    "exclude_score_less_2": "exclude_score_less_2",
    "exclude_score_less_1": "exclude_score_less_1",
    "exclude_low_ranked": "exclude_low_ranked",
    "best_match": "best_match",
    "none": "none",
}


def default_cache_dir() -> Path:
    root = os.environ.get("FLAI_CACHE_DIR")
    if root:
        return Path(root).expanduser() / "diopt"
    return Path.home() / ".cache" / "flai-gene-classification" / "diopt"


def _error_cache_ttl_seconds() -> int:
    raw = os.environ.get("FLAI_DIOPT_ERROR_CACHE_TTL_SECONDS")
    if not raw:
        return DEFAULT_ERROR_CACHE_TTL_SECONDS
    try:
        return int(raw)
    except ValueError:
        return DEFAULT_ERROR_CACHE_TTL_SECONDS


def _clean_text(value: Any) -> str:
    if isinstance(value, (list, tuple, set)):
        return ";".join(str(v).strip() for v in value if str(v).strip())
    return re.sub(r"\s+", " ", str(value or "")).strip()


def _first_value(row: dict[str, Any], *keys: str) -> str:
    lowered = {str(k).lower(): v for k, v in row.items()}
    for key in keys:
        if key in row and row[key] not in (None, ""):
            return _clean_text(row[key])
        low = key.lower()
        if low in lowered and lowered[low] not in (None, ""):
            return _clean_text(lowered[low])
    return ""


def _score_from_row(row: dict[str, Any]) -> int:
    raw = _first_value(
        row,
        "weighted_score",
        "weighted score",
        "score",
        "diopt_score",
        "DIOPT score",
    )
    match = re.search(r"\d+", raw)
    return int(match.group(0)) if match else 0


def _normalize_row(row: dict[str, Any], output_taxon: int) -> dict[str, Any] | None:
    symbol = _first_value(
        row,
        "symbol",
        "gene_symbol",
        "ortholog_symbol",
        "out_symbol",
        "to_symbol",
        "target_symbol",
        "human_symbol",
        "mouse_symbol",
    )
    entrez = _first_value(
        row,
        "entrez_gene_id",
        "entrez_id",
        "entrez",
        "geneid",
        "gene_id",
        "out_gene_id",
        "to_gene_id",
        "target_gene_id",
        "NCBI Gene ID",
    )
    hgnc_id = _first_value(row, "hgnc_id", "HGNC ID", "hgnc")
    mgi_id = _first_value(row, "mgi_id", "MGI ID", "mgi")
    species_specific_id = _first_value(row, "species_specific_geneid", "species_specific_gene_id")
    species_specific_type = _first_value(row, "species_specific_geneid_type", "species_specific_gene_id_type").upper()
    if species_specific_id and species_specific_type == "HGNC":
        hgnc_id = hgnc_id or f"HGNC:{species_specific_id}"
    if species_specific_id and species_specific_type == "MGI":
        mgi_id = mgi_id or f"MGI:{species_specific_id}"
    authority_id = hgnc_id if output_taxon == DIOPT_TAXON_HUMAN else mgi_id
    best_raw = _first_value(row, "best_score", "best_score_yn", "best match", "best_match")
    algorithms = _first_value(
        row,
        "supporting_algorithms",
        "ortholog_sources",
        "sources",
        "algorithms",
        "methods",
    )
    score = _score_from_row(row)

    if not symbol and not entrez and not authority_id:
        return None

    return {
        "gene_symbol": symbol,
        "entrez_gene_id": entrez,
        "authority_id": authority_id,
        "hgnc_id": hgnc_id,
        "mgi_id": mgi_id,
        "diopt_score": score,
        "diopt_best_match": best_raw,
        "diopt_supporting_algorithms": algorithms,
        "raw": row,
    }


def _extract_rows(payload: Any) -> list[dict[str, Any]]:
    if isinstance(payload, list):
        return [x for x in payload if isinstance(x, dict)]
    if not isinstance(payload, dict):
        return []
    for key in ("data", "results", "orthologs", "rows", "result"):
        value = payload.get(key)
        if isinstance(value, list):
            return [x for x in value if isinstance(x, dict)]
        if isinstance(value, dict):
            nested = _extract_rows(value)
            if nested:
                return nested
    if {"geneid", "symbol", "score"} & set(str(k).lower() for k in payload.keys()):
        return [payload]
    rows = []
    for value in payload.values():
        if isinstance(value, dict):
            rows.extend(_extract_rows(value))
        elif isinstance(value, list):
            rows.extend([item for item in value if isinstance(item, dict)])
    return rows


class DIOPTClient:
    """Small resilient DIOPT client.

    DIOPT has changed versioned API paths over time, so the client tries a small
    ordered set of candidate endpoints and caches the first successful payload.
    """

    def __init__(
        self,
        cache_dir: str | Path | None = None,
        timeout: int = 30,
        error_cache_ttl_seconds: int | None = None,
    ):
        self.cache_dir = Path(cache_dir) if cache_dir else default_cache_dir()
        self.timeout = timeout
        self.error_cache_ttl_seconds = (
            error_cache_ttl_seconds
            if error_cache_ttl_seconds is not None
            else _error_cache_ttl_seconds()
        )
        self._thread_local = threading.local()
        self._preferred_endpoint_lock = threading.Lock()
        self._preferred_endpoint_by_key: dict[tuple[int, str], str] = {}

    def _session(self) -> requests.Session:
        session = getattr(self._thread_local, "session", None)
        if session is None:
            session = requests.Session()
            self._thread_local.session = session
        return session

    def _cache_path(self, fbgn_id: str, output_taxon: int, diopt_filter: str) -> Path:
        safe_fbgn = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(fbgn_id).strip())
        safe_filter = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(diopt_filter).strip())
        return self.cache_dir / f"{safe_fbgn}_{int(output_taxon)}_{safe_filter}.json"

    def _candidate_urls(
        self,
        fbgn_id: str,
        output_taxon: int,
        diopt_filter: str,
        query_id: str | None = None,
    ) -> list[tuple[str, str]]:
        base = "https://www.flyrnai.org/tools/diopt/web/diopt_api"
        entrez_query = str(query_id or fbgn_id).strip()
        fbgn_query = str(fbgn_id).strip()
        filt = FILTER_ALIASES.get(str(diopt_filter).strip(), str(diopt_filter).strip())
        candidates = [
            (
                "10_entrez",
                f"{base}/10/get_orthologs_from_entrez/"
                f"{DIOPT_INPUT_TAXON_FLY}/{entrez_query}/{int(output_taxon)}/{filt}",
            ),
            (
                "9_entrez",
                f"{base}/9/get_orthologs_from_entrez/"
                f"{DIOPT_INPUT_TAXON_FLY}/{entrez_query}/{int(output_taxon)}/{filt}",
            ),
            (
                "10_flybase",
                f"{base}/10/get_orthologs_from_flybase/"
                f"{DIOPT_INPUT_TAXON_FLY}/{fbgn_query}/{int(output_taxon)}/{filt}",
            ),
            (
                "9_flybase",
                f"{base}/9/get_orthologs_from_flybase/"
                f"{DIOPT_INPUT_TAXON_FLY}/{fbgn_query}/{int(output_taxon)}/{filt}",
            ),
            (
                "v10_flybase",
                f"{base}/v10/get_orthologs_from_flybase/"
                f"{DIOPT_INPUT_TAXON_FLY}/{fbgn_query}/{int(output_taxon)}/{filt}",
            ),
            (
                "v9_flybase",
                f"{base}/v9/get_orthologs_from_flybase/"
                f"{DIOPT_INPUT_TAXON_FLY}/{fbgn_query}/{int(output_taxon)}/{filt}",
            ),
        ]
        preference_key = (int(output_taxon), filt)
        with self._preferred_endpoint_lock:
            preferred = self._preferred_endpoint_by_key.get(preference_key)
        if preferred:
            candidates.sort(key=lambda item: 0 if item[0] == preferred else 1)
        return candidates

    def _load_cache(self, path: Path) -> Any | None:
        if not path.exists():
            return None
        try:
            with open(path, "r", encoding="utf-8") as f:
                cached = json.load(f)
        except Exception:
            return None
        if isinstance(cached, dict) and cached.get("error"):
            try:
                cached_at = float(cached.get("cached_at", 0))
            except (TypeError, ValueError):
                cached_at = 0
            if self.error_cache_ttl_seconds <= 0:
                return None
            if not cached_at or time.time() - cached_at > self.error_cache_ttl_seconds:
                return None
        return cached

    def _save_cache(self, path: Path, payload: Any):
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(f".{path.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp")
        try:
            with open(tmp, "w", encoding="utf-8") as f:
                json.dump(payload, f, ensure_ascii=True, indent=2, sort_keys=True)
            tmp.replace(path)
        finally:
            tmp.unlink(missing_ok=True)

    def fetch_raw(
        self,
        fbgn_id: str,
        output_taxon: int,
        diopt_filter: str = "exclude_low_score_2",
        query_id: str | None = None,
    ) -> Any:
        path = self._cache_path(fbgn_id, output_taxon, diopt_filter)
        cached = self._load_cache(path)
        if cached is not None:
            return cached

        errors: list[str] = []
        for url in self._candidate_urls(fbgn_id, output_taxon, diopt_filter, query_id=query_id):
            endpoint_id, endpoint_url = url
            try:
                response = self._session().get(endpoint_url, timeout=self.timeout)
                if response.status_code != 200:
                    errors.append(f"{response.status_code} {endpoint_url}")
                    continue
                try:
                    payload = response.json()
                except Exception:
                    errors.append(f"non-json {endpoint_url}")
                    continue
                filt = FILTER_ALIASES.get(str(diopt_filter).strip(), str(diopt_filter).strip())
                with self._preferred_endpoint_lock:
                    self._preferred_endpoint_by_key[(int(output_taxon), filt)] = endpoint_id
                self._save_cache(path, payload)
                return payload
            except Exception as exc:
                errors.append(f"{type(exc).__name__}: {endpoint_url}")
            time.sleep(0.2)

        payload = {
            "error": "DIOPT lookup failed",
            "fbgn_id": fbgn_id,
            "errors": errors,
            "cached_at": time.time(),
        }
        self._save_cache(path, payload)
        return payload

    def orthologs(
        self,
        fbgn_id: str,
        *,
        output_taxon: int,
        diopt_filter: str = "exclude_low_score_2",
        query_id: str | None = None,
    ) -> list[dict[str, Any]]:
        payload = self.fetch_raw(fbgn_id, output_taxon, diopt_filter, query_id=query_id)
        rows = _extract_rows(payload)
        normalized = []
        for row in rows:
            norm = _normalize_row(row, int(output_taxon))
            if norm:
                normalized.append(norm)
        return normalized


def fbgn_to_orthologs(
    fbgn_id: str,
    *,
    output_taxon: int,
    diopt_filter: str = "exclude_low_score_2",
    query_id: str | None = None,
    cache_dir: str | Path | None = None,
) -> list[dict[str, Any]]:
    return DIOPTClient(cache_dir=cache_dir).orthologs(
        fbgn_id,
        output_taxon=output_taxon,
        diopt_filter=diopt_filter,
        query_id=query_id,
    )
