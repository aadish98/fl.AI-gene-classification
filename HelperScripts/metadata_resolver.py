#!/usr/bin/env python3
"""Shared PMCID/PMID metadata resolution helpers for pipeline and maintenance scripts."""

from __future__ import annotations

from datetime import datetime, timezone
import re
from typing import Any, Callable, Dict, Optional

import requests


CacheGetter = Callable[[str], Optional[dict]]
CacheSetter = Callable[..., None]


def normalize_pmid(value: Any) -> str:
    text = str(value or "").strip()
    return text if text.isdigit() else ""


def normalize_pmcid(value: Any) -> str:
    text = str(value or "").strip().upper()
    if not text:
        return ""
    m = re.fullmatch(r"(PMC)?(\d+)", text)
    if not m:
        return ""
    return f"PMC{m.group(2)}"


def normalize_doi(value: Any) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    text = re.sub(r"^https?://(dx\.)?doi\.org/", "", text, flags=re.IGNORECASE)
    return text.strip().lower()


def normalize_authors(raw_authors: Any) -> list[str]:
    if raw_authors is None:
        return []

    if isinstance(raw_authors, str):
        text = raw_authors.strip()
        if not text:
            return []
        parts = [p.strip() for p in text.split(";")]
        return [p for p in parts if p]

    cleaned: list[str] = []
    try:
        for item in list(raw_authors):
            value = str(item or "").strip()
            if value:
                cleaned.append(value)
    except Exception:
        return []
    return cleaned


def authors_display(authors: Any) -> str:
    values = normalize_authors(authors)
    return "; ".join(values)


def metadata_is_reasonable(meta: dict) -> bool:
    journal = str(meta.get("journal", "") or "").strip()
    authors_str = str(meta.get("authors_display", "") or "").strip()
    if journal.startswith("[") and journal.endswith("]"):
        return False
    if authors_str.startswith("[") and authors_str.endswith("]"):
        return False
    return True


def _article_to_metadata(article: Any) -> dict:
    year_value = str(getattr(article, "year", "") or "").strip()
    year_digits = "".join(ch for ch in year_value if ch.isdigit())
    year = year_digits[:4] if len(year_digits) >= 4 else ""
    authors = normalize_authors(getattr(article, "authors", []))
    return {
        "pmid": normalize_pmid(getattr(article, "pmid", "")),
        "pmcid": normalize_pmcid(getattr(article, "pmc", "")),
        "title": str(getattr(article, "title", "") or "").strip(),
        "abstract": str(getattr(article, "abstract", "") or "").strip(),
        "year": year,
        "journal": str(getattr(article, "journal", "") or "").strip(),
        "authors": authors,
        "authors_display": authors_display(authors),
        "doi": normalize_doi(getattr(article, "doi", "")),
        "source": "metapub",
        "updated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }


def pmcid_to_pmid(pmcid: str, api_key: str = "", email: str = "aadish98@gmail.com") -> str:
    pmcid_norm = normalize_pmcid(pmcid)
    if not pmcid_norm:
        return ""

    params = {
        "tool": "fly_gene_classifier",
        "email": email,
        "ids": pmcid_norm,
        "format": "json",
    }
    if api_key:
        params["api_key"] = api_key

    try:
        r = requests.get(
            "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/",
            params=params,
            timeout=20,
        )
        if r.status_code != 200:
            return ""
        data = r.json() or {}
        recs = data.get("records", []) or []
        if not recs:
            return ""
        return normalize_pmid(recs[0].get("pmid", ""))
    except Exception:
        return ""


def resolve_reference_metadata(
    paper_id: str,
    fetcher: Any,
    *,
    api_key: str = "",
    pmid_hint: str = "",
    cache_getter: Optional[CacheGetter] = None,
    cache_setter: Optional[CacheSetter] = None,
    email: str = "aadish98@gmail.com",
) -> dict:
    """Resolve normalized metadata from Paper_ID (PMCID/PMID)."""
    pid = str(paper_id or "").strip()
    pmid = normalize_pmid(pid) or normalize_pmid(pmid_hint)
    pmcid = normalize_pmcid(pid)

    if not pmid and pmcid:
        pmid = pmcid_to_pmid(pmcid, api_key=api_key, email=email)

    cached: dict = {}
    if pmid and cache_getter:
        cached = cache_getter(pmid) or {}

    result = {
        "pmid": pmid,
        "pmcid": pmcid,
        "title": str(cached.get("title", "") or "").strip(),
        "abstract": str(cached.get("abstract", "") or "").strip(),
        "year": str(cached.get("year", "") or "").strip(),
        "journal": str(cached.get("journal", "") or "").strip(),
        "authors": normalize_authors(cached.get("authors", [])),
        "authors_display": authors_display(cached.get("authors", [])),
        "doi": normalize_doi(cached.get("doi", "")),
        "source": str(cached.get("source", "cache") or "cache"),
        "updated_at": str(cached.get("updated_at", "") or "").strip(),
    }

    needs_fetch = (
        not result["title"]
        or not result["abstract"]
        or not result["journal"]
        or not result["authors"]
        or not result["year"]
    )

    if needs_fetch:
        article = None
        try:
            if pmid:
                article = fetcher.article_by_pmid(pmid)
            elif pmcid:
                article = fetcher.article_by_pmcid(pmcid)
        except Exception:
            article = None

        if article is not None:
            art = _article_to_metadata(article)
            if art.get("pmid"):
                pmid = art["pmid"]
                result["pmid"] = pmid
            if art.get("pmcid"):
                result["pmcid"] = art["pmcid"]
            for field in ["title", "abstract", "year", "journal", "doi", "source", "updated_at"]:
                if art.get(field):
                    result[field] = art[field]
            if art.get("authors"):
                result["authors"] = art["authors"]
                result["authors_display"] = art["authors_display"]

            if pmid and cache_setter:
                cache_setter(
                    pmid,
                    result.get("title", ""),
                    result.get("abstract", ""),
                    year=result.get("year", ""),
                    journal=result.get("journal", ""),
                    authors=result.get("authors", []),
                    doi=result.get("doi", ""),
                    pmcid=result.get("pmcid", ""),
                    source=result.get("source", "metapub"),
                    updated_at=result.get("updated_at", ""),
                )

    result["authors"] = normalize_authors(result.get("authors", []))
    result["authors_display"] = authors_display(result["authors"])
    result["doi"] = normalize_doi(result.get("doi", ""))
    if not metadata_is_reasonable(result):
        result["journal"] = str(result.get("journal", "") or "").strip("[]")
        result["authors"] = normalize_authors(result.get("authors", []))
        result["authors_display"] = authors_display(result["authors"])
    if not result.get("updated_at"):
        result["updated_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    return result
