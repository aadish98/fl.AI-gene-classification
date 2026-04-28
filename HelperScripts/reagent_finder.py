"""Fly reagent extraction helpers, isolated from the classification pipeline."""

from __future__ import annotations

import re
from collections.abc import Callable
from typing import Any, Optional

from pydantic import BaseModel, Field


StructuredParser = Callable[..., Any]


class ReferenceReagent(BaseModel):
    stock_id: str = ""
    collection: str = ""
    reagent_type: str = ""
    evidence_snippet: str = ""
    functional_validity: str = ""
    reagent_name: str = ""


class ReagentExtraction(BaseModel):
    reagents: list[ReferenceReagent] = Field(default_factory=list)


def _clean_text(value: Any) -> str:
    """Collapse repeated whitespace and coerce to string."""
    return re.sub(r"\s+", " ", str(value or "")).strip()


def _model_dump(value: Any) -> dict[str, Any]:
    """Return a plain dict for a parsed Pydantic model."""
    if value is None:
        return {}
    if hasattr(value, "model_dump"):
        return value.model_dump()
    if hasattr(value, "dict"):
        return value.dict()
    if isinstance(value, dict):
        return dict(value)
    return {}


def _normalize_collection_name(collection: str) -> str:
    """Normalize collection names into stable export values."""
    cleaned = _clean_text(collection)
    lowered = cleaned.lower()
    if not cleaned:
        return ""
    if "bloomington" in lowered or "bdsc" in lowered:
        return "BDSC"
    if "vienna" in lowered or "vdrc" in lowered:
        return "VDRC"
    if "national institute of genetics" in lowered or re.search(r"\bnig\b", lowered):
        return "NIG"
    if "kyoto" in lowered or "dgrc" in lowered or "drosophila genetic resource center" in lowered:
        return "Kyoto"
    return cleaned


def _normalize_stock_id(stock_id: str, collection: str) -> str:
    """Normalize stock IDs to stable, dedupe-friendly forms."""
    cleaned = _clean_text(stock_id)
    cleaned = cleaned.replace("#", "")
    cleaned = re.sub(r"\s+", "", cleaned)
    if not cleaned:
        return ""

    if collection == "BDSC":
        match = re.fullmatch(r"(?:BL|BDSC)?(\d+)", cleaned, flags=re.IGNORECASE)
        if match:
            return f"BL{match.group(1)}"
    if collection == "VDRC":
        match = re.fullmatch(r"(?:VDRC)?v?(\d+)", cleaned, flags=re.IGNORECASE)
        if match:
            return f"v{match.group(1)}"
    return cleaned


def _normalize_reagent_record(raw_reagent: Any) -> Optional[dict[str, str]]:
    """Normalize one reagent record and enforce canonical (stock_id, collection) pairs."""
    parsed = _model_dump(raw_reagent)
    collection = _normalize_collection_name(parsed.get("collection", ""))
    stock_id = _normalize_stock_id(parsed.get("stock_id", ""), collection)
    if not stock_id or not collection:
        return None
    return {
        "stock_id": stock_id,
        "collection": collection,
        "reagent_type": _clean_text(parsed.get("reagent_type", "")),
        "evidence_snippet": _clean_text(parsed.get("evidence_snippet", "")),
        "functional_validity": _clean_text(parsed.get("functional_validity", "")),
        "reagent_name": _clean_text(parsed.get("reagent_name", "")),
    }


def extract_reference_reagents(
    text_chunk: str,
    gene_symbol: str,
    fbgn_id: str,
    synonyms: list[str],
    *,
    parse_structured_completion: StructuredParser,
    model_name: Optional[str] = None,
    title: str = "",
    abstract: str = "",
    chunk_index: int = 1,
    total_chunks: int = 1,
) -> list[dict[str, str]]:
    """Extract canonical reagent pairs plus metadata from one paper chunk."""
    synonyms_str = ", ".join(sorted(set(synonyms))) if synonyms else ""
    sys = (
        "You are an expert biomedical assistant. Use ONLY the provided title, abstract, "
        "and text chunk. Extract reagent records only when they pertain to the target gene."
    )
    user = f"""Extract reagents for gene {gene_symbol} (FBgn: {fbgn_id}; Synonyms: {synonyms_str}).

Return JSON with a single field:
- reagents: list of objects with stock_id, collection, reagent_type, evidence_snippet, functional_validity, reagent_name

Rules:
- Only include a reagent if both stock_id and collection can be identified
- Focus only on reagents that pertain to {gene_symbol}
- Deduplicate exact repeats within this chunk
- If no qualifying reagents are present, return an empty list
- Use concise evidence_snippet text quoted or paraphrased from the chunk
- Normalize collection names when clear, including BDSC, VDRC, NIG, and Kyoto

Chunk {chunk_index} of {total_chunks}
Title: {title or ''}
Abstract: {abstract or ''}

Text chunk:
{text_chunk or ''}"""
    try:
        out = parse_structured_completion(
            [{"role": "system", "content": sys}, {"role": "user", "content": user}],
            response_format=ReagentExtraction,
            model_name=model_name,
            max_output_tokens=1800,
        )
        parsed = _model_dump(out)
        normalized = []
        for reagent in parsed.get("reagents", []) or []:
            record = _normalize_reagent_record(reagent)
            if record:
                normalized.append(record)
        return normalized
    except Exception:
        return []


def _merge_reagent_records(base: dict[str, str], candidate: dict[str, str]) -> dict[str, str]:
    """Merge duplicate reagent records while keeping the richest metadata."""
    merged = dict(base or {})
    for field in ("reagent_name", "reagent_type", "functional_validity"):
        if not merged.get(field) and candidate.get(field):
            merged[field] = candidate[field]
    if len(candidate.get("evidence_snippet", "")) > len(merged.get("evidence_snippet", "")):
        merged["evidence_snippet"] = candidate.get("evidence_snippet", "")
    return merged


def _merge_deduplicated_reagents(reagent_records: list[dict[str, str]]) -> list[dict[str, str]]:
    """Merge reagent records across chunks using the canonical pair as the key."""
    deduped: dict[tuple[str, str], dict[str, str]] = {}
    for record in reagent_records or []:
        stock_id = _clean_text(record.get("stock_id", ""))
        collection = _normalize_collection_name(record.get("collection", ""))
        if not stock_id or not collection:
            continue
        key = (stock_id.lower(), collection.lower())
        normalized_record = {
            "stock_id": stock_id,
            "collection": collection,
            "reagent_type": _clean_text(record.get("reagent_type", "")),
            "evidence_snippet": _clean_text(record.get("evidence_snippet", "")),
            "functional_validity": _clean_text(record.get("functional_validity", "")),
            "reagent_name": _clean_text(record.get("reagent_name", "")),
        }
        if key in deduped:
            deduped[key] = _merge_reagent_records(deduped[key], normalized_record)
        else:
            deduped[key] = normalized_record
    return sorted(
        deduped.values(),
        key=lambda item: (item.get("collection", ""), item.get("stock_id", "")),
    )


def _format_reagent_pairs(reagent_records: list[dict[str, str]]) -> str:
    """Format canonical reagent pairs for sheet export."""
    return "; ".join(
        f"({record.get('stock_id', '')}, {record.get('collection', '')})"
        for record in reagent_records or []
        if record.get("stock_id") and record.get("collection")
    )


__all__ = [
    "ReferenceReagent",
    "ReagentExtraction",
    "extract_reference_reagents",
    "_format_reagent_pairs",
    "_merge_deduplicated_reagents",
    "_normalize_collection_name",
    "_normalize_reagent_record",
    "_normalize_stock_id",
]
