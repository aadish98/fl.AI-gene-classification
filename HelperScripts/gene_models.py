"""Species-neutral models for gene classification pipelines."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class GeneRecord:
    """Normalized gene identity used by the classifier core."""

    species: str
    gene_id: str
    symbol: str
    authority_id: str = ""
    synonyms: set[str] = field(default_factory=set)
    trace: dict[str, Any] = field(default_factory=dict)

    @property
    def display_id(self) -> str:
        return self.authority_id or self.gene_id


@dataclass
class GeneCatalog:
    """Lookup tables for one species."""

    species: str
    genes: dict[str, GeneRecord] = field(default_factory=dict)
    symbol_to_gene_id: dict[str, str] = field(default_factory=dict)

    def get(self, gene_id: str) -> GeneRecord | None:
        return self.genes.get(str(gene_id or "").strip())

    def from_symbol(self, symbol: str) -> GeneRecord | None:
        gene_id = self.symbol_to_gene_id.get(str(symbol or "").strip())
        return self.get(gene_id) if gene_id else None


@dataclass
class SourceHit:
    """A source label attached to one gene/reference association."""

    source_key: str
    source_label: str
    priority: int = 0


@dataclass
class ReferenceCandidate:
    """One candidate literature reference for one gene."""

    gene_id: str
    paper_id: str
    source_key: str
    source_label: str
    pmid: str = ""
    pmcid: str = ""
    doi: str = ""
    year: int = 0
    snippet: str = ""
    source_labels: set[str] = field(default_factory=set)
    source_hits: list[SourceHit] = field(default_factory=list)

    def __post_init__(self):
        self.gene_id = str(self.gene_id or "").strip()
        self.paper_id = str(self.paper_id or "").strip()
        self.source_key = str(self.source_key or "").strip()
        self.source_label = str(self.source_label or "").strip() or self.source_key
        self.pmid = str(self.pmid or "").strip()
        self.pmcid = str(self.pmcid or "").strip().upper()
        if self.source_label:
            self.source_labels.add(self.source_label)
        if self.source_key or self.source_label:
            self.source_hits.append(SourceHit(self.source_key, self.source_label))


@dataclass
class GeneClassificationResult:
    """Classification result and reference bookkeeping for one gene."""

    gene: GeneRecord
    category: str = "None"
    confidence: int = 0
    rationale: str = ""
    classified_by: str = "Self"
    pmcids: set[str] = field(default_factory=set)
    supporting_refs: set[str] = field(default_factory=set)
    full_text_refs: set[str] = field(default_factory=set)
