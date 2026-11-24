"""Resolve a FlyBase symbol shared by several genes using publication counts.

The counts come from the local index built by build_synonym_evidence.py (from the
FlyBase Chado XML gene file): for each gene, how many distinct publications used a
given symbol as one of that gene's symbols. No network access and no LLM is involved.
"""

from __future__ import annotations

import csv
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from HelperScripts.build_synonym_evidence import EVIDENCE_PATH, FLYBASE_RELEASE

# A candidate wins only with at least this many publications ...
MIN_WINNER_PUBS = 3
# ... and at least this many times the runner-up's publications.
WIN_RATIO = 2.0


@dataclass(frozen=True)
class Candidate:
    fbgn: str
    current: str
    n_pubs: int


class SynonymEvidence:
    def __init__(self, by_synonym: dict[str, list[Candidate]]):
        self._by_synonym = by_synonym

    def candidates(self, symbol: str) -> list[Candidate]:
        """Genes that use `symbol` as a symbol, most-cited first."""
        found = self._by_synonym.get(symbol, [])
        return sorted(found, key=lambda item: (-item.n_pubs, item.current, item.fbgn))


def load_evidence(path: Path = EVIDENCE_PATH) -> Optional[SynonymEvidence]:
    """Load the evidence index, or return None when it has not been built."""
    if not Path(path).exists():
        return None
    by_synonym: dict[str, list[Candidate]] = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_synonym.setdefault(row["synonym"], []).append(
                Candidate(row["fbgn"], row["current_symbol"], int(row["n_pubs"]))
            )
    return SynonymEvidence(by_synonym)


def pick_winner(candidates: list[Candidate]) -> Optional[Candidate]:
    """Return the clear leader among candidates (sorted most-cited first), else None."""
    if not candidates:
        return None
    leader = candidates[0]
    if leader.n_pubs < MIN_WINNER_PUBS:
        return None
    if len(candidates) == 1:
        return leader
    runner_up = candidates[1]
    if leader.n_pubs >= WIN_RATIO * runner_up.n_pubs:
        return leader
    return None


def describe_counts(candidates: list[Candidate]) -> str:
    return ", ".join(f"{item.current} {item.n_pubs}" for item in candidates)


def describe_candidates(candidates: list[Candidate]) -> str:
    return ", ".join(f"{item.current} {item.fbgn} ({item.n_pubs})" for item in candidates)


__all__ = [
    "Candidate",
    "FLYBASE_RELEASE",
    "MIN_WINNER_PUBS",
    "SynonymEvidence",
    "WIN_RATIO",
    "describe_candidates",
    "describe_counts",
    "load_evidence",
    "pick_winner",
]
