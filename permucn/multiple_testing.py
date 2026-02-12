"""Multiple-testing correction utilities."""

from __future__ import annotations

from typing import List, Optional, Sequence


def bh_adjust_with_none(pvalues: Sequence[Optional[float]]) -> List[Optional[float]]:
    """Benjamini-Hochberg adjustment preserving None entries."""
    indexed = [(i, p) for i, p in enumerate(pvalues) if p is not None]
    m = len(indexed)
    out: List[Optional[float]] = [None] * len(pvalues)
    if m == 0:
        return out

    indexed.sort(key=lambda x: x[1])

    qvals = [0.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):
        idx, p = indexed[rank - 1]
        q = (p * m) / rank
        if q > 1.0:
            q = 1.0
        if q > prev:
            q = prev
        prev = q
        qvals[rank - 1] = q

    for (orig_idx, _), q in zip(indexed, qvals):
        out[orig_idx] = q
    return out
