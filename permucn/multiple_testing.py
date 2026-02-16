"""Multiple-testing correction utilities."""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence


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


def tarone_screen_min_pvalues(
    min_pvalues: Sequence[Optional[float]],
    alpha: float,
) -> Dict[str, object]:
    """Compute Tarone testability mask and Bonferroni denominator.

    Uses the discrete Bonferroni screening rule:
    find the smallest ``k`` such that ``m(k) <= k``, where
    ``m(k) = #{i | p_min_i <= alpha / k}``.
    """
    if alpha <= 0.0 or alpha >= 1.0:
        raise ValueError("alpha must be in (0, 1)")

    valid_indices = [i for i, p in enumerate(min_pvalues) if p is not None]
    m_total = len(valid_indices)
    testable = [False] * len(min_pvalues)
    if m_total == 0:
        return {
            "m_total": 0,
            "bonferroni_denom": 0,
            "threshold": None,
            "testable_mask": testable,
            "m_testable": 0,
        }

    denom = m_total
    eps = 1e-15
    for k in range(1, m_total + 1):
        cutoff = alpha / k
        m_k = sum(1 for i in valid_indices if float(min_pvalues[i]) <= cutoff + eps)
        if m_k <= k:
            denom = k
            break

    threshold = alpha / denom
    for i in valid_indices:
        if float(min_pvalues[i]) <= threshold + eps:
            testable[i] = True

    return {
        "m_total": m_total,
        "bonferroni_denom": denom,
        "threshold": threshold,
        "testable_mask": testable,
        "m_testable": sum(1 for v in testable if v),
    }


def bonferroni_adjust_selected(
    pvalues: Sequence[Optional[float]],
    selected_mask: Sequence[bool],
    denom: int,
) -> List[Optional[float]]:
    """Bonferroni-adjust selected hypotheses; unselected entries remain None."""
    if len(pvalues) != len(selected_mask):
        raise ValueError("pvalues and selected_mask must have the same length")
    if denom < 0:
        raise ValueError("denom must be >= 0")

    out: List[Optional[float]] = [None] * len(pvalues)
    if denom == 0:
        return out

    for i, (p, selected) in enumerate(zip(pvalues, selected_mask)):
        if p is None or not selected:
            continue
        q = float(p) * float(denom)
        out[i] = 1.0 if q > 1.0 else q
    return out
