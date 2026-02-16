"""Fisher exact-test utilities for binary concordance tables."""

from __future__ import annotations

import math


def fisher_exact_one_sided_from_counts(
    *,
    fg_concordant_count: int,
    fg_total: int,
    bg_concordant_count: int,
    bg_total: int,
) -> float:
    """Return one-sided Fisher exact p-value for foreground enrichment.

    The tested 2x2 table is:

    - row 1 (foreground): concordant / non-concordant
    - row 2 (background): concordant / non-concordant
    """
    _validate_2x2_counts(
        fg_concordant_count=fg_concordant_count,
        fg_total=fg_total,
        bg_concordant_count=bg_concordant_count,
        bg_total=bg_total,
    )
    n_fg = fg_total
    n_bg = bg_total
    n_concordant = fg_concordant_count + bg_concordant_count

    lower, upper = _support_bounds(n_fg=n_fg, n_bg=n_bg, n_concordant=n_concordant)
    obs = fg_concordant_count
    if obs < lower or obs > upper:
        raise ValueError("Observed concordant count is outside hypergeometric support")

    tail_logps = [_hypergeom_log_pmf(x=x, n_fg=n_fg, n_bg=n_bg, n_concordant=n_concordant) for x in range(obs, upper + 1)]
    return math.exp(_logsumexp(tail_logps))


def fisher_min_attainable_pvalue(
    *,
    fg_total: int,
    bg_total: int,
    total_concordant: int,
) -> float:
    """Return the minimal attainable one-sided p-value for fixed margins."""
    if fg_total < 0 or bg_total < 0 or total_concordant < 0:
        raise ValueError("Margins must be non-negative")
    if total_concordant > fg_total + bg_total:
        raise ValueError("total_concordant cannot exceed fg_total + bg_total")

    _, upper = _support_bounds(n_fg=fg_total, n_bg=bg_total, n_concordant=total_concordant)
    return math.exp(_hypergeom_log_pmf(x=upper, n_fg=fg_total, n_bg=bg_total, n_concordant=total_concordant))


def _validate_2x2_counts(
    *,
    fg_concordant_count: int,
    fg_total: int,
    bg_concordant_count: int,
    bg_total: int,
) -> None:
    if fg_total < 0 or bg_total < 0:
        raise ValueError("Totals must be non-negative")
    if fg_concordant_count < 0 or bg_concordant_count < 0:
        raise ValueError("Concordant counts must be non-negative")
    if fg_concordant_count > fg_total:
        raise ValueError("fg_concordant_count cannot exceed fg_total")
    if bg_concordant_count > bg_total:
        raise ValueError("bg_concordant_count cannot exceed bg_total")


def _support_bounds(*, n_fg: int, n_bg: int, n_concordant: int) -> tuple[int, int]:
    lower = max(0, n_concordant - n_bg)
    upper = min(n_fg, n_concordant)
    return lower, upper


def _hypergeom_log_pmf(*, x: int, n_fg: int, n_bg: int, n_concordant: int) -> float:
    total = n_fg + n_bg
    if x < 0 or x > n_fg:
        return float("-inf")
    if n_concordant - x < 0 or n_concordant - x > n_bg:
        return float("-inf")

    return (
        _log_choose(n_concordant, x)
        + _log_choose(total - n_concordant, n_fg - x)
        - _log_choose(total, n_fg)
    )


def _log_choose(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1.0) - math.lgamma(k + 1.0) - math.lgamma(n - k + 1.0)


def _logsumexp(values: list[float]) -> float:
    if not values:
        return float("-inf")
    vmax = max(values)
    if math.isinf(vmax):
        return vmax
    return vmax + math.log(sum(math.exp(v - vmax) for v in values))
