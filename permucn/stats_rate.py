"""Rate-mode association statistics."""

from __future__ import annotations

import statistics
from typing import Iterable, List, Sequence


def build_rates(deltas: Sequence[int], branch_lengths: Sequence[float]) -> List[float]:
    if len(deltas) != len(branch_lengths):
        raise ValueError("Delta and branch-length vectors must have the same length")

    out: List[float] = []
    for i, (d, l) in enumerate(zip(deltas, branch_lengths)):
        if l <= 0.0:
            raise ValueError(f"Non-positive branch length at branch index {i}: {l}")
        out.append(float(d) / float(l))
    return out


def observed_rate_stat(
    rates: Sequence[float],
    fg_01_mask: int,
    fg_10_mask: int,
    direction_sign: int,
) -> float:
    n = fg_01_mask.bit_count() + fg_10_mask.bit_count()
    if n == 0:
        return float("nan")

    s01 = _sum_mask(rates, fg_01_mask)
    s10 = _sum_mask(rates, fg_10_mask)
    return direction_sign * (s01 - s10) / n


def permutation_rate_stats(
    rates: Sequence[float],
    perm_masks_01: Sequence[int],
    perm_masks_10: Sequence[int],
    direction_sign: int,
) -> List[float]:
    out: List[float] = []
    for m01, m10 in zip(perm_masks_01, perm_masks_10):
        n = m01.bit_count() + m10.bit_count()
        if n == 0:
            out.append(float("nan"))
            continue
        s01 = _sum_mask(rates, m01)
        s10 = _sum_mask(rates, m10)
        out.append(direction_sign * (s01 - s10) / n)
    return out


def rate_summary(
    rates: Sequence[float],
    fg_01_mask: int,
    fg_10_mask: int,
    direction_sign: int,
    all_mask: int,
) -> dict:
    fg_vals = _signed_fg_values(rates, fg_01_mask, fg_10_mask, direction_sign)

    fg_mask = fg_01_mask | fg_10_mask
    bg_mask = all_mask & ~fg_mask
    bg_vals = [direction_sign * rates[idx] for idx in _iter_mask_indices(bg_mask)]

    fg_mean = (sum(fg_vals) / len(fg_vals)) if fg_vals else None
    bg_mean = (sum(bg_vals) / len(bg_vals)) if bg_vals else None

    fg_median = statistics.median(fg_vals) if fg_vals else None
    bg_median = statistics.median(bg_vals) if bg_vals else None

    return {
        "fg_mean_signed_rate": fg_mean,
        "bg_mean_signed_rate": bg_mean,
        "fg_median_signed_rate": fg_median,
        "bg_median_signed_rate": bg_median,
    }


def _sum_mask(values: Sequence[float], mask: int) -> float:
    s = 0.0
    m = mask
    while m:
        lsb = m & -m
        idx = lsb.bit_length() - 1
        s += values[idx]
        m ^= lsb
    return s


def _iter_mask_indices(mask: int) -> Iterable[int]:
    m = mask
    while m:
        lsb = m & -m
        idx = lsb.bit_length() - 1
        yield idx
        m ^= lsb


def _signed_fg_values(
    rates: Sequence[float],
    fg_01_mask: int,
    fg_10_mask: int,
    direction_sign: int,
) -> List[float]:
    out: List[float] = []
    for idx in _iter_mask_indices(fg_01_mask):
        out.append(direction_sign * rates[idx])
    for idx in _iter_mask_indices(fg_10_mask):
        out.append(-direction_sign * rates[idx])
    return out
