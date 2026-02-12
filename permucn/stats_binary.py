"""Binary-mode association statistics."""

from __future__ import annotations

from typing import Iterable, List, Sequence, Tuple


def sign_masks(deltas: Sequence[int]) -> Tuple[int, int]:
    """Return bitmasks for positive and negative deltas."""
    pos = 0
    neg = 0
    for i, d in enumerate(deltas):
        if d > 0:
            pos |= 1 << i
        elif d < 0:
            neg |= 1 << i
    return pos, neg


def observed_binary_stat(
    pos_mask: int,
    neg_mask: int,
    fg_01_mask: int,
    fg_10_mask: int,
    direction_sign: int,
    sig_mask: int | None = None,
) -> int:
    """Observed binary statistic using shared sign convention."""
    use_pos = pos_mask if sig_mask is None else (pos_mask & sig_mask)
    use_neg = neg_mask if sig_mask is None else (neg_mask & sig_mask)

    if direction_sign > 0:
        return (fg_01_mask & use_pos).bit_count() + (fg_10_mask & use_neg).bit_count()
    return (fg_01_mask & use_neg).bit_count() + (fg_10_mask & use_pos).bit_count()


def permutation_binary_stats(
    perm_masks_01: Sequence[int],
    perm_masks_10: Sequence[int],
    pos_mask: int,
    neg_mask: int,
    direction_sign: int,
    sig_mask: int | None = None,
) -> List[int]:
    use_pos = pos_mask if sig_mask is None else (pos_mask & sig_mask)
    use_neg = neg_mask if sig_mask is None else (neg_mask & sig_mask)

    out: List[int] = []
    if direction_sign > 0:
        for m01, m10 in zip(perm_masks_01, perm_masks_10):
            s = (m01 & use_pos).bit_count() + (m10 & use_neg).bit_count()
            out.append(s)
    else:
        for m01, m10 in zip(perm_masks_01, perm_masks_10):
            s = (m01 & use_neg).bit_count() + (m10 & use_pos).bit_count()
            out.append(s)
    return out


def binary_summary(
    pos_mask: int,
    neg_mask: int,
    fg_01_mask: int,
    fg_10_mask: int,
    direction_sign: int,
    all_mask: int,
    sig_mask: int | None = None,
) -> dict:
    fg_mask = fg_01_mask | fg_10_mask
    bg_mask = all_mask & ~fg_mask

    fg_conc = observed_binary_stat(
        pos_mask=pos_mask,
        neg_mask=neg_mask,
        fg_01_mask=fg_01_mask,
        fg_10_mask=fg_10_mask,
        direction_sign=direction_sign,
        sig_mask=sig_mask,
    )

    use_pos = pos_mask if sig_mask is None else (pos_mask & sig_mask)
    use_neg = neg_mask if sig_mask is None else (neg_mask & sig_mask)

    if direction_sign > 0:
        bg_conc = (bg_mask & use_pos).bit_count()
    else:
        bg_conc = (bg_mask & use_neg).bit_count()

    fg_total = fg_mask.bit_count()
    bg_total = bg_mask.bit_count()

    return {
        "fg_concordant_count": fg_conc,
        "fg_total": fg_total,
        "bg_concordant_count": bg_conc,
        "bg_total": bg_total,
        "fg_concordance_rate": (fg_conc / fg_total) if fg_total else None,
        "bg_concordance_rate": (bg_conc / bg_total) if bg_total else None,
    }


def empirical_pvalue_one_sided(obs: float, perm_stats: Sequence[float]) -> float:
    k = sum(1 for s in perm_stats if s >= obs)
    n = len(perm_stats)
    return (k + 1.0) / (n + 1.0)
