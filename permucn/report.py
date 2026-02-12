"""Reporting helpers."""

from __future__ import annotations

from typing import List


def family_fieldnames(mode: str) -> List[str]:
    base = [
        "family_id",
        "mode",
        "direction",
        "include_trait_loss",
        "n_fg_01",
        "n_fg_10",
        "stat_obs",
        "p_empirical",
        "q_bh",
        "n_perm_used",
        "refined",
        "status",
    ]

    if mode == "binary":
        base.extend(
            [
                "fg_concordant_count",
                "fg_total",
                "bg_concordant_count",
                "bg_total",
                "fg_concordance_rate",
                "bg_concordance_rate",
            ]
        )
    elif mode == "rate":
        base.extend(
            [
                "fg_mean_signed_rate",
                "bg_mean_signed_rate",
                "fg_median_signed_rate",
                "bg_median_signed_rate",
            ]
        )

    return base
