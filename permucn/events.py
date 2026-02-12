"""Per-family branch event helpers."""

from __future__ import annotations

import math
from typing import Iterable, List, Optional, Sequence


def build_significance_mask(prob_vec: Sequence[float] | None, alpha: float) -> int:
    """Return bitmask of branches significant at alpha.

    If prob_vec is None, returns 0 (no significant branches).
    """
    if prob_vec is None:
        return 0
    mask = 0
    for i, p in enumerate(prob_vec):
        if math.isnan(p):
            continue
        if p < alpha:
            mask |= 1 << i
    return mask
