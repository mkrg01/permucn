"""Maximum-likelihood ancestral reconstruction for binary traits."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from .tree import CanonicalTree


@dataclass(frozen=True)
class ASRResult:
    q01: float
    q10: float
    log_likelihood: float
    posterior_by_node: List[Tuple[float, float]]
    hard_state_by_node: List[Optional[int]]
    fg_01_mask: int
    fg_10_mask: int
    n_fg_01: int
    n_fg_10: int


def run_trait_asr_ml(
    tree: CanonicalTree,
    species_to_state: Dict[str, int],
    posterior_hi: float = 0.6,
    posterior_lo: float = 0.4,
) -> ASRResult:
    """Infer ancestral states and foreground transition branches under an ML Mk2 model."""
    _validate_species_match(tree, species_to_state)
    _validate_branch_lengths(tree)

    tip_state_by_node: List[Optional[int]] = [None] * len(tree.labels)
    for node_id, species in enumerate(tree.tip_species_by_node):
        if species is None:
            continue
        tip_state_by_node[node_id] = species_to_state[species]

    q01, q10, best_ll = _fit_rates_ml(tree, tip_state_by_node)
    ll, upward, edge_msg, posterior = _evaluate_model(tree, tip_state_by_node, q01, q10)

    hard: List[Optional[int]] = [None] * len(tree.labels)
    for node_id, (p0, p1) in enumerate(posterior):
        if p1 >= posterior_hi:
            hard[node_id] = 1
        elif p1 <= posterior_lo:
            hard[node_id] = 0
        else:
            hard[node_id] = None

    fg_01_mask = 0
    fg_10_mask = 0
    for bidx, child in enumerate(tree.node_by_branch_index):
        parent = tree.parent_by_node[child]
        if parent is None:
            continue
        p = hard[parent]
        c = hard[child]
        if p is None or c is None:
            continue
        if p == 0 and c == 1:
            fg_01_mask |= 1 << bidx
        elif p == 1 and c == 0:
            fg_10_mask |= 1 << bidx

    return ASRResult(
        q01=q01,
        q10=q10,
        log_likelihood=max(ll, best_ll),
        posterior_by_node=posterior,
        hard_state_by_node=hard,
        fg_01_mask=fg_01_mask,
        fg_10_mask=fg_10_mask,
        n_fg_01=fg_01_mask.bit_count(),
        n_fg_10=fg_10_mask.bit_count(),
    )


def _validate_species_match(tree: CanonicalTree, species_to_state: Dict[str, int]) -> None:
    tree_species = {s for s in tree.tip_species_by_node if s is not None}
    trait_species = set(species_to_state)

    missing = sorted(tree_species - trait_species)
    extra = sorted(trait_species - tree_species)
    if missing or extra:
        msg = ["Species mismatch between trait table and tree tips."]
        if missing:
            preview = ", ".join(missing[:8])
            msg.append(f"Missing in trait table ({len(missing)}): {preview}")
        if extra:
            preview = ", ".join(extra[:8])
            msg.append(f"Extra in trait table ({len(extra)}): {preview}")
        raise ValueError(" ".join(msg))


def _validate_branch_lengths(tree: CanonicalTree) -> None:
    bad: List[str] = []
    for bidx, node_id in enumerate(tree.node_by_branch_index):
        length = tree.branch_length_by_node[node_id]
        if not math.isfinite(length) or length < 0.0:
            bad.append(f"{tree.branch_key_by_index[bidx]}={length}")

    if bad:
        preview = ", ".join(bad[:6])
        raise ValueError(
            "Invalid branch lengths in tree for ASR; require finite and >= 0: "
            f"{preview}{' ...' if len(bad) > 6 else ''}"
        )


def _fit_rates_ml(
    tree: CanonicalTree,
    tip_state_by_node: Sequence[Optional[int]],
) -> Tuple[float, float, float]:
    """Simple robust ML fitting for q01/q10 via multi-stage grid search."""
    grid = [10 ** x for x in [-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0]]

    best_q01 = 1.0
    best_q10 = 1.0
    best_ll = float("-inf")

    for q01 in grid:
        for q10 in grid:
            ll, _, _, _ = _evaluate_model(tree, tip_state_by_node, q01, q10)
            if ll > best_ll:
                best_ll = ll
                best_q01 = q01
                best_q10 = q10

    # Local refinements around best in log-space.
    log_q01 = math.log10(best_q01)
    log_q10 = math.log10(best_q10)
    for radius in [0.5, 0.3, 0.15, 0.08]:
        points = [
            -radius,
            -radius * 0.5,
            0.0,
            radius * 0.5,
            radius,
        ]
        improved = False
        local_best = (best_ll, best_q01, best_q10)

        for dx in points:
            for dy in points:
                q01 = 10 ** (log_q01 + dx)
                q10 = 10 ** (log_q10 + dy)
                ll, _, _, _ = _evaluate_model(tree, tip_state_by_node, q01, q10)
                if ll > local_best[0]:
                    local_best = (ll, q01, q10)
                    improved = True

        if improved:
            best_ll, best_q01, best_q10 = local_best
            log_q01 = math.log10(best_q01)
            log_q10 = math.log10(best_q10)

    return best_q01, best_q10, best_ll


def _trans_probs(t: float, q01: float, q10: float) -> Tuple[float, float, float, float]:
    qsum = q01 + q10
    if qsum <= 0.0:
        return (1.0, 0.0, 0.0, 1.0)

    pi0 = q10 / qsum
    pi1 = q01 / qsum
    e = math.exp(-qsum * t)

    p00 = pi0 + pi1 * e
    p01 = pi1 - pi1 * e
    p10 = pi0 - pi0 * e
    p11 = pi1 + pi0 * e
    return (p00, p01, p10, p11)


def _evaluate_model(
    tree: CanonicalTree,
    tip_state_by_node: Sequence[Optional[int]],
    q01: float,
    q10: float,
) -> Tuple[float, List[Tuple[float, float]], Dict[Tuple[int, int], Tuple[float, float]], List[Tuple[float, float]]]:
    """Return log-likelihood and per-node posterior under given q rates."""
    n = len(tree.labels)
    root = tree.root

    postorder = _postorder(tree.root, tree.children_by_node)

    upward: List[Tuple[float, float]] = [(1.0, 1.0)] * n
    subtree_scale: List[float] = [0.0] * n
    edge_msg: Dict[Tuple[int, int], Tuple[float, float]] = {}

    for node in postorder:
        children = tree.children_by_node[node]
        if not children:
            obs = tip_state_by_node[node]
            if obs == 0:
                upward[node] = (1.0, 0.0)
            elif obs == 1:
                upward[node] = (0.0, 1.0)
            else:
                upward[node] = (1.0, 1.0)
            subtree_scale[node] = 0.0
            continue

        l0 = 1.0
        l1 = 1.0
        scale = 0.0
        for child in children:
            t = tree.branch_length_by_node[child]
            p00, p01, p10, p11 = _trans_probs(t, q01, q10)
            c0, c1 = upward[child]
            m0 = p00 * c0 + p01 * c1
            m1 = p10 * c0 + p11 * c1
            edge_msg[(node, child)] = (m0, m1)
            l0 *= m0
            l1 *= m1
            scale += subtree_scale[child]

        norm = max(l0, l1)
        if norm <= 0.0:
            return float("-inf"), upward, edge_msg, [(0.5, 0.5)] * n

        upward[node] = (l0 / norm, l1 / norm)
        subtree_scale[node] = scale + math.log(norm)

    qsum = q01 + q10
    pi0 = q10 / qsum
    pi1 = q01 / qsum
    r0, r1 = upward[root]
    root_lik = pi0 * r0 + pi1 * r1
    if root_lik <= 0.0:
        return float("-inf"), upward, edge_msg, [(0.5, 0.5)] * n
    loglik = math.log(root_lik) + subtree_scale[root]

    # Downward pass for marginals.
    preorder = _preorder(root, tree.children_by_node)
    down: List[Tuple[float, float]] = [(0.5, 0.5)] * n
    down[root] = _normalize((pi0, pi1))

    for parent in preorder:
        children = tree.children_by_node[parent]
        if not children:
            continue

        # Product of child messages by parent state.
        all0 = 1.0
        all1 = 1.0
        for child in children:
            m0, m1 = edge_msg[(parent, child)]
            all0 *= m0
            all1 *= m1

        for child in children:
            m0, m1 = edge_msg[(parent, child)]
            # Remove this child's contribution to keep only outside-of-child subtree.
            excl0 = all0 / m0 if m0 != 0.0 else _product_excluding(parent, child, edge_msg, tree, state=0)
            excl1 = all1 / m1 if m1 != 0.0 else _product_excluding(parent, child, edge_msg, tree, state=1)

            up0, up1 = down[parent]
            base0 = up0 * excl0
            base1 = up1 * excl1

            t = tree.branch_length_by_node[child]
            p00, p01, p10, p11 = _trans_probs(t, q01, q10)
            c0 = base0 * p00 + base1 * p10
            c1 = base0 * p01 + base1 * p11
            down[child] = _normalize((c0, c1))

    posterior: List[Tuple[float, float]] = [(0.5, 0.5)] * n
    for node in range(n):
        u0, u1 = down[node]
        l0, l1 = upward[node]
        posterior[node] = _normalize((u0 * l0, u1 * l1))

    return loglik, upward, edge_msg, posterior


def _product_excluding(
    parent: int,
    excluded_child: int,
    edge_msg: Dict[Tuple[int, int], Tuple[float, float]],
    tree: CanonicalTree,
    state: int,
) -> float:
    out = 1.0
    for child in tree.children_by_node[parent]:
        if child == excluded_child:
            continue
        m0, m1 = edge_msg[(parent, child)]
        out *= m0 if state == 0 else m1
    return out


def _normalize(v: Tuple[float, float]) -> Tuple[float, float]:
    a, b = v
    s = a + b
    if s <= 0.0:
        return (0.5, 0.5)
    return (a / s, b / s)


def _postorder(root: int, children_by_node: Sequence[Sequence[int]]) -> List[int]:
    out: List[int] = []

    def rec(n: int) -> None:
        for ch in children_by_node[n]:
            rec(ch)
        out.append(n)

    rec(root)
    return out


def _preorder(root: int, children_by_node: Sequence[Sequence[int]]) -> List[int]:
    out: List[int] = []

    def rec(n: int) -> None:
        out.append(n)
        for ch in children_by_node[n]:
            rec(ch)

    rec(root)
    return out
