"""Foreground branch permutation generation with topology constraints."""

from __future__ import annotations

import multiprocessing as mp
import os
import random
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import dataclass
from itertools import repeat
from typing import Dict, List, Sequence, Tuple

from .tree import CanonicalTree


@dataclass(frozen=True)
class PermutationCache:
    masks_01: List[int]
    masks_10: List[int]
    total_attempts: int
    total_restarts: int


class _PermutationSampler:
    def __init__(
        self,
        tree: CanonicalTree,
        obs_mask_01: int,
        obs_mask_10: int,
        include_trait_loss: bool,
        max_perm_attempts: int,
        max_set_attempts: int,
    ) -> None:
        self.tree = tree
        self.include_trait_loss = include_trait_loss
        self.max_perm_attempts = max_perm_attempts
        self.max_set_attempts = max_set_attempts

        self.obs_mask_01 = obs_mask_01
        self.obs_mask_10 = obs_mask_10 if include_trait_loss else 0

        self.bin_counts_01 = _bin_counts(self.obs_mask_01, tree.clade_bin_by_branch_index)
        self.bin_counts_10 = _bin_counts(self.obs_mask_10, tree.clade_bin_by_branch_index)

        self.candidates_by_bin = [[] for _ in range(8)]
        for idx, b in enumerate(tree.clade_bin_by_branch_index):
            self.candidates_by_bin[b].append(idx)

    def generate_one(self, seed: int) -> Tuple[int, int, int, int]:
        attempts = 0
        restarts = 0
        rng = random.Random(seed)

        for _ in range(self.max_perm_attempts):
            attempts += 1
            try:
                m01 = self._sample_set(self.bin_counts_01, allowed_mask=self.tree.all_mask, rng=rng)

                if not self.include_trait_loss or self.obs_mask_10 == 0:
                    m10 = 0
                else:
                    if m01 != 0:
                        allowed_10 = _descendants_of_mask(m01, self.tree.desc_mask_by_branch_index, strict=True)
                        allowed_10 &= ~m01
                        # If impossible under descendant constraint, fallback to independent 1->0 sampling.
                        if not _has_bin_capacity(allowed_10, self.bin_counts_10, self.candidates_by_bin):
                            allowed_10 = self.tree.all_mask
                    else:
                        # Edge case policy: n_fg_01 == 0 and n_fg_10 > 0
                        allowed_10 = self.tree.all_mask

                    m10 = self._sample_set(self.bin_counts_10, allowed_mask=allowed_10, rng=rng)

                return m01, m10, attempts, restarts
            except RuntimeError:
                restarts += 1
                continue

        raise RuntimeError(
            "Failed to generate a valid permutation under constraints. "
            "Try reducing n-perm or relaxing constraints."
        )

    def _sample_set(self, bin_counts: Sequence[int], allowed_mask: int, rng: random.Random) -> int:
        target_total = sum(bin_counts)
        if target_total == 0:
            return 0

        if not _has_bin_capacity(allowed_mask, bin_counts, self.candidates_by_bin):
            raise RuntimeError("Insufficient candidates for requested bin counts")

        anc = self.tree.anc_mask_by_branch_index
        desc = self.tree.desc_mask_by_branch_index

        # Harder bins first: low capacity relative to demand.
        demand_bins = [b for b, c in enumerate(bin_counts) if c > 0]
        order = sorted(
            demand_bins,
            key=lambda b: (
                _allowed_count(self.candidates_by_bin[b], allowed_mask),
                -bin_counts[b],
            ),
        )

        for _ in range(self.max_set_attempts):
            selected_mask = 0
            ok = True

            for b in order:
                need = bin_counts[b]
                if need <= 0:
                    continue

                pool = [idx for idx in self.candidates_by_bin[b] if (allowed_mask >> idx) & 1]
                rng.shuffle(pool)

                picked = 0
                # Greedy with random traversal.
                for idx in pool:
                    if picked >= need:
                        break
                    conflict = ((anc[idx] | desc[idx]) & selected_mask) != 0
                    if conflict:
                        continue
                    selected_mask |= 1 << idx
                    picked += 1

                if picked < need:
                    ok = False
                    break

            if ok and selected_mask.bit_count() == target_total:
                return selected_mask

        raise RuntimeError("Could not sample a valid set for requested bin composition")


class PermutationGenerator:
    def __init__(
        self,
        tree: CanonicalTree,
        obs_mask_01: int,
        obs_mask_10: int,
        include_trait_loss: bool,
        seed: int | None = None,
        max_perm_attempts: int = 200,
        max_set_attempts: int = 200,
    ) -> None:
        self.tree = tree
        self.include_trait_loss = include_trait_loss
        self.seed = seed
        self.max_perm_attempts = max_perm_attempts
        self.max_set_attempts = max_set_attempts
        self.obs_mask_01 = obs_mask_01
        self.obs_mask_10 = obs_mask_10 if include_trait_loss else 0

    def generate(self, n_perm: int, jobs: int = 1) -> PermutationCache:
        if n_perm <= 0:
            raise ValueError("n_perm must be > 0")

        jobs_eff = _effective_jobs(jobs)
        jobs_eff = min(jobs_eff, n_perm)

        base_seed = self.seed if self.seed is not None else random.SystemRandom().randrange(1 << 63)
        payload = self._payload()
        indices = list(range(n_perm))

        if jobs_eff <= 1:
            rows, total_attempts, total_restarts = _generate_indices_chunk(payload, indices, base_seed)
            return _cache_from_rows(rows, n_perm, total_attempts, total_restarts)

        chunks = _split_indices(indices, jobs_eff)
        try:
            rows, total_attempts, total_restarts = _generate_parallel_chunks(payload, chunks, base_seed, jobs_eff)
        except (PermissionError, OSError):
            rows, total_attempts, total_restarts = _generate_thread_chunks(payload, chunks, base_seed, jobs_eff)

        return _cache_from_rows(rows, n_perm, total_attempts, total_restarts)

    def _payload(self) -> Dict[str, object]:
        return {
            "tree": self.tree,
            "obs_mask_01": self.obs_mask_01,
            "obs_mask_10": self.obs_mask_10,
            "include_trait_loss": self.include_trait_loss,
            "max_perm_attempts": self.max_perm_attempts,
            "max_set_attempts": self.max_set_attempts,
        }


def _effective_jobs(jobs: int) -> int:
    if jobs < 0:
        raise ValueError("jobs must be >= 0")
    if jobs == 0:
        return os.cpu_count() or 1
    return jobs


def _split_indices(indices: List[int], jobs: int) -> List[List[int]]:
    chunk_size = max(1, len(indices) // jobs)
    chunks: List[List[int]] = []
    i = 0
    while i < len(indices):
        chunks.append(indices[i : i + chunk_size])
        i += chunk_size
    return chunks


def _seed_for_index(base_seed: int, idx: int) -> int:
    # 64-bit mix so each permutation index maps to a stable seed.
    return (base_seed + (idx + 1) * 0x9E3779B97F4A7C15) & ((1 << 64) - 1)


def _build_sampler(payload: Dict[str, object]) -> _PermutationSampler:
    tree = payload["tree"]
    if not isinstance(tree, CanonicalTree):
        raise RuntimeError("Invalid worker payload: tree")
    return _PermutationSampler(
        tree=tree,
        obs_mask_01=int(payload["obs_mask_01"]),
        obs_mask_10=int(payload["obs_mask_10"]),
        include_trait_loss=bool(payload["include_trait_loss"]),
        max_perm_attempts=int(payload["max_perm_attempts"]),
        max_set_attempts=int(payload["max_set_attempts"]),
    )


def _generate_indices_chunk(
    payload: Dict[str, object],
    indices: List[int],
    base_seed: int,
) -> Tuple[List[Tuple[int, int, int]], int, int]:
    sampler = _build_sampler(payload)

    rows: List[Tuple[int, int, int]] = []
    total_attempts = 0
    total_restarts = 0
    for idx in indices:
        m01, m10, attempts, restarts = sampler.generate_one(_seed_for_index(base_seed, idx))
        rows.append((idx, m01, m10))
        total_attempts += attempts
        total_restarts += restarts
    return rows, total_attempts, total_restarts


def _generate_parallel_chunks(
    payload: Dict[str, object],
    chunks: List[List[int]],
    base_seed: int,
    jobs: int,
) -> Tuple[List[Tuple[int, int, int]], int, int]:
    out_rows: List[Tuple[int, int, int]] = []
    total_attempts = 0
    total_restarts = 0

    if os.name == "posix":
        mp_ctx = mp.get_context("fork")
        with ProcessPoolExecutor(max_workers=jobs, mp_context=mp_ctx) as executor:
            for rows, attempts, restarts in executor.map(
                _generate_indices_chunk, repeat(payload), chunks, repeat(base_seed)
            ):
                out_rows.extend(rows)
                total_attempts += attempts
                total_restarts += restarts
        return out_rows, total_attempts, total_restarts

    with ProcessPoolExecutor(max_workers=jobs) as executor:
        for rows, attempts, restarts in executor.map(
            _generate_indices_chunk, repeat(payload), chunks, repeat(base_seed)
        ):
            out_rows.extend(rows)
            total_attempts += attempts
            total_restarts += restarts
    return out_rows, total_attempts, total_restarts


def _generate_thread_chunks(
    payload: Dict[str, object],
    chunks: List[List[int]],
    base_seed: int,
    jobs: int,
) -> Tuple[List[Tuple[int, int, int]], int, int]:
    out_rows: List[Tuple[int, int, int]] = []
    total_attempts = 0
    total_restarts = 0

    with ThreadPoolExecutor(max_workers=jobs) as executor:
        for rows, attempts, restarts in executor.map(
            _generate_indices_chunk, repeat(payload), chunks, repeat(base_seed)
        ):
            out_rows.extend(rows)
            total_attempts += attempts
            total_restarts += restarts
    return out_rows, total_attempts, total_restarts


def _cache_from_rows(
    rows: List[Tuple[int, int, int]],
    n_perm: int,
    total_attempts: int,
    total_restarts: int,
) -> PermutationCache:
    masks_01 = [0] * n_perm
    masks_10 = [0] * n_perm
    for idx, m01, m10 in rows:
        masks_01[idx] = m01
        masks_10[idx] = m10
    return PermutationCache(
        masks_01=masks_01,
        masks_10=masks_10,
        total_attempts=total_attempts,
        total_restarts=total_restarts,
    )


def _bin_counts(mask: int, bin_by_idx: Sequence[int]) -> List[int]:
    out = [0] * 8
    m = mask
    while m:
        lsb = m & -m
        idx = lsb.bit_length() - 1
        out[bin_by_idx[idx]] += 1
        m ^= lsb
    return out


def _descendants_of_mask(mask: int, desc_masks: Sequence[int], strict: bool = True) -> int:
    out = 0
    m = mask
    while m:
        lsb = m & -m
        idx = lsb.bit_length() - 1
        out |= desc_masks[idx]
        m ^= lsb
    if strict:
        out &= ~mask
    return out


def _allowed_count(pool: Sequence[int], allowed_mask: int) -> int:
    return sum(1 for idx in pool if (allowed_mask >> idx) & 1)


def _has_bin_capacity(
    allowed_mask: int,
    bin_counts: Sequence[int],
    candidates_by_bin: Sequence[Sequence[int]],
) -> bool:
    for b, need in enumerate(bin_counts):
        if need <= 0:
            continue
        capacity = _allowed_count(candidates_by_bin[b], allowed_mask)
        if capacity < need:
            return False
    return True

