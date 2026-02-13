"""Permutation cache load/save utilities."""

from __future__ import annotations

import gzip
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional

from .permutation import PermutationCache
from .tree import CanonicalTree


CACHE_VERSION = 1


@dataclass(frozen=True)
class CacheMatchSpec:
    tree_fingerprint: str
    include_trait_loss: bool
    fg_01_mask: int
    fg_10_mask: int


def tree_fingerprint(tree: CanonicalTree) -> str:
    """Compute a stable fingerprint for branch ordering/identity."""
    payload = "\n".join(tree.branch_key_by_index).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def make_cache_spec(
    tree: CanonicalTree,
    include_trait_loss: bool,
    fg_01_mask: int,
    fg_10_mask: int,
) -> CacheMatchSpec:
    return CacheMatchSpec(
        tree_fingerprint=tree_fingerprint(tree),
        include_trait_loss=include_trait_loss,
        fg_01_mask=fg_01_mask,
        fg_10_mask=fg_10_mask,
    )


def load_cache_bundle(path: str | Path) -> Dict[str, Any]:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    if str(path).endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return json.load(handle)

    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def save_cache_bundle(path: str | Path, bundle: Dict[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if str(path).endswith(".gz"):
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            json.dump(bundle, handle, ensure_ascii=True, sort_keys=True)
            handle.write("\n")
        return

    with path.open("w", encoding="utf-8") as handle:
        json.dump(bundle, handle, ensure_ascii=True, sort_keys=True)
        handle.write("\n")


def is_bundle_compatible(bundle: Dict[str, Any], spec: CacheMatchSpec) -> bool:
    if int(bundle.get("version", -1)) != CACHE_VERSION:
        return False
    if bundle.get("tree_fingerprint") != spec.tree_fingerprint:
        return False
    if bool(bundle.get("include_trait_loss", False)) != bool(spec.include_trait_loss):
        return False
    if bundle.get("fg_01_mask_hex") != _mask_to_hex(spec.fg_01_mask):
        return False
    if bundle.get("fg_10_mask_hex") != _mask_to_hex(spec.fg_10_mask):
        return False
    return True


def empty_bundle(spec: CacheMatchSpec) -> Dict[str, Any]:
    return {
        "version": CACHE_VERSION,
        "tree_fingerprint": spec.tree_fingerprint,
        "include_trait_loss": spec.include_trait_loss,
        "fg_01_mask_hex": _mask_to_hex(spec.fg_01_mask),
        "fg_10_mask_hex": _mask_to_hex(spec.fg_10_mask),
        "initial": None,
        "refine": None,
    }


def get_stage_cache(bundle: Dict[str, Any], stage: str, n_perm_required: int) -> Optional[PermutationCache]:
    raw = bundle.get(stage)
    if not isinstance(raw, dict):
        return None

    try:
        n_perm = int(raw.get("n_perm", 0))
    except (TypeError, ValueError):
        return None
    if n_perm < n_perm_required:
        return None

    masks_01_hex = raw.get("masks_01_hex", [])
    masks_10_hex = raw.get("masks_10_hex", [])
    if not isinstance(masks_01_hex, list) or not isinstance(masks_10_hex, list):
        return None
    if len(masks_01_hex) < n_perm_required or len(masks_10_hex) < n_perm_required:
        return None

    try:
        masks_01 = [_hex_to_mask(h) for h in masks_01_hex[:n_perm_required]]
        masks_10 = [_hex_to_mask(h) for h in masks_10_hex[:n_perm_required]]
    except (TypeError, ValueError):
        return None

    try:
        total_attempts = int(raw.get("total_attempts", 0))
    except (TypeError, ValueError):
        total_attempts = 0
    try:
        total_restarts = int(raw.get("total_restarts", 0))
    except (TypeError, ValueError):
        total_restarts = 0

    return PermutationCache(
        masks_01=masks_01,
        masks_10=masks_10,
        total_attempts=total_attempts,
        total_restarts=total_restarts,
    )


def put_stage_cache(bundle: Dict[str, Any], stage: str, cache: PermutationCache) -> None:
    bundle[stage] = {
        "n_perm": len(cache.masks_01),
        "masks_01_hex": [_mask_to_hex(m) for m in cache.masks_01],
        "masks_10_hex": [_mask_to_hex(m) for m in cache.masks_10],
        "total_attempts": cache.total_attempts,
        "total_restarts": cache.total_restarts,
    }


def _mask_to_hex(mask: int) -> str:
    return format(mask, "x")


def _hex_to_mask(text: str) -> int:
    return int(text, 16) if text else 0
