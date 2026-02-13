"""Command-line entrypoint for permucn."""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional

from . import __version__
from .cache import (
    empty_bundle,
    get_stage_cache,
    is_bundle_compatible,
    load_cache_bundle,
    make_cache_spec,
    put_stage_cache,
    save_cache_bundle,
)
from .events import build_significance_mask
from .io import (
    load_change_matrix,
    load_probability_map,
    load_trait_table,
    write_json,
    write_tsv,
)
from .multiple_testing import bh_adjust_with_none
from .permutation import PermutationCache, PermutationGenerator
from .report import family_fieldnames
from .stats_binary import (
    binary_summary,
    empirical_pvalue_one_sided,
    observed_binary_stat,
    permutation_binary_stats,
    sign_masks,
)
from .stats_rate import build_rates, observed_rate_stat, permutation_rate_stats, rate_summary
from .trait_ml import run_trait_asr_ml
from .tree import load_canonical_tree
from .viz import generate_visual_outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="permucn",
        description="Permutation-based copy number / trait association testing",
    )

    parser.add_argument(
        "--cafe-dir",
        default=None,
        help="Directory containing CAFE output files",
    )
    parser.add_argument(
        "--trait-tsv",
        default=None,
        help="Trait TSV path",
    )
    parser.add_argument(
        "--trait-column",
        default=None,
        help="Trait column name in trait TSV (optional; auto-detected if omitted)",
    )
    parser.add_argument("--mode", choices=["binary", "rate"], default="binary")
    parser.add_argument("--direction", choices=["gain", "loss"], default="gain")

    incl = parser.add_mutually_exclusive_group()
    incl.add_argument("--include-trait-loss", dest="include_trait_loss", action="store_true")
    incl.add_argument("--no-include-trait-loss", dest="include_trait_loss", action="store_false")
    parser.set_defaults(include_trait_loss=True)

    parser.add_argument("--asr-method", choices=["ml"], default="ml")
    parser.add_argument("--asr-posterior-hi", type=float, default=0.6)
    parser.add_argument("--asr-posterior-lo", type=float, default=0.4)

    parser.add_argument("--cafe-significant-only", action="store_true")
    parser.add_argument("--cafe-alpha", type=float, default=0.05)

    parser.add_argument("--n-perm-initial", type=int, default=1000)
    parser.add_argument("--n-perm-refine", type=int, default=1000000)
    parser.add_argument("--refine-p-threshold", type=float, default=0.01)

    parser.add_argument(
        "--clade-bin-scheme",
        choices=["log2"],
        default="log2",
        help="Clade-size binning scheme",
    )
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of worker processes (1 = sequential, 0 = auto)",
    )
    parser.add_argument("--out-prefix", default="permucn_results")
    parser.add_argument(
        "--perm-cache",
        default=None,
        help="Optional JSON/JSON.GZ file to load/save permutation cache",
    )
    parser.add_argument(
        "--make-plots",
        action="store_true",
        help="Generate QQ/histogram PDF outputs (requires matplotlib)",
    )
    parser.add_argument(
        "--qvalue-threshold",
        type=float,
        default=0.05,
        help="Write families with q_bh <= threshold to <out-prefix>.top_hits.tsv",
    )
    parser.add_argument(
        "--pvalue-top-n",
        type=int,
        default=0,
        help=(
            "Write top N families ranked by smallest p_empirical to "
            "<out-prefix>.top_pvalues.tsv (0 disables)"
        ),
    )
    parser.add_argument(
        "--hist-bins",
        type=int,
        default=20,
        help="Number of bins for p-value histogram output",
    )

    return parser


def _required_paths(cafe_dir: Path) -> Dict[str, Path]:
    return {
        "change": cafe_dir / "Gamma_change.tab",
        "prob": cafe_dir / "Gamma_branch_probabilities.tab",
        "asr_tree": cafe_dir / "Gamma_asr.tre",
        "family_meta": cafe_dir / "Gamma_family_results.txt",
    }


def _validate_args(args: argparse.Namespace) -> None:
    if hasattr(args, "cafe_dir") and not args.cafe_dir:
        raise ValueError("--cafe-dir is required")
    if hasattr(args, "trait_tsv") and not args.trait_tsv:
        raise ValueError("--trait-tsv is required")

    if args.asr_posterior_lo < 0 or args.asr_posterior_hi > 1 or args.asr_posterior_lo >= args.asr_posterior_hi:
        raise ValueError("Invalid ASR posterior thresholds: require 0 <= lo < hi <= 1")

    if args.n_perm_initial <= 0:
        raise ValueError("--n-perm-initial must be > 0")
    if args.n_perm_refine <= 0:
        raise ValueError("--n-perm-refine must be > 0")

    if args.refine_p_threshold <= 0 or args.refine_p_threshold >= 1:
        raise ValueError("--refine-p-threshold must be in (0, 1)")

    if args.cafe_alpha <= 0 or args.cafe_alpha >= 1:
        raise ValueError("--cafe-alpha must be in (0, 1)")

    qvalue_threshold = getattr(args, "qvalue_threshold", 0.05)
    pvalue_top_n = getattr(args, "pvalue_top_n", 0)
    hist_bins = getattr(args, "hist_bins", 20)
    jobs = getattr(args, "jobs", 1)
    if qvalue_threshold < 0 or qvalue_threshold > 1:
        raise ValueError("--qvalue-threshold must be in [0, 1]")
    if pvalue_top_n < 0:
        raise ValueError("--pvalue-top-n must be >= 0")
    if hist_bins <= 0:
        raise ValueError("--hist-bins must be > 0")
    if jobs < 0:
        raise ValueError("--jobs must be >= 0")

    if args.mode == "rate" and args.cafe_significant_only:
        raise ValueError("--cafe-significant-only is valid only in binary mode")


def _family_result_base(
    family_id: str,
    args: argparse.Namespace,
    n_fg_01: int,
    n_fg_10: int,
) -> Dict[str, object]:
    return {
        "family_id": family_id,
        "mode": args.mode,
        "direction": args.direction,
        "include_trait_loss": args.include_trait_loss,
        "n_fg_01": n_fg_01,
        "n_fg_10": n_fg_10,
        "stat_obs": None,
        "p_empirical": None,
        "q_bh": None,
        "n_perm_used": 0,
        "refined": False,
        "status": "not_tested",
    }


def _log_progress(message: str) -> None:
    print(f"[permucn] {message}", file=sys.stderr, flush=True)


def _log_warning(message: str) -> None:
    print(f"[permucn][warning] {message}", file=sys.stderr, flush=True)


def _compute_family_binary(
    deltas: List[int],
    perm: PermutationCache,
    direction_sign: int,
    fg_01_mask: int,
    fg_10_mask: int,
    all_mask: int,
    cafe_sig_mask: int | None,
) -> Dict[str, object]:
    pos_mask, neg_mask = sign_masks(deltas)

    stat_obs = observed_binary_stat(
        pos_mask=pos_mask,
        neg_mask=neg_mask,
        fg_01_mask=fg_01_mask,
        fg_10_mask=fg_10_mask,
        direction_sign=direction_sign,
        sig_mask=cafe_sig_mask,
    )

    perm_stats = permutation_binary_stats(
        perm_masks_01=perm.masks_01,
        perm_masks_10=perm.masks_10,
        pos_mask=pos_mask,
        neg_mask=neg_mask,
        direction_sign=direction_sign,
        sig_mask=cafe_sig_mask,
    )

    p = empirical_pvalue_one_sided(stat_obs, perm_stats)
    summary = binary_summary(
        pos_mask=pos_mask,
        neg_mask=neg_mask,
        fg_01_mask=fg_01_mask,
        fg_10_mask=fg_10_mask,
        direction_sign=direction_sign,
        all_mask=all_mask,
        sig_mask=cafe_sig_mask,
    )

    out: Dict[str, object] = {
        "stat_obs": stat_obs,
        "p_empirical": p,
    }
    out.update(summary)
    return out


def _compute_family_rate(
    deltas: List[int],
    branch_lengths: List[float],
    perm: PermutationCache,
    direction_sign: int,
    fg_01_mask: int,
    fg_10_mask: int,
    all_mask: int,
) -> Dict[str, object]:
    rates = build_rates(deltas, branch_lengths)

    stat_obs = observed_rate_stat(
        rates=rates,
        fg_01_mask=fg_01_mask,
        fg_10_mask=fg_10_mask,
        direction_sign=direction_sign,
    )

    perm_stats = permutation_rate_stats(
        rates=rates,
        perm_masks_01=perm.masks_01,
        perm_masks_10=perm.masks_10,
        direction_sign=direction_sign,
    )

    p = empirical_pvalue_one_sided(stat_obs, perm_stats)
    summary = rate_summary(
        rates=rates,
        fg_01_mask=fg_01_mask,
        fg_10_mask=fg_10_mask,
        direction_sign=direction_sign,
        all_mask=all_mask,
    )

    out: Dict[str, object] = {
        "stat_obs": stat_obs,
        "p_empirical": p,
    }
    out.update(summary)
    return out


_WORKER_CONTEXT: Dict[str, object] | None = None


def _effective_jobs(raw_jobs: int) -> int:
    if raw_jobs < 0:
        raise ValueError("--jobs must be >= 0")
    if raw_jobs == 0:
        return os.cpu_count() or 1
    return raw_jobs


def _build_worker_context(
    *,
    mode: str,
    family_values: List[List[int]],
    perm: PermutationCache,
    direction_sign: int,
    fg_01_mask: int,
    fg_10_mask: int,
    all_mask: int,
    branch_lengths: List[float] | None,
    cafe_sig_masks: List[int] | None,
) -> Dict[str, object]:
    return {
        "mode": mode,
        "family_values": family_values,
        "perm": perm,
        "direction_sign": direction_sign,
        "fg_01_mask": fg_01_mask,
        "fg_10_mask": fg_10_mask,
        "all_mask": all_mask,
        "branch_lengths": branch_lengths,
        "cafe_sig_masks": cafe_sig_masks,
    }


def _init_worker(context: Dict[str, object]) -> None:
    global _WORKER_CONTEXT
    _WORKER_CONTEXT = context


def _evaluate_family_with_context(index: int, context: Dict[str, object]) -> tuple[int, Dict[str, object]]:
    mode = str(context["mode"])
    family_values = context["family_values"]
    perm = context["perm"]
    direction_sign = int(context["direction_sign"])
    fg_01_mask = int(context["fg_01_mask"])
    fg_10_mask = int(context["fg_10_mask"])
    all_mask = int(context["all_mask"])

    deltas = family_values[index]

    cafe_sig_mask = None
    cafe_sig_masks = context.get("cafe_sig_masks")
    if cafe_sig_masks is not None:
        cafe_sig_mask = cafe_sig_masks[index]

    if mode == "binary":
        out = _compute_family_binary(
            deltas=deltas,
            perm=perm,
            direction_sign=direction_sign,
            fg_01_mask=fg_01_mask,
            fg_10_mask=fg_10_mask,
            all_mask=all_mask,
            cafe_sig_mask=cafe_sig_mask,
        )
    else:
        branch_lengths = context["branch_lengths"]
        if branch_lengths is None:
            raise RuntimeError("Missing branch lengths in rate-mode worker context")
        out = _compute_family_rate(
            deltas=deltas,
            branch_lengths=branch_lengths,
            perm=perm,
            direction_sign=direction_sign,
            fg_01_mask=fg_01_mask,
            fg_10_mask=fg_10_mask,
            all_mask=all_mask,
        )

    return index, out


def _worker_evaluate_family(index: int) -> tuple[int, Dict[str, object]]:
    if _WORKER_CONTEXT is None:
        raise RuntimeError("Worker context is not initialized")
    return _evaluate_family_with_context(index, _WORKER_CONTEXT)


def _evaluate_families(
    *,
    indices: List[int],
    jobs: int,
    context: Dict[str, object],
) -> Dict[int, Dict[str, object]]:
    out: Dict[int, Dict[str, object]] = {}

    if jobs <= 1 or len(indices) <= 1:
        for index in indices:
            i, row = _evaluate_family_with_context(index, context)
            out[i] = row
        return out

    chunksize = max(1, len(indices) // (jobs * 4))
    try:
        if os.name == "posix":
            mp_ctx = mp.get_context("fork")
            with ProcessPoolExecutor(
                max_workers=jobs,
                mp_context=mp_ctx,
                initializer=_init_worker,
                initargs=(context,),
            ) as executor:
                for i, row in executor.map(_worker_evaluate_family, indices, chunksize=chunksize):
                    out[i] = row
            return out

        with ProcessPoolExecutor(
            max_workers=jobs,
            initializer=_init_worker,
            initargs=(context,),
        ) as executor:
            for i, row in executor.map(_worker_evaluate_family, indices, chunksize=chunksize):
                out[i] = row
        return out
    except (PermissionError, OSError):
        out.clear()
        with ThreadPoolExecutor(max_workers=jobs) as executor:
            for i, row in executor.map(lambda idx: _evaluate_family_with_context(idx, context), indices):
                out[i] = row

    return out


def run(args: argparse.Namespace) -> int:
    _log_progress("[1/8] Validating arguments and inputs")
    _validate_args(args)

    cafe_dir = Path(args.cafe_dir)
    trait_tsv = Path(args.trait_tsv)
    jobs = _effective_jobs(getattr(args, "jobs", 1))

    paths = _required_paths(cafe_dir)
    for key in ["change", "asr_tree"]:
        if not paths[key].exists():
            raise FileNotFoundError(f"Required input file missing: {paths[key]}")

    _log_progress("[2/8] Loading tree/trait data and running trait ASR")
    tree = load_canonical_tree(paths["asr_tree"])
    root_key = tree.branch_keys_by_node[tree.root]

    trait = load_trait_table(trait_tsv, args.trait_column)

    asr = run_trait_asr_ml(
        tree=tree,
        species_to_state=trait.species_to_state,
        posterior_hi=args.asr_posterior_hi,
        posterior_lo=args.asr_posterior_lo,
    )
    map_state_by_node = [1 if p1 >= 0.5 else 0 for _, p1 in asr.posterior_by_node]
    potential_transition_branch_keys: List[str] = []
    potential_01 = 0
    potential_10 = 0
    for bidx, child in enumerate(tree.node_by_branch_index):
        parent = tree.parent_by_node[child]
        if parent is None:
            continue

        hard_parent = asr.hard_state_by_node[parent]
        hard_child = asr.hard_state_by_node[child]
        if hard_parent is not None and hard_child is not None:
            continue

        parent_state = map_state_by_node[parent]
        child_state = map_state_by_node[child]
        if parent_state == child_state:
            continue

        potential_transition_branch_keys.append(tree.branch_key_by_index[bidx])
        if parent_state == 0:
            potential_01 += 1
        else:
            potential_10 += 1

    if potential_transition_branch_keys:
        preview = ", ".join(potential_transition_branch_keys[:8])
        _log_warning(
            "ASR posterior thresholding skipped potential phenotype-transition branches: "
            f"{len(potential_transition_branch_keys)} branch(es) would be transitions under "
            f"posterior>=0.5 binarization (0->1={potential_01}, 1->0={potential_10}; "
            f"branch keys: {preview})"
        )

    fg_01_mask = asr.fg_01_mask
    fg_10_mask = asr.fg_10_mask if args.include_trait_loss else 0
    n_fg_01 = fg_01_mask.bit_count()
    n_fg_10 = fg_10_mask.bit_count()
    fg_total = n_fg_01 + n_fg_10
    _log_progress(
        f"[2/8] Foreground branches detected: 0->1={n_fg_01}, 1->0={n_fg_10}, total={fg_total}"
    )

    _log_progress("[3/8] Preparing permutation cache and initial permutations")
    cache_spec = make_cache_spec(
        tree=tree,
        include_trait_loss=args.include_trait_loss,
        fg_01_mask=fg_01_mask,
        fg_10_mask=fg_10_mask,
    )
    cache_bundle = None
    cache_loaded = False
    initial_source = "generated"
    refine_source = "not_used"
    perm_initial = None
    if args.perm_cache:
        cache_path = Path(args.perm_cache)
        if cache_path.exists():
            loaded = load_cache_bundle(cache_path)
            if is_bundle_compatible(loaded, cache_spec):
                cache_bundle = loaded
                cache_loaded = True
            else:
                cache_bundle = empty_bundle(cache_spec)
        else:
            cache_bundle = empty_bundle(cache_spec)

    if fg_total == 0:
        initial_source = "skipped_no_foreground"
        refine_source = "skipped_no_foreground"
        _log_progress("[3/8] Skipping permutation generation because no foreground branches were found")
    else:
        # Initial permutations.
        perm_seed = args.seed
        if cache_bundle is not None:
            perm_initial = get_stage_cache(cache_bundle, "initial", args.n_perm_initial)
            if perm_initial is not None:
                initial_source = "cache"
        if perm_initial is None:
            _log_progress(f"[3/8] Generating initial permutations (n={args.n_perm_initial}, jobs={jobs})")
            perm_gen = PermutationGenerator(
                tree=tree,
                obs_mask_01=fg_01_mask,
                obs_mask_10=fg_10_mask,
                include_trait_loss=args.include_trait_loss,
                seed=perm_seed,
            )
            perm_initial = perm_gen.generate(args.n_perm_initial, jobs=jobs)
            initial_source = "generated"
            if cache_bundle is not None:
                put_stage_cache(cache_bundle, "initial", perm_initial)
        else:
            _log_progress(f"[3/8] Reusing initial permutations from cache (n={args.n_perm_initial})")

    # Family delta matrix.
    _log_progress("[4/8] Loading family change matrix")
    change = load_change_matrix(
        path=paths["change"],
        branch_to_index=tree.branch_index_by_key,
        ignored_branch_keys={root_key},
    )

    prob_map: Dict[str, List[float]] = {}
    if args.cafe_significant_only:
        if not paths["prob"].exists():
            raise FileNotFoundError(
                "--cafe-significant-only requires Gamma_branch_probabilities.tab, "
                f"but file is missing: {paths['prob']}"
            )
        _log_progress(f"[4/8] Loading branch probabilities for significance masking (alpha={args.cafe_alpha})")
        prob_map = load_probability_map(
            path=paths["prob"],
            branch_to_index=tree.branch_index_by_key,
            ignored_branch_keys={root_key},
        )

    direction_sign = 1 if args.direction == "gain" else -1

    branch_lengths = [tree.branch_length_by_node[n] for n in tree.node_by_branch_index]
    if args.mode == "rate":
        bad = [i for i, l in enumerate(branch_lengths) if l <= 0.0]
        if bad:
            preview = ", ".join(str(i) for i in bad[:8])
            raise ValueError(
                "Non-positive branch lengths found in canonical tree for rate mode; "
                f"branch indices: {preview}"
            )

    rows: List[Dict[str, object]] = [
        _family_result_base(fam_id, args, n_fg_01, n_fg_10) for fam_id in change.family_ids
    ]
    pvalues: List[Optional[float]] = [None] * len(rows)
    cafe_sig_masks: List[int] | None = None

    if fg_total == 0:
        _log_progress("[5/8] Skipping family tests because no valid foreground branches were found")
        for row in rows:
            row["status"] = "no_valid_foreground"
    else:
        _log_progress(
            f"[5/8] Running initial family tests for {len(rows)} families (n_perm={args.n_perm_initial})"
        )
        if perm_initial is None:
            raise RuntimeError("Initial permutations are unavailable")
        if args.cafe_significant_only:
            cafe_sig_masks = [build_significance_mask(prob_map.get(fam_id), args.cafe_alpha) for fam_id in change.family_ids]

        context_initial = _build_worker_context(
            mode=args.mode,
            family_values=change.values,
            perm=perm_initial,
            direction_sign=direction_sign,
            fg_01_mask=fg_01_mask,
            fg_10_mask=fg_10_mask,
            all_mask=tree.all_mask,
            branch_lengths=branch_lengths if args.mode == "rate" else None,
            cafe_sig_masks=cafe_sig_masks,
        )
        initial_indices = list(range(len(rows)))
        initial_results = _evaluate_families(indices=initial_indices, jobs=jobs, context=context_initial)

        for i, row in enumerate(rows):
            row.update(initial_results[i])
            row["n_perm_used"] = args.n_perm_initial
            row["status"] = "ok"
            pvalues[i] = row["p_empirical"]
        _log_progress("[5/8] Initial family tests completed")

    # Optional refinement stage.
    refined_indices = [
        i
        for i, row in enumerate(rows)
        if row["status"] == "ok"
        and row["p_empirical"] is not None
        and float(row["p_empirical"]) <= args.refine_p_threshold
        and args.n_perm_refine > args.n_perm_initial
    ]

    perm_refine = None
    if refined_indices:
        _log_progress(
            f"[6/8] Running refinement for {len(refined_indices)} families (n_perm={args.n_perm_refine})"
        )
        if cache_bundle is not None:
            perm_refine = get_stage_cache(cache_bundle, "refine", args.n_perm_refine)
            if perm_refine is not None:
                refine_source = "cache"
        if perm_refine is None:
            _log_progress(f"[6/8] Generating refine permutations (n={args.n_perm_refine}, jobs={jobs})")
            refine_seed = None if args.seed is None else args.seed + 7919
            refine_gen = PermutationGenerator(
                tree=tree,
                obs_mask_01=fg_01_mask,
                obs_mask_10=fg_10_mask,
                include_trait_loss=args.include_trait_loss,
                seed=refine_seed,
            )
            perm_refine = refine_gen.generate(args.n_perm_refine, jobs=jobs)
            refine_source = "generated"
            if cache_bundle is not None:
                put_stage_cache(cache_bundle, "refine", perm_refine)
        else:
            _log_progress(f"[6/8] Reusing refine permutations from cache (n={args.n_perm_refine})")

        context_refine = _build_worker_context(
            mode=args.mode,
            family_values=change.values,
            perm=perm_refine,
            direction_sign=direction_sign,
            fg_01_mask=fg_01_mask,
            fg_10_mask=fg_10_mask,
            all_mask=tree.all_mask,
            branch_lengths=branch_lengths if args.mode == "rate" else None,
            cafe_sig_masks=cafe_sig_masks,
        )
        refined_results = _evaluate_families(indices=refined_indices, jobs=jobs, context=context_refine)

        for i in refined_indices:
            rows[i].update(refined_results[i])
            rows[i]["n_perm_used"] = args.n_perm_refine
            rows[i]["refined"] = True
            pvalues[i] = rows[i]["p_empirical"]
        _log_progress("[6/8] Refinement completed")
    else:
        _log_progress("[6/8] Refinement skipped (no families passed refine criteria)")

    _log_progress("[7/8] Applying multiple-testing correction and writing result files")
    qvals = bh_adjust_with_none(pvalues)
    for row, q in zip(rows, qvals):
        row["q_bh"] = q

    out_prefix = Path(args.out_prefix)
    out_tsv = Path(str(out_prefix) + ".family_results.tsv")
    out_json = Path(str(out_prefix) + ".run_metadata.json")

    write_tsv(out_tsv, rows, family_fieldnames(args.mode))

    viz_outputs = generate_visual_outputs(
        rows=rows,
        out_prefix=out_prefix,
        qvalue_threshold=args.qvalue_threshold,
        pvalue_top_n=args.pvalue_top_n,
        hist_bins=args.hist_bins,
        make_plots=args.make_plots,
    )

    if args.perm_cache and cache_bundle is not None:
        save_cache_bundle(args.perm_cache, cache_bundle)
        _log_progress(f"[7/8] Updated permutation cache: {args.perm_cache}")

    metadata = {
        "tool": "permucn",
        "version": __version__,
        "inputs": {
            "cafe_dir": str(cafe_dir),
            "trait_tsv": str(trait_tsv),
            "change_table": str(paths["change"]),
            "branch_prob_table": str(paths["prob"]) if paths["prob"].exists() else None,
            "asr_tree": str(paths["asr_tree"]),
        },
        "parameters": {
            "mode": args.mode,
            "direction": args.direction,
            "include_trait_loss": args.include_trait_loss,
            "asr_method": args.asr_method,
            "asr_posterior_hi": args.asr_posterior_hi,
            "asr_posterior_lo": args.asr_posterior_lo,
            "cafe_significant_only": args.cafe_significant_only,
            "cafe_alpha": args.cafe_alpha,
            "n_perm_initial": args.n_perm_initial,
            "n_perm_refine": args.n_perm_refine,
            "refine_p_threshold": args.refine_p_threshold,
            "qvalue_threshold": args.qvalue_threshold,
            "pvalue_top_n": args.pvalue_top_n,
            "clade_bin_scheme": args.clade_bin_scheme,
            "seed": args.seed,
            "jobs_requested": args.jobs,
            "jobs_effective": jobs,
        },
        "trait_columns": {
            "species_column": trait.species_column,
            "trait_column_used": trait.trait_column,
            "trait_column_source": trait.trait_column_source,
            "row_count": trait.row_count,
        },
        "tree": {
            "n_nodes": len(tree.labels),
            "n_non_root_branches": len(tree.branch_key_by_index),
            "root_branch_key": root_key,
            "clade_bins": {
                "1": 0,
                "2": 1,
                "3-4": 2,
                "5-8": 3,
                "9-16": 4,
                "17-32": 5,
                "33-64": 6,
                "65+": 7,
            },
        },
        "asr": {
            "q01": asr.q01,
            "q10": asr.q10,
            "log_likelihood": asr.log_likelihood,
            "n_fg_01": n_fg_01,
            "n_fg_10": n_fg_10,
        },
        "permutation": {
            "cache": {
                "path": args.perm_cache,
                "cache_loaded": cache_loaded,
                "initial_source": initial_source,
                "refine_source": refine_source,
            },
            "initial": {
                "n_perm": args.n_perm_initial if perm_initial is not None else 0,
                "total_attempts": perm_initial.total_attempts if perm_initial is not None else 0,
                "total_restarts": perm_initial.total_restarts if perm_initial is not None else 0,
            },
            "refine": {
                "n_refined_families": len(refined_indices),
                "n_perm": args.n_perm_refine if refined_indices else 0,
                "total_attempts": perm_refine.total_attempts if perm_refine else 0,
                "total_restarts": perm_refine.total_restarts if perm_refine else 0,
            },
        },
        "results": {
            "n_families": len(rows),
            "n_tested": sum(1 for r in rows if r["status"] == "ok"),
            "n_refined": len(refined_indices),
            "output_tsv": str(out_tsv),
            "output_metadata_json": str(out_json),
            "visual_outputs": viz_outputs,
        },
    }
    write_json(out_json, metadata)
    _log_progress("[8/8] Run complete; outputs were written successfully")

    print(f"Wrote family results: {out_tsv}")
    print(f"Wrote metadata: {out_json}")
    print(f"Families analyzed: {len(rows)}")
    print(f"Families tested: {metadata['results']['n_tested']}")
    if refined_indices:
        print(f"Families refined: {len(refined_indices)}")

    return 0


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(args)
    except Exception as exc:  # pragma: no cover - top-level UX
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
