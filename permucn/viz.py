"""Visualization and summary output helpers."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


def generate_visual_outputs(
    rows: Sequence[Dict[str, object]],
    out_prefix: str | Path,
    qvalue_threshold: float = 0.05,
    pvalue_top_n: int = 0,
    hist_bins: int = 20,
    make_plots: bool = False,
) -> Dict[str, object]:
    out_prefix = Path(out_prefix)

    valid = [r for r in rows if r.get("status") == "ok" and r.get("p_empirical") is not None]
    pvals = [float(r["p_empirical"]) for r in valid]

    outputs: Dict[str, object] = {
        "top_hits_tsv": None,
        "top_pvalues_tsv": None,
        "pvalue_hist_tsv": None,
        "qq_tsv": None,
        "pvalue_hist_pdf": None,
        "qq_pdf": None,
        "plot_warnings": [],
    }

    top_path = Path(str(out_prefix) + ".top_hits.tsv")
    _write_top_hits(valid, top_path, qvalue_threshold)
    outputs["top_hits_tsv"] = str(top_path)

    if int(pvalue_top_n) > 0:
        top_p_path = Path(str(out_prefix) + ".top_pvalues.tsv")
        _write_top_pvalues(valid, top_p_path, pvalue_top_n)
        outputs["top_pvalues_tsv"] = str(top_p_path)

    if not pvals:
        return outputs

    hist_rows = _histogram_rows(pvals, bins=hist_bins)
    hist_tsv = Path(str(out_prefix) + ".pvalue_hist.tsv")
    _write_rows(hist_tsv, ["bin_start", "bin_end", "count"], hist_rows)
    outputs["pvalue_hist_tsv"] = str(hist_tsv)

    qq_rows = _qq_rows(pvals)
    qq_tsv = Path(str(out_prefix) + ".qq.tsv")
    _write_rows(qq_tsv, ["rank", "observed_p", "expected_p", "minus_log10_observed", "minus_log10_expected"], qq_rows)
    outputs["qq_tsv"] = str(qq_tsv)

    if make_plots:
        warning = _maybe_plot_hist(pvals, Path(str(out_prefix) + ".pvalue_hist.pdf"), bins=hist_bins)
        if warning is None:
            outputs["pvalue_hist_pdf"] = str(Path(str(out_prefix) + ".pvalue_hist.pdf"))
        else:
            outputs["plot_warnings"].append(warning)

        warning = _maybe_plot_qq(pvals, Path(str(out_prefix) + ".qq.pdf"))
        if warning is None:
            outputs["qq_pdf"] = str(Path(str(out_prefix) + ".qq.pdf"))
        else:
            outputs["plot_warnings"].append(warning)

    return outputs


def _write_top_hits(rows: Sequence[Dict[str, object]], path: Path, qvalue_threshold: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    ranked = sorted(
        rows,
        key=lambda r: (
            _none_last(r.get("q_bh")),
            _none_last(r.get("p_empirical")),
            -float(r.get("stat_obs", 0.0) or 0.0),
        ),
    )

    threshold = float(qvalue_threshold)
    keep = [row for row in ranked if row.get("q_bh") is not None and float(row["q_bh"]) <= threshold]
    fields = ["rank", "family_id", "q_bh", "p_empirical", "stat_obs", "mode", "direction", "status"]

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for i, row in enumerate(keep, start=1):
            writer.writerow(
                {
                    "rank": i,
                    "family_id": row.get("family_id"),
                    "q_bh": row.get("q_bh"),
                    "p_empirical": row.get("p_empirical"),
                    "stat_obs": row.get("stat_obs"),
                    "mode": row.get("mode"),
                    "direction": row.get("direction"),
                    "status": row.get("status"),
                }
            )


def _histogram_rows(values: Sequence[float], bins: int) -> List[Dict[str, object]]:
    bins = max(1, int(bins))
    counts = [0] * bins
    for p in values:
        p = min(max(p, 0.0), 1.0)
        idx = min(int(p * bins), bins - 1)
        counts[idx] += 1

    out: List[Dict[str, object]] = []
    width = 1.0 / bins
    for i, c in enumerate(counts):
        start = i * width
        end = (i + 1) * width
        out.append({"bin_start": start, "bin_end": end, "count": c})
    return out


def _write_top_pvalues(rows: Sequence[Dict[str, object]], path: Path, top_n: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    ranked = sorted(
        rows,
        key=lambda r: (
            _none_last(r.get("p_empirical")),
            _none_last(r.get("q_bh")),
            -float(r.get("stat_obs", 0.0) or 0.0),
        ),
    )

    keep = ranked[: max(0, int(top_n))]
    fields = ["rank", "family_id", "p_empirical", "q_bh", "stat_obs", "mode", "direction", "status"]

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for i, row in enumerate(keep, start=1):
            writer.writerow(
                {
                    "rank": i,
                    "family_id": row.get("family_id"),
                    "p_empirical": row.get("p_empirical"),
                    "q_bh": row.get("q_bh"),
                    "stat_obs": row.get("stat_obs"),
                    "mode": row.get("mode"),
                    "direction": row.get("direction"),
                    "status": row.get("status"),
                }
            )


def _qq_rows(values: Sequence[float]) -> List[Dict[str, object]]:
    n = len(values)
    obs = sorted(min(max(p, 1e-300), 1.0) for p in values)

    out: List[Dict[str, object]] = []
    for i, p in enumerate(obs, start=1):
        exp_p = i / (n + 1.0)
        out.append(
            {
                "rank": i,
                "observed_p": p,
                "expected_p": exp_p,
                "minus_log10_observed": -math.log10(p),
                "minus_log10_expected": -math.log10(exp_p),
            }
        )
    return out


def _write_rows(path: Path, fieldnames: Sequence[str], rows: Sequence[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def _maybe_plot_hist(values: Sequence[float], out_pdf: Path, bins: int) -> str | None:
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception:
        return "matplotlib is not available; skipped histogram PDF"

    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(values, bins=max(1, int(bins)), range=(0, 1), color="#3A78C2", edgecolor="black")
    ax.set_xlabel("Empirical p-value")
    ax.set_ylabel("Count")
    ax.set_title("p-value histogram")
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)
    return None


def _maybe_plot_qq(values: Sequence[float], out_pdf: Path) -> str | None:
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception:
        return "matplotlib is not available; skipped QQ PDF"

    obs = sorted(min(max(p, 1e-300), 1.0) for p in values)
    n = len(obs)
    exp = [i / (n + 1.0) for i in range(1, n + 1)]

    x = [-math.log10(p) for p in exp]
    y = [-math.log10(p) for p in obs]

    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(x, y, s=10, color="#2F7D4A", alpha=0.8)
    lim = max(max(x), max(y)) if x and y else 1.0
    ax.plot([0, lim], [0, lim], color="black", linewidth=1)
    ax.set_xlabel("Expected -log10(p)")
    ax.set_ylabel("Observed -log10(p)")
    ax.set_title("QQ plot")
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)
    return None


def _none_last(v: object) -> float:
    if v is None:
        return float("inf")
    return float(v)
