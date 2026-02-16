# Output Interpretation

## Output Files

Given `--out-prefix results/run1`, `permucn` writes:

- `results/run1.family_results.tsv`
- `results/run1.run_metadata.json`
- `results/run1.top_hits.tsv`

If at least one tested family has p-values:

- `results/run1.pvalue_hist.tsv`
- `results/run1.qq.tsv`

If `--make-plots` and `matplotlib` is available:

- `results/run1.pvalue_hist.pdf`
- `results/run1.qq.pdf`

## `family_results.tsv`

One row per family.

### Core Columns

- `family_id`: family identifier from CAFE table
- `mode`: `binary` or `rate`
- `direction`: `gain` or `loss`
- `include_trait_loss`: whether `1->0` foreground was included
- `n_fg_01`, `n_fg_10`: number of inferred foreground branches
- `stat_obs`: observed test statistic
- `p_empirical`: primary p-value used for ranking
  - permutation mode: one-sided empirical p-value
  - fisher-tarone mode: one-sided Fisher exact p-value
- `q_bh`: BH-adjusted q-value (blank in fisher-tarone mode)
- `n_perm_used`: permutation count used for the final p-value
- `refined`: whether this family was recomputed in refinement stage
- `status`: `ok`, `untestable_tarone`, or `no_valid_foreground`

### Binary Mode Additions

- `fg_concordant_count`, `fg_total`
- `bg_concordant_count`, `bg_total`
- `fg_concordance_rate`, `bg_concordance_rate`
- `p_fisher`: one-sided Fisher exact p-value (set in fisher-tarone mode)
- `p_min_attainable`: minimum attainable Fisher p-value for fixed margins
- `tarone_testable`: whether the family is testable under Tarone threshold
- `p_bonf_tarone`: Tarone-Bonferroni adjusted p-value (testable families only)
- `reject_tarone`: whether Fisher p-value passes Tarone-Bonferroni threshold

Interpretation:

- Higher `stat_obs` means stronger concordance in the specified direction.
- Compare foreground vs background concordance rates for effect context.

### Rate Mode Additions

- `fg_mean_signed_rate`, `bg_mean_signed_rate`
- `fg_median_signed_rate`, `bg_median_signed_rate`

Interpretation:

- Signed rates align with selected `--direction`.
- Foreground vs background summaries help judge practical effect size.

## `top_hits.tsv`

Contains families with `q_bh <= --qvalue-threshold`.

Ranking order:

1. smaller `q_bh`
2. smaller `p_empirical`
3. larger `stat_obs`

Columns:

- `rank`
- `family_id`
- `q_bh`
- `p_empirical`
- `stat_obs`
- `mode`
- `direction`
- `status`

## `run_metadata.json`

Structured run metadata for auditability:

- `tool`, `version`
- `inputs`
- `parameters`
- `trait_columns`
- `tree`
- `asr`
- `permutation`
- `tarone`
- `results`

Useful fields for reproducibility:

- `parameters.seed`
- `parameters.jobs_requested`, `parameters.jobs_effective`
- `permutation.cache` section
- `permutation.initial` / `permutation.refine` attempt statistics
- `tarone.m_total`, `tarone.m_testable`, `tarone.bonferroni_denom`, `tarone.threshold`

## Diagnostic Tables

`pvalue_hist.tsv`:

- histogram bins (`bin_start`, `bin_end`, `count`) over empirical p-values

`qq.tsv`:

- expected vs observed p-values and `-log10` transforms for QQ diagnostics

## Interpreting Empty/Reduced Outputs

- If all families are `no_valid_foreground`, p-values are absent and histogram/QQ files are not written.
- In fisher-tarone mode, families marked `untestable_tarone` are excluded from `n_tested` and p-value diagnostics.
- `top_hits.tsv` is still written, but can contain only the header row.
- In metadata, check `results.n_tested` and `asr.n_fg_01` / `asr.n_fg_10`.
