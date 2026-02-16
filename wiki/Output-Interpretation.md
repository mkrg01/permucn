# Output Interpretation

## Output Files

Given `--out-prefix results/run1`, `permucn` writes:

- `results/run1.family_results.tsv`
- `results/run1.run_metadata.json`
- `results/run1.top_hits.tsv`
- `results/run1.top_pvalues.tsv` (default: `--pvalue-top-n 100`; set `--pvalue-top-n 0` to disable)

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
- `p_empirical`: one-sided empirical p-value (permutation mode; blank in fisher-tarone mode)
- `q_bh`: BH-adjusted q-value (permutation mode; blank in fisher-tarone mode)
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

Contains families passing adjusted threshold.

Ranking order:

1. smaller adjusted value
2. smaller p-value
3. larger `stat_obs`

Columns by mode:

- permutation mode: `rank`, `family_id`, `q_bh`, `p_empirical`, `stat_obs`, `mode`, `direction`, `status`
- fisher-tarone mode: `rank`, `family_id`, `p_bonf_tarone`, `p_fisher`, `stat_obs`, `mode`, `direction`, `status`

## `top_pvalues.tsv`

Contains the top N families ranked by primary p-value:

1. smaller primary p-value
2. smaller adjusted value
3. larger `stat_obs`

Primary/adjusted column pairs by mode:

- permutation mode: `p_empirical`, `q_bh`
- fisher-tarone mode: `p_fisher`, `p_bonf_tarone`

Use `--pvalue-top-n` to change N (`0` disables this file).

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

- histogram bins (`bin_start`, `bin_end`, `count`) over primary p-values
  - permutation mode: `p_empirical`
  - fisher-tarone mode: `p_fisher`

`qq.tsv`:

- expected vs observed p-values and `-log10` transforms for QQ diagnostics

## Interpreting Empty/Reduced Outputs

- If all families are `no_valid_foreground`, p-values are absent and histogram/QQ files are not written.
- In fisher-tarone mode, families marked `untestable_tarone` are excluded from `n_tested` and p-value diagnostics.
- `top_hits.tsv` and `top_pvalues.tsv` (if enabled) are still written, but can contain only the header row.
- In metadata, check `results.n_tested` and `asr.n_fg_01` / `asr.n_fg_10`.
