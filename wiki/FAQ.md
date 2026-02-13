# FAQ

## `Species mismatch between trait table and tree tips.`

Trait species and tree tip species must match exactly.

Check:

- spelling/case differences
- underscores vs spaces
- missing species rows
- extra species rows in trait table

## `Multiple binary trait columns detected (...)`

Auto-detection found multiple valid `0/1` columns.

Fix:

- specify `--trait-column <name>`

## `No binary trait column detected automatically.`

No candidate column satisfied strict binary value rules.

Fix:

- ensure intended trait column contains only `0`/`1`
- then pass `--trait-column`

## `--cafe-significant-only requires Gamma_branch_probabilities.tab`

Binary significance filtering needs branch probability input.

Fix:

- provide `Gamma_branch_probabilities.tab`, or
- remove `--cafe-significant-only`

## `Non-positive branch lengths found ... rate mode`

Rate mode divides branch deltas by branch lengths, so all non-root lengths must be positive.

Fix:

- use a tree with positive branch lengths, or
- run `--mode binary`

## All families show `status = no_valid_foreground`

No valid foreground branches were inferred from ASR under current thresholds.

Try:

- adjust `--asr-posterior-hi` / `--asr-posterior-lo`
- verify trait coding and species matching
- decide whether `--include-trait-loss` should be enabled

## Why are `*.pvalue_hist.tsv` and `*.qq.tsv` missing?

These files are generated only when at least one tested family has `p_empirical`.

Check:

- `run_metadata.json` -> `results.n_tested`
- `family_results.tsv` -> `status` column

## `top_hits.tsv` is empty

An empty `top_hits.tsv` usually means no families passed `q_bh <= --qvalue-threshold`.

Try:

- increase permutation counts for stability
- adjust biological filtering strategy
- use a more permissive `--qvalue-threshold` for exploratory review

## How should I choose permutation counts?

Practical pattern:

1. start with `--n-perm-initial 1000`
2. refine with larger `--n-perm-refine` (for example `100000` or `1000000`)
3. keep `--refine-p-threshold` small (for example `0.01`)

Use `--perm-cache` to avoid recomputing unchanged permutation sets.

## Is parallel mode reproducible?

Yes, with fixed `--seed`, outputs are deterministic across sequential and parallel runs (same settings and inputs).

## PDF plots were not created

`--make-plots` requires `matplotlib`.

If unavailable, run completes and records warning messages in:

- `run_metadata.json` -> `results.visual_outputs.plot_warnings`
