# Algorithm Notes

## Pipeline Summary

`permucn` runs the following high-level steps:

1. Parse canonical tree from `Gamma_asr.tre`.
2. Load and validate trait table.
3. Infer ancestral trait states by ML (`Mk2`) and define foreground branches.
4. Run binary statistics by either:
   - constrained permutations, or
   - Fisher exact test + Tarone screening.
5. Optionally refine low-p families with more permutations (permutation mode only).
6. Apply multiple-testing correction and write outputs.

## Tree Canonicalization

- Branch keys come from node labels after removing state suffixes (e.g., `_0`, `_1`).
- Non-root branches are indexed and represented as bitmasks.
- Per-branch metadata caches include:
  - ancestor masks
  - descendant masks
  - clade sizes
  - log2 clade-size bins

## Trait ASR (ML)

- Binary trait model with transition rates `q01` and `q10`
- Rates fitted by coarse-to-fine grid search in log-space
- Posterior per node computed by upward/downward dynamic programming
- Hard states assigned by thresholds:
  - state `1` if posterior >= `--asr-posterior-hi`
  - state `0` if posterior <= `--asr-posterior-lo`
  - otherwise ambiguous

Foreground transitions:

- `0->1` branches (`fg_01`)
- optional `1->0` branches (`fg_10`) unless `--no-include-trait-loss`

## Constrained Permutation

Permutation samples preserve structure using branch bitmasks:

- preserve observed clade-bin composition (`log2` scheme)
- avoid ancestor/descendant conflicts within the same sampled set
- when trait loss is included:
  - sample `1->0` from descendants of sampled `0->1` when feasible
  - fallback to independent sampling if descendant constraint is impossible

This keeps null permutations topology-aware while matching foreground size profile.

## Statistics

## Binary Mode

Per family, create branch sign masks from copy-number deltas:

- positive (`delta > 0`)
- negative (`delta < 0`)

Observed statistic is directional concordance count on foreground branches.

### Permutation path (`--binary-test permutation`)

Empirical p-value is one-sided:

`p = (k + 1) / (n + 1)`

where:

- `k` = number of permutation stats `>= observed`
- `n` = number of permutations

### Fisher + Tarone path (`--binary-test fisher-tarone`)

- Build a 2x2 table from:
  - foreground concordant / non-concordant
  - background concordant / non-concordant
- Compute one-sided Fisher exact p-value (foreground enrichment in selected direction).
- Compute each family's minimum attainable p-value from fixed margins.
- Apply Tarone screening to mark untestable families and derive a Tarone-Bonferroni threshold.
- In this path, `p_empirical` and `q_bh` are left empty.

## Rate Mode

Per branch rate:

- `rate_i = delta_i / branch_length_i`

Observed statistic is directional mean contrast over foreground masks (`fg_01`, `fg_10`).
Empirical p-value uses the same one-sided formula as binary mode.

## Two-Stage Permutation

Stage 1:

- compute all families with `--n-perm-initial`

Stage 2 (refinement):

- families with `p_empirical <= --refine-p-threshold`
- recompute with `--n-perm-refine` (if larger than initial)

This concentrates compute on candidate families.

## Multiple Testing

- Permutation path: BH correction is applied across tested families (`p_empirical` available), producing `q_bh`.
- Fisher + Tarone path: `q_bh` is not used; use `p_bonf_tarone` for adjusted significance.
- Families not tested keep `q_bh = None`.

## Parallelization and Determinism

- `--jobs 1`: sequential
- `--jobs 0`: auto CPU count
- `--jobs >= 2`: process pool (thread fallback on restricted systems)

Determinism:

- with fixed `--seed`, permutation masks are reproducible
- parallel and sequential runs produce identical statistical outputs for same seed and settings

## Cache Compatibility

Permutation cache reuse requires matching:

- tree fingerprint (branch order/identity)
- `include_trait_loss`
- foreground masks (`fg_01`, `fg_10`)

If incompatible, cache is ignored and regenerated.
