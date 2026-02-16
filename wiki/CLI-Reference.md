# CLI Reference

## Usage

```bash
permucn [OPTIONS]
```

```bash
permucn get-test-data [OPTIONS]
```

## Test Data Command

| Option | Default | Description |
| --- | --- | --- |
| `--dataset {toy_example,polar_fish,all}` | `toy_example` | Which sample dataset to fetch. |
| `--out-dir` | `test_data` | Destination directory for fetched data. |
| `--source {auto,local,github}` | `auto` | `local`: copy from repository root, `github`: download from GitHub, `auto`: prefer local then github. |
| `--ref` | `main` | Git ref used when source is `github`. |
| `--force` | `False` | Overwrite existing dataset directory in `--out-dir`. |

## Core Inputs

| Option | Default | Description |
| --- | --- | --- |
| `--cafe-dir` | `None` | Directory containing CAFE output files. Required. |
| `--trait-tsv` | `None` | Trait TSV file path. Required. |
| `--trait-column` | `None` | Trait column name. If omitted, auto-detection is used. |
| `--out-prefix` | `permucn_results` | Output file prefix. |

## Analysis Mode

| Option | Default | Description |
| --- | --- | --- |
| `--mode {binary,rate}` | `binary` | Statistic mode. |
| `--direction {gain,loss}` | `gain` | Directional sign convention for statistics. |
| `--binary-test {permutation,fisher-tarone}` | `permutation` | Binary mode test engine. `fisher-tarone` runs one-sided Fisher exact tests and Tarone screening. |
| `--fwer-alpha` | `0.05` | Family-wise error rate used by Tarone-Bonferroni in `fisher-tarone` mode. |
| `--include-trait-loss` / `--no-include-trait-loss` | `include` | Whether to include inferred `1->0` foreground transitions. |
| `--asr-method {ml}` | `ml` | Trait ASR method. Current implementation supports ML only. |
| `--asr-posterior-hi` | `0.6` | Posterior threshold for hard state `1`. |
| `--asr-posterior-lo` | `0.4` | Posterior threshold for hard state `0`. |
| `--cafe-significant-only` | `False` | Binary mode only. Restrict event counting to branches significant in CAFE probabilities. |
| `--cafe-alpha` | `0.05` | Significance threshold used with `--cafe-significant-only`. |

## Permutation and Compute

| Option | Default | Description |
| --- | --- | --- |
| `--n-perm-initial` | `1000` | Number of permutations for initial stage. |
| `--n-perm-refine` | `1000000` | Number of permutations for refinement stage. |
| `--refine-p-threshold` | `0.01` | Families with initial `p <= threshold` are recomputed in refinement stage. |
| `--clade-bin-scheme {log2}` | `log2` | Clade-size binning scheme for constrained sampling. |
| `--seed` | `None` | Random seed for reproducibility. |
| `--jobs` | `1` | Worker count (`1` sequential, `0` auto CPU count). |
| `--perm-cache` | `None` | Optional permutation cache path (`.json` / `.json.gz`). |

Note: when `--binary-test fisher-tarone` is selected, permutation generation/refinement and permutation cache are not used.

## Output and Visualization

| Option | Default | Description |
| --- | --- | --- |
| `--qvalue-threshold` | `0.05` | Threshold for including families in `*.top_hits.tsv`. |
| `--hist-bins` | `20` | Number of bins for p-value histogram TSV/PDF. |
| `--make-plots` | `False` | Generate histogram and QQ PDFs if `matplotlib` is available. |

## Validation Rules

- `0 <= --asr-posterior-lo < --asr-posterior-hi <= 1`
- `--n-perm-initial > 0`
- `--n-perm-refine > 0`
- `0 < --refine-p-threshold < 1`
- `0 < --cafe-alpha < 1`
- `0 <= --qvalue-threshold <= 1`
- `--hist-bins > 0`
- `--jobs >= 0`
- `--cafe-significant-only` is allowed only in `binary` mode
- `0 < --fwer-alpha < 1`
- `--binary-test fisher-tarone` is allowed only in `binary` mode

## Common Commands

Binary run:

```bash
permucn \
  --cafe-dir <cafe_output_dir> \
  --trait-tsv <trait.tsv> \
  --mode binary \
  --out-prefix results/binary
```

Rate run:

```bash
permucn \
  --cafe-dir <cafe_output_dir> \
  --trait-tsv <trait.tsv> \
  --mode rate \
  --out-prefix results/rate
```

Reproducible run with cache:

```bash
permucn \
  --cafe-dir <cafe_output_dir> \
  --trait-tsv <trait.tsv> \
  --mode binary \
  --seed 123 \
  --perm-cache results/perm_cache.json.gz \
  --out-prefix results/repro
```
