# Getting Started

## Prerequisites

- Python `>=3.12`
- CAFE output directory with required files
- Trait table (`.tsv`) with species and binary trait values (`0`/`1`)

## Install

```bash
pip install permucn
permucn --help
```

Conda / Bioconda:

```bash
conda install -c conda-forge -c bioconda permucn
permucn --help
```

Optional plotting support:

```bash
pip install "permucn[plots]"
```

Local development install:

```bash
pip install -e .
permucn --help
```

## Fetch Sample Data

```bash
permucn get-test-data --out-dir permucn_test_data
```

## Prepare Inputs

`--cafe-dir` must include:

- `Gamma_change.tab`
- `Gamma_asr.tre`

Optional/conditional files:

- `Gamma_branch_probabilities.tab` (required only when `--cafe-significant-only` is used)
- `Gamma_family_results.txt` (optional companion file)

Trait table requirements:

- TSV with header
- one species column
- one binary trait column
- species names must match tree tips exactly

See [Input Format](Input-Format.md) for strict rules and examples.

## First Run (Binary Mode)

```bash
permucn \
  --cafe-dir permucn_test_data/toy_example/cafe_output \
  --trait-tsv permucn_test_data/toy_example/species_trait.tsv \
  --no-include-trait-loss \
  --n-perm-initial 20 \
  --n-perm-refine 50 \
  --seed 7 \
  --out-prefix results/toy_binary
```

## First Run (Rate Mode)

```bash
permucn \
  --cafe-dir permucn_test_data/toy_example/cafe_output \
  --trait-tsv permucn_test_data/toy_example/species_trait.tsv \
  --mode rate \
  --seed 7 \
  --out-prefix results/toy_rate
```

## Larger Sample Dataset

```bash
permucn get-test-data --dataset polar_fish --out-dir permucn_test_data
```

```bash
permucn \
  --cafe-dir permucn_test_data/polar_fish/cafe_output \
  --trait-tsv permucn_test_data/polar_fish/species_trait.tsv \
  --jobs 4 \
  --perm-cache results/perm_cache.json.gz \
  --out-prefix results/polar_fish
```

## Check Outputs

Primary outputs:

- `results/*.family_results.tsv`
- `results/*.run_metadata.json`
- `results/*.top_hits.tsv`
- `results/*.top_pvalues.tsv` (default: `--pvalue-top-n 100`; set `--pvalue-top-n 0` to disable)

Diagnostics (at least one tested family has p-values):

- `results/*.pvalue_hist.tsv`
- `results/*.qq.tsv`

See [Output Interpretation](Output-Interpretation.md) for details.

## Recommended Next Steps

1. Increase permutation counts (`--n-perm-initial`, `--n-perm-refine`).
2. Enable cache reuse with `--perm-cache`.
3. Use `--cafe-significant-only` if you want branch-level CAFE significance filtering in binary mode.
4. Adjust ranked output size with `--pvalue-top-n`.
5. Generate PDFs with `--make-plots` when `matplotlib` is available.
