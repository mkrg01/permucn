# Getting Started

## Prerequisites

- Python `>=3.12`
- CAFE output directory with required files
- Trait table (`.tsv`) with species and binary trait values (`0`/`1`)

## Install

```bash
pip install "permucn @ git+https://github.com/mkrg01/permucn.git"
permucn --help
```

After publishing to PyPI:

```bash
pip install permucn
```

Optional plotting support:

```bash
pip install "permucn[plots] @ git+https://github.com/mkrg01/permucn.git"
```

Local development install:

```bash
pip install -e .
permucn --help
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

See [Input Format](Input-Format) for strict rules and examples.

## First Run (Binary Mode)

```bash
permucn \
  --cafe-dir gfe_data/genome_evolution/cafe/cafe_output \
  --trait-tsv gfe_data/species_trait/species_trait.tsv \
  --mode binary \
  --direction gain \
  --seed 42 \
  --jobs 4 \
  --out-prefix results/first_binary
```

## First Run (Rate Mode)

```bash
permucn \
  --cafe-dir gfe_data/genome_evolution/cafe/cafe_output \
  --trait-tsv gfe_data/species_trait/species_trait.tsv \
  --mode rate \
  --direction gain \
  --seed 42 \
  --jobs 4 \
  --out-prefix results/first_rate
```

## Check Outputs

Primary outputs:

- `results/first_*.family_results.tsv`
- `results/first_*.run_metadata.json`
- `results/first_*.top_hits.tsv`

Diagnostic outputs:

- `results/first_*.pvalue_hist.tsv`
- `results/first_*.qq.tsv`

See [Output Interpretation](Output-Interpretation) for details.

## Recommended Next Steps

1. Increase permutation counts (`--n-perm-initial`, `--n-perm-refine`).
2. Enable cache reuse with `--perm-cache`.
3. Use `--cafe-significant-only` if you want branch-level CAFE significance filtering in binary mode.
4. Generate PDFs with `--make-plots` when `matplotlib` is available.
