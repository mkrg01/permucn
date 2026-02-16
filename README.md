# permucn

[![PyPI](https://img.shields.io/pypi/v/permucn.svg)](https://pypi.org/project/permucn/) [![Python](https://img.shields.io/pypi/pyversions/permucn.svg)](https://pypi.org/project/permucn/) [![Tests](https://github.com/mkrg01/permucn/actions/workflows/tests.yml/badge.svg)](https://github.com/mkrg01/permucn/actions/workflows/tests.yml) [![Conda Recipe](https://github.com/mkrg01/permucn/actions/workflows/conda-recipe.yml/badge.svg)](https://github.com/mkrg01/permucn/actions/workflows/conda-recipe.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

`permucn` is a command-line tool for testing whether trait transitions (`0->1`, optionally `1->0`) are associated with gene-family copy-number evolution from CAFE outputs.
In binary mode, you can choose either permutation-based testing or Fisher's exact test with Tarone screening.

## Install

Requirements: Python `>=3.12`

```bash
pip install permucn
```

Conda / Bioconda:

```bash
conda install -c conda-forge -c bioconda permucn
```

Optional plot support:

```bash
pip install "permucn[plots]"
```

## Quick Start

1. Fetch sample data (works after `pip install`):

```bash
permucn get-test-data --out-dir permucn_test_data
```

2. Run `permucn`:

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

3. Check outputs:

- `results/toy_binary.family_results.tsv`
- `results/toy_binary.run_metadata.json`
- `results/toy_binary.top_hits.tsv`
- `results/toy_binary.top_pvalues.tsv` (default: top 100 by p-value)

If at least one tested family has p-values, these are also written:

- `results/toy_binary.pvalue_hist.tsv`
- `results/toy_binary.qq.tsv`

For a larger sample dataset:

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

## Required Inputs

`--cafe-dir` must include:

- `Gamma_change.tab` (required)
- `Gamma_asr.tre` (required)
- `Gamma_branch_probabilities.tab` (required only with `--cafe-significant-only`)
- `Gamma_family_results.txt` (optional)

`--trait-tsv` must be a TSV with:

- one species column (`species`, `taxon`, `name`, etc.; first column fallback)
- one binary trait column (`0/1`)

Important rules:

- Species names in trait TSV must match tree tip names.
- If multiple binary trait columns exist, specify `--trait-column`.
- `--cafe-significant-only` is valid only in `binary` mode.
- `rate` mode requires strictly positive non-root branch lengths.

## Common Commands

Default binary mode:

```bash
permucn --cafe-dir <cafe_output_dir> --trait-tsv <trait.tsv> --out-prefix results/binary
```

Binary mode with CAFE-significant events only:

```bash
permucn --cafe-dir <cafe_output_dir> --trait-tsv <trait.tsv> --cafe-significant-only --out-prefix results/binary_sig
```

Binary mode with Fisher + Tarone:

```bash
permucn \
  --cafe-dir <cafe_output_dir> \
  --trait-tsv <trait.tsv> \
  --binary-test fisher-tarone \
  --out-prefix results/binary_fisher
```

Rate mode:

```bash
permucn --cafe-dir <cafe_output_dir> --trait-tsv <trait.tsv> --mode rate --out-prefix results/rate
```

## Reproducibility and Performance

- `--seed`: reproducible permutations
- `--jobs`: parallelism (`1` sequential, `0` auto CPU)
- `--perm-cache`: reuse permutations across runs (`.json` / `.json.gz`)
- `--pvalue-top-n`: size of `*.top_pvalues.tsv` (`0` disables)

## Documentation

- [wiki/Getting-Started.md](wiki/Getting-Started.md)
- [wiki/CLI-Reference.md](wiki/CLI-Reference.md)
- [wiki/Input-Format.md](wiki/Input-Format.md)
- [wiki/Output-Interpretation.md](wiki/Output-Interpretation.md)
- [wiki/Algorithm-Notes.md](wiki/Algorithm-Notes.md)
- [wiki/FAQ.md](wiki/FAQ.md)

## License

MIT ([LICENSE](LICENSE))
