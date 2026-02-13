# permucn

`permucn` is a CLI tool for testing whether trait transitions (`0->1`, optionally `1->0`) are associated with gene-family copy-number evolution, using CAFE outputs.

## What You Can Do

- Run association tests per gene family with permutation-based empirical p-values
- Control FDR using BH-adjusted q-values (`q_bh`)
- Focus on trait-gain branches (`0->1`) with optional trait-loss inclusion
- Export reproducible result tables and diagnostics

## Install

Requirements:

- Python `>=3.12`

Install from source:

```bash
pip install -e .
```

Check CLI:

```bash
permucn --help
```

## Quick Start

Run with bundled toy data (fast check):

```bash
permucn \
  --cafe-dir test_data/toy_example/cafe_output \
  --trait-tsv test_data/toy_example/species_trait.tsv \
  --mode binary \
  --no-include-trait-loss \
  --asr-posterior-hi 0.6 \
  --asr-posterior-lo 0.4 \
  --n-perm-initial 20 \
  --n-perm-refine 50 \
  --seed 7 \
  --out-prefix results/toy_binary
```

Run with the larger polar fish dataset:

```bash
permucn \
  --cafe-dir test_data/polar_fish/cafe_output \
  --trait-tsv test_data/polar_fish/species_trait.tsv \
  --mode binary \
  --direction gain \
  --jobs 4 \
  --n-perm-initial 1000 \
  --n-perm-refine 1000000 \
  --refine-p-threshold 0.01 \
  --qvalue-threshold 0.05 \
  --pvalue-top-n 100 \
  --perm-cache results/perm_cache.json.gz \
  --out-prefix results/polar_fish
```

## Required Inputs

`--cafe-dir` must contain:

- `Gamma_change.tab` (required)
- `Gamma_asr.tre` (required; first `TREE ...;` entry is used)
- `Gamma_branch_probabilities.tab` (required only with `--cafe-significant-only`)
- `Gamma_family_results.txt` (optional)

Tree branch length rules:

- ASR requires finite branch lengths (`NaN`/`inf` not allowed).
- Negative branch lengths are not allowed.
- `rate` mode additionally requires strictly positive non-root branch lengths.

`--trait-tsv` must be a tab-separated table with:

- one species column (auto-detected from common names like `species`, `taxon`, `name`; fallback is first column)
- one binary trait column (`0/1`)

Trait column rules:

- If exactly one binary column exists, it is auto-selected.
- If multiple binary columns exist, specify `--trait-column`.
- Missing or invalid trait values stop the run with an error.

Species names in trait TSV must exactly match tree tip species names (for example, `Acanthochromis_polyacanthus<66>` maps to `Acanthochromis_polyacanthus`).

## Common Usage

Default binary mode:

```bash
permucn \
  --cafe-dir <cafe_output_dir> \
  --trait-tsv <trait.tsv> \
  --mode binary \
  --out-prefix results/binary_run
```

Binary mode using only CAFE-significant branch events:

```bash
permucn \
  --cafe-dir <cafe_output_dir> \
  --trait-tsv <trait.tsv> \
  --mode binary \
  --cafe-significant-only \
  --cafe-alpha 0.05 \
  --out-prefix results/binary_sig
```

Rate mode:

```bash
permucn \
  --cafe-dir <cafe_output_dir> \
  --trait-tsv <trait.tsv> \
  --mode rate \
  --direction gain \
  --out-prefix results/rate_run
```

Notes:

- `--cafe-significant-only` is valid only in `binary` mode.
- `rate` mode requires strictly positive branch lengths in the canonical tree.

## Outputs

Always written:

- `<out-prefix>.family_results.tsv` (main per-family results)
- `<out-prefix>.run_metadata.json` (run settings and metadata)
- `<out-prefix>.top_hits.tsv` (families passing `q_bh <= --qvalue-threshold`)
- `<out-prefix>.top_pvalues.tsv` (top `--pvalue-top-n` families by smallest `p_empirical`, only when `--pvalue-top-n > 0`)

Written when at least one tested family has p-values:

- `<out-prefix>.pvalue_hist.tsv`
- `<out-prefix>.qq.tsv`

Optional PDFs (`--make-plots`, requires `matplotlib`):

- `<out-prefix>.pvalue_hist.pdf`
- `<out-prefix>.qq.pdf`

### `family_results.tsv` Key Columns

- `family_id`: gene family identifier
- `p_empirical`: empirical one-sided p-value from permutations
- `q_bh`: BH-adjusted p-value
- `stat_obs`: observed test statistic
- `n_perm_used`: number of permutations used for that family
- `status`: test status (`ok`, `no_valid_foreground`)

Binary-mode extra columns:

- `fg_concordance_rate`
- `bg_concordance_rate`

Rate-mode extra columns:

- `fg_mean_signed_rate`
- `bg_mean_signed_rate`

## Reproducibility and Performance

- Use `--seed` to make results reproducible.
- Use `--jobs`:
  - `1`: sequential
  - `0`: auto-detect CPU count
  - `>=2`: parallel processing
- Use `--perm-cache` (`.json` or `.json.gz`) to reuse compatible permutation sets across runs.

## More Documentation

For detailed references:

- `wiki/Home.md`
- `wiki/Getting-Started.md`
- `wiki/Input-Format.md`
- `wiki/CLI-Reference.md`
- `wiki/Output-Interpretation.md`
- `wiki/Algorithm-Notes.md`
- `wiki/FAQ.md`

## Troubleshooting

- `Species mismatch between trait table and tree tips.`  
  Make species names identical between your trait table and tree tips.
- `Multiple binary trait columns detected ...`  
  Specify `--trait-column`.
- `No binary trait column detected automatically ...`  
  Ensure target trait column is strictly `0/1`, then pass `--trait-column`.
- `--cafe-significant-only requires Gamma_branch_probabilities.tab`  
  Add that file or disable `--cafe-significant-only`.
- `Non-positive branch lengths found ... rate mode`  
  Use a tree with positive branch lengths or switch to `binary` mode.
- `Invalid branch lengths in tree for ASR ...`  
  Ensure all non-root branches in `Gamma_asr.tre` have finite, non-negative lengths.

## License

MIT (see `LICENSE`).
