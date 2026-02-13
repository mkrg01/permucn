# permucn

`permucn` is a permutation-based CLI for testing associations between binary trait transitions and gene-family copy-number evolution from CAFE outputs.

## What This Tool Does

- Reads CAFE outputs (`Gamma_change.tab`, `Gamma_asr.tre`, optional `Gamma_branch_probabilities.tab`)
- Infers trait ancestral states by ML (`Mk2`) from a species trait table
- Defines foreground branches (`0->1`, optional `1->0`)
- Computes per-family empirical one-sided p-values by constrained branch permutation
- Applies BH correction (`q_bh`) and writes summary/diagnostic outputs

## README and Wiki Split

This `README.md` is the "start in 5 minutes" document.

Use the repository-local docs under `wiki/` for deeper references and long-form docs:

- Wiki Home: `wiki/Home.md`
- Getting Started: `wiki/Getting-Started.md`
- Input Format: `wiki/Input-Format.md`
- CLI Reference: `wiki/CLI-Reference.md`
- Output Interpretation: `wiki/Output-Interpretation.md`
- Algorithm Notes: `wiki/Algorithm-Notes.md`
- FAQ / Troubleshooting: `wiki/FAQ.md`

Suggested maintenance policy:

- Keep `README.md` focused on overview, install, quick start, and minimal troubleshooting.
- Move long explanations (input edge cases, algorithm details, output interpretation examples) to Wiki pages.
- When adding/changing CLI options, update both:
  - README quick examples
  - Wiki `CLI-Reference` and relevant deep-dive pages

## Installation

Requirements:

- Python `>=3.12` (`pyproject.toml`)

Install from source:

```bash
pip install -e .
```

Check CLI:

```bash
permucn --help
```

## Quick Start

Run with bundled example data:

```bash
permucn \
  --cafe-dir gfe_data/genome_evolution/cafe/cafe_output \
  --trait-tsv gfe_data/species_trait/species_trait.tsv \
  --mode binary \
  --direction gain \
  --jobs 4 \
  --n-perm-initial 1000 \
  --n-perm-refine 1000000 \
  --refine-p-threshold 0.01 \
  --qvalue-threshold 0.05 \
  --perm-cache results/perm_cache.json.gz \
  --out-prefix results/permucn
```

## Required Inputs

`--cafe-dir` must contain:

- `Gamma_change.tab` (required)
- `Gamma_asr.tre` (required, first `TREE ...;` entry is used)
- `Gamma_branch_probabilities.tab` (required only with `--cafe-significant-only`)
- `Gamma_family_results.txt` (optional companion file; not required by current CLI execution path)

`--trait-tsv` must be a tab-separated table with:

- one species column (auto-detected from common names like `species`, `taxon`, `name`; fallback is first column)
- one binary trait column (`0/1`)

Trait-column behavior:

- If exactly one binary column is found, it is auto-selected.
- If multiple binary columns are found, set `--trait-column`.
- Missing/invalid trait values fail fast.

Species names in trait TSV must exactly match tree tip species names (derived from tip labels such as `Acanthochromis_polyacanthus<66>` -> `Acanthochromis_polyacanthus`).

## Common Usage Patterns

Binary mode (default):

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

- `<out-prefix>.family_results.tsv`
- `<out-prefix>.run_metadata.json`
- `<out-prefix>.top_hits.tsv` (`q_bh <= --qvalue-threshold`)

Written when at least one tested family has p-values:

- `<out-prefix>.pvalue_hist.tsv`
- `<out-prefix>.qq.tsv`

Optional PDFs (`--make-plots`, requires `matplotlib`):

- `<out-prefix>.pvalue_hist.pdf`
- `<out-prefix>.qq.pdf`

### `family_results.tsv` Columns

Base columns:

- `family_id`
- `mode`
- `direction`
- `include_trait_loss`
- `n_fg_01`
- `n_fg_10`
- `stat_obs`
- `p_empirical`
- `q_bh`
- `n_perm_used`
- `refined`
- `status`

Binary-mode extra columns:

- `fg_concordant_count`
- `fg_total`
- `bg_concordant_count`
- `bg_total`
- `fg_concordance_rate`
- `bg_concordance_rate`

Rate-mode extra columns:

- `fg_mean_signed_rate`
- `bg_mean_signed_rate`
- `fg_median_signed_rate`
- `bg_median_signed_rate`

Status values:

- `ok`: tested normally
- `no_valid_foreground`: no inferred foreground branches

## Reproducibility and Performance

- Set `--seed` for deterministic permutations and reproducible outputs.
- Set `--jobs`:
  - `1`: sequential
  - `0`: auto (CPU count)
  - `>=2`: parallel permutation generation and family scoring
- Use `--perm-cache` (`.json` or `.json.gz`) to reuse compatible permutation sets across runs.

Cache compatibility checks include tree fingerprint, foreground masks, and `include_trait_loss`.

## Troubleshooting

- `Species mismatch between trait table and tree tips.`  
  Align species names exactly between trait table and tip labels.
- `Multiple binary trait columns detected ...`  
  Specify `--trait-column`.
- `No binary trait column detected automatically ...`  
  Ensure target trait column is strictly `0/1`, then pass `--trait-column`.
- `--cafe-significant-only requires Gamma_branch_probabilities.tab`  
  Add that file or disable `--cafe-significant-only`.
- `Non-positive branch lengths found ... rate mode`  
  Use a tree with positive branch lengths or switch to `binary` mode.

## Development

Run tests:

```bash
python -m unittest discover -s tests -v
```

## License

MIT (see `LICENSE`).
