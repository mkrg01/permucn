# permucn

`permucn` is a permutation-based CLI for testing association between binary trait evolution and gene-family copy-number evolution from CAFE outputs.

## Install

```bash
pip install -e .
```

## Python support

- Python 3.12, 3.13, 3.14

## Example

```bash
permucn \
  --cafe-dir gfe_data/genome_evolution/cafe/cafe_output \
  --trait-tsv gfe_data/species_trait/species_trait.tsv \
  --jobs 4 \
  --mode binary \
  --refine-p-threshold 0.01 \
  --qvalue-threshold 0.05 \
  --perm-cache results/perm_cache.json.gz \
  --out-prefix results/permucn
```

## Outputs

- `<out-prefix>.family_results.tsv`
- `<out-prefix>.run_metadata.json`
- `<out-prefix>.top_hits.tsv` (families with `q_bh <= --qvalue-threshold`)
- `<out-prefix>.pvalue_hist.tsv`
- `<out-prefix>.qq.tsv`

Optional plot outputs (when `--make-plots` is enabled and matplotlib is installed):

- `<out-prefix>.pvalue_hist.pdf`
- `<out-prefix>.qq.pdf`

`--jobs` parallelizes both permutation generation and per-family scoring.

Use `permucn --help` for all options.

## Run tests

```bash
python -m unittest discover -s tests -v
```
