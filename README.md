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
  --mode binary \
  --cafe-dir gfe_data/genome_evolution/cafe/cafe_output \
  --trait-tsv gfe_data/species_trait/species_trait.tsv \
  --jobs 4 \
  --n-perm-initial 1000 \
  --n-perm-refine 1000000 \
  --refine-p-threshold 0.01 \
  --perm-cache results/perm_cache.json.gz \
  --out-prefix results/permucn
```

## Outputs

- `<out-prefix>.family_results.tsv`
- `<out-prefix>.run_metadata.json`
- `<out-prefix>.top_hits.tsv`
- `<out-prefix>.pvalue_hist.tsv`
- `<out-prefix>.qq.tsv`

Optional plot outputs (when `--make-plots` is enabled and matplotlib is installed):

- `<out-prefix>.pvalue_hist.png`
- `<out-prefix>.qq.png`

`--jobs` parallelizes both permutation generation and per-family scoring.

Use `permucn --help` for all options.

## Run tests

```bash
python -m unittest discover -s tests -v
```
