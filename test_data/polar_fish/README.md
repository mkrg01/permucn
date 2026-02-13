# polar_fish

This dataset mirrors the current polar fish example inputs.

Contents:

- `test_data/polar_fish/cafe_output/`
- `test_data/polar_fish/species_trait.tsv`

Example run:

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
  --perm-cache results/perm_cache.json.gz \
  --out-prefix results/polar_fish
```
