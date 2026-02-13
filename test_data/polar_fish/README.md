# polar_fish

Sample dataset for running `permucn` under conditions closer to real data.

Contents:

- `test_data/polar_fish/cafe_output/`: CAFE outputs
- `test_data/polar_fish/species_trait.tsv`: species-level binary trait table

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
  --pvalue-top-n 100 \
  --perm-cache results/perm_cache.json.gz \
  --out-prefix results/polar_fish
```

This dataset takes longer to run than `toy_example`. If you only need a quick CLI check, use `test_data/toy_example/`.
