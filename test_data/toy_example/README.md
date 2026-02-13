# toy_example

This is a minimal, reproducible input dataset for quickly validating `permucn`.

Contents:

- `test_data/toy_example/cafe_output/Gamma_asr.tre`
- `test_data/toy_example/cafe_output/Gamma_change.tab`
- `test_data/toy_example/cafe_output/Gamma_branch_probabilities.tab`
- `test_data/toy_example/cafe_output/Gamma_family_results.txt`
- `test_data/toy_example/species_trait.tsv`

Example run:

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
