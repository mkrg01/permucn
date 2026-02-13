# permucn Wiki

This Wiki contains detailed operational and technical documentation for `permucn`.

For a quick 5-minute introduction, start with `README.md` in the repository root.

## Start Here

1. [Getting Started](Getting-Started)
2. [Input Format](Input-Format)
3. [CLI Reference](CLI-Reference)
4. [Output Interpretation](Output-Interpretation)
5. [Algorithm Notes](Algorithm-Notes)
6. [FAQ](FAQ)

## Typical Workflow

1. Prepare CAFE output files and a species trait TSV.
2. Run `permucn` in `binary` or `rate` mode.
3. Inspect `*.family_results.tsv` and `*.run_metadata.json`.
4. Review `*.top_hits.tsv` and diagnostics (`*.pvalue_hist.tsv`, `*.qq.tsv`).
5. Tune permutation counts and rerun with `--perm-cache`.

## Minimal Example

```bash
permucn \
  --cafe-dir gfe_data/genome_evolution/cafe/cafe_output \
  --trait-tsv gfe_data/species_trait/species_trait.tsv \
  --mode binary \
  --out-prefix results/permucn
```

## Documentation Ownership

- Keep repository `README.md` short and onboarding-focused.
- Put long-form details and edge-case behavior in this Wiki.
- When CLI options change, update both README examples and Wiki references.
