# Bioconda Recipe Notes

This directory contains the recipe template for `permucn`.

## Install command

```bash
conda install -c conda-forge -c bioconda permucn
```

## Release flow (automated)

1. Merge releasable commits into `main` (for `release-please`, use commit types like `feat:` / `fix:`).
2. `.github/workflows/release.yml` runs and opens/updates a release PR.
3. Merge the release PR. This automatically:
   - creates a Git tag and GitHub Release
   - publishes to PyPI
   - opens a PR in this repository to update `conda-recipe/meta.yaml` (`version` + `sha256`)
4. Merge the generated recipe PR in this repository.

## Local recipe check

Run a local recipe build check:

```bash
conda build conda-recipe
```

## Manual fallback (recipe update only)

If you need to update the recipe outside the release workflow:

```bash
python scripts/update_conda_recipe.py --version <VERSION>
```

## Submit to Bioconda

1. Fork `bioconda/bioconda-recipes`.
2. Copy this repository's latest `conda-recipe/meta.yaml` to `recipes/permucn/meta.yaml` in that fork.
3. Run Bioconda checks:

```bash
bioconda-utils lint recipes config.yml --packages permucn
bioconda-utils build recipes config.yml --packages permucn --docker
```

4. Open a PR to `bioconda/bioconda-recipes`.
