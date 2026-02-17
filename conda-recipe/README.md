# Bioconda Recipe Notes

This directory contains the recipe template for `permucn`.

## Install command

```bash
conda install -c conda-forge -c bioconda permucn
```

## Release flow (automated after initial registration)

1. Merge releasable commits into `main` (for `release-please`, use commit types like `feat:` / `fix:`).
2. `.github/workflows/release.yml` runs and opens/updates a release PR.
3. Merge the release PR. This automatically:
   - creates a Git tag and GitHub Release
   - publishes to PyPI
4. Bioconda autobump detects the new upstream release and opens/updates a PR in `bioconda/bioconda-recipes`.
5. Review and merge the Bioconda PR.

## One-time setup (manual)

1. Open the first PR to `bioconda/bioconda-recipes` with `recipes/permucn/meta.yaml`.
2. Merge that PR so `permucn` is registered in Bioconda recipes.
3. After that, regular version bumps are expected to be handled by Bioconda autobump.

## Local recipe check

Run a local recipe build check:

```bash
conda build conda-recipe
```

## Manual fallback (update recipe file)

If you need to update the recipe outside the release workflow:

```bash
python scripts/update_conda_recipe.py --version <VERSION>
```

## Manual fallback: submit to Bioconda

1. Fork `bioconda/bioconda-recipes`.
2. Copy this repository's latest `conda-recipe/meta.yaml` to `recipes/permucn/meta.yaml` in that fork.
3. Run Bioconda checks:

```bash
bioconda-utils lint recipes config.yml --packages permucn
bioconda-utils build recipes config.yml --packages permucn --docker
```

4. Open a PR to `bioconda/bioconda-recipes`.
