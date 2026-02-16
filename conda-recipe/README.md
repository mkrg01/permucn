# Bioconda Recipe Notes

This directory contains the recipe template for `permucn`.

## Install command

```bash
conda install -c conda-forge -c bioconda permucn
```

## Release update checklist

1. Bump version in:
   - `pyproject.toml`
   - `permucn/__init__.py`
   - `conda-recipe/meta.yaml`
2. Build source distribution and compute hash:

```bash
python -m build --sdist
sha256sum dist/permucn-<VERSION>.tar.gz
```

3. Update `sha256` in `conda-recipe/meta.yaml`.

## Local recipe check

`meta.yaml` supports local source mode for validation:

```bash
PERMUCN_CONDA_LOCAL_SOURCE=1 conda build conda-recipe
```

## Submit to Bioconda

1. Fork `bioconda/bioconda-recipes`.
2. Copy `conda-recipe/meta.yaml` to `recipes/permucn/meta.yaml` in that fork.
3. Run Bioconda checks:

```bash
bioconda-utils lint recipes config.yml --packages permucn
bioconda-utils build recipes config.yml --packages permucn --docker
```

4. Open a PR to `bioconda/bioconda-recipes`.
