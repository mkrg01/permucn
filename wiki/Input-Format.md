# Input Format

## Overview

`permucn` consumes:

- a CAFE output directory (`--cafe-dir`)
- a trait table (`--trait-tsv`)

The tool validates both structure and value constraints strictly.

## CAFE Directory (`--cafe-dir`)

### Required Files

- `Gamma_change.tab`
- `Gamma_asr.tre`

### Conditional File

- `Gamma_branch_probabilities.tab` is required only when `--cafe-significant-only` is enabled.

### Optional Companion

- `Gamma_family_results.txt` may exist in standard CAFE outputs but is not required by the current core run path.

## `Gamma_asr.tre` Requirements

- NEXUS format with at least one `TREE ... = ...;` entry
- only the first tree entry is read
- branch lengths are read from the tree
- in `rate` mode, all non-root branch lengths must be `> 0`

Branch keys are derived from node labels by stripping state suffixes:

- `Acanthochromis_polyacanthus<66>_1` -> `Acanthochromis_polyacanthus<66>`
- `<123>_0` -> `<123>`

## `Gamma_change.tab` Requirements

- tab-separated text
- first column is family ID
- remaining columns are branch keys
- each row contains per-branch copy-number changes (e.g., `+1`, `0`, `-2`)

Notes:

- unknown branch keys (not present in canonical tree) cause an error
- missing numeric entries in branch columns are interpreted as `0`
- root branch is ignored

## `Gamma_branch_probabilities.tab` Requirements

- tab-separated text
- first column is family ID (often `#FamilyID`)
- remaining columns are branch keys consistent with the tree
- values are branch-level probabilities from CAFE

This table is used only when `--cafe-significant-only` is set. Branches are filtered by `--cafe-alpha`.

## Trait Table (`--trait-tsv`) Requirements

### Structure

- tab-separated text with header
- one species column
- one binary trait column (`0` or `1`)

### Species Column Detection

Auto-detection priority (case-insensitive):

- `species`
- `taxon`
- `taxon_id`
- `tip`
- `label`
- `name`
- `scientific_name`

If none match, the first column is used as species column.

### Trait Column Detection

- If `--trait-column` is provided, that column is used.
- Otherwise, columns except species column are scanned.
- A column is considered binary if all non-missing values are `0`/`1`.
- Exactly one binary candidate is required for auto-detection.

Failure cases:

- no binary candidate -> error
- multiple binary candidates -> error (`--trait-column` required)

### Missing Values

Recognized missing tokens include:

- empty string
- `NA`, `N/A`, `na`, `n/a`, `NaN`, `nan`

For the selected trait column, missing values are not allowed and produce an error.

### Species Matching Rule

Trait species set must match tree tip species set exactly:

- missing species in trait table -> error
- extra species not found in tree -> error

## Quick Validation Checklist

1. Branch key sets in CAFE tables and tree are consistent.
2. Trait table has exactly one intended binary trait column.
3. Species names match tree tips exactly.
4. For `rate` mode, all branch lengths are positive.
