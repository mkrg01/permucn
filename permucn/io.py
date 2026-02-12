"""I/O utilities for permucn."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


_MISSING = {"", "NA", "N/A", "na", "n/a", "NaN", "nan"}
_SPECIES_CANDIDATES = [
    "species",
    "taxon",
    "taxon_id",
    "tip",
    "label",
    "name",
    "scientific_name",
]


@dataclass(frozen=True)
class TraitTable:
    species_to_state: Dict[str, int]
    species_column: str
    trait_column: str
    trait_column_source: str
    row_count: int


@dataclass(frozen=True)
class FamilyMatrix:
    family_ids: List[str]
    values: List[List[int]]


def _normalize_header(h: str) -> str:
    return h.strip().lower()


def _first_data_header(raw: str) -> str:
    if raw.startswith("#"):
        return raw[1:]
    return raw


def load_trait_table(path: str | Path, trait_column: str | None = None) -> TraitTable:
    """Load and validate species trait table with strict binary trait values."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Trait TSV not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Trait TSV has no header: {path}")

        headers = [h.strip() for h in reader.fieldnames]
        header_norm = {_normalize_header(h): h for h in headers}

        species_col = None
        for key in _SPECIES_CANDIDATES:
            if key in header_norm:
                species_col = header_norm[key]
                break
        if species_col is None:
            species_col = headers[0]

        rows = list(reader)

    if trait_column is not None:
        if trait_column not in headers:
            raise ValueError(
                f"Trait column '{trait_column}' not found in trait file headers: {headers}"
            )
        chosen_trait = trait_column
        source = "manual"
    else:
        candidate_cols: List[str] = []
        for col in headers:
            if col == species_col:
                continue
            ok = True
            for row in rows:
                value = (row.get(col) or "").strip()
                if value in _MISSING:
                    continue
                if value not in {"0", "1"}:
                    ok = False
                    break
            if ok:
                candidate_cols.append(col)

        if len(candidate_cols) == 1:
            chosen_trait = candidate_cols[0]
            source = "auto"
        elif len(candidate_cols) == 0:
            raise ValueError(
                "No binary trait column detected automatically. "
                "Please provide --trait-column explicitly."
            )
        else:
            raise ValueError(
                "Multiple binary trait columns detected "
                f"({candidate_cols}). Please provide --trait-column explicitly."
            )

    species_to_state: Dict[str, int] = {}
    for idx, row in enumerate(rows, start=2):
        species = (row.get(species_col) or "").strip()
        if not species:
            raise ValueError(f"Empty species value at {path}:{idx}")

        trait_raw = (row.get(chosen_trait) or "").strip()
        if trait_raw in _MISSING:
            raise ValueError(
                f"Missing trait value at {path}:{idx} for species '{species}' "
                f"in column '{chosen_trait}'"
            )
        if trait_raw not in {"0", "1"}:
            raise ValueError(
                f"Trait value must be 0/1 at {path}:{idx}; got '{trait_raw}'"
            )

        trait = int(trait_raw)
        if species in species_to_state and species_to_state[species] != trait:
            raise ValueError(
                f"Conflicting trait assignments for species '{species}' in {path}"
            )
        species_to_state[species] = trait

    return TraitTable(
        species_to_state=species_to_state,
        species_column=species_col,
        trait_column=chosen_trait,
        trait_column_source=source,
        row_count=len(rows),
    )


def read_cafe_header(path: str | Path) -> Tuple[str, List[str]]:
    """Read CAFE tab file header and return (family_col, branch_cols)."""
    path = Path(path)
    with path.open("r", encoding="utf-8") as handle:
        line = handle.readline().rstrip("\n")
    if not line:
        raise ValueError(f"Empty CAFE table: {path}")

    parts = line.split("\t")
    if not parts:
        raise ValueError(f"Invalid CAFE header: {path}")
    parts[0] = _first_data_header(parts[0])

    # Some CAFE tables include trailing tab/blank column.
    while parts and parts[-1] == "":
        parts.pop()

    family_col = parts[0]
    branch_cols = parts[1:]
    return family_col, branch_cols


def _safe_int(v: str) -> int:
    v = v.strip()
    if v in _MISSING:
        return 0
    if v.startswith("+"):
        v = v[1:]
    return int(v)


def _safe_float(v: str) -> float:
    v = v.strip()
    if v in _MISSING:
        return float("nan")
    return float(v)


def load_change_matrix(
    path: str | Path,
    branch_to_index: Dict[str, int],
    ignored_branch_keys: Iterable[str] | None = None,
) -> FamilyMatrix:
    """Load branch change matrix as dense per-family vectors in branch-index order."""
    path = Path(path)
    _, branch_cols = read_cafe_header(path)

    ignored = set(ignored_branch_keys or [])
    col_to_idx: List[int | None] = []
    unknown: List[str] = []
    for col in branch_cols:
        if col in branch_to_index:
            col_to_idx.append(branch_to_index[col])
        elif col in ignored or col == "":
            col_to_idx.append(None)
        else:
            col_to_idx.append(None)
            unknown.append(col)
    if unknown:
        preview = ", ".join(unknown[:6])
        raise ValueError(
            "Branch keys from change table not found in canonical tree: "
            f"{preview}{' ...' if len(unknown) > 6 else ''}"
        )

    family_ids: List[str] = []
    values: List[List[int]] = []
    width = len(branch_to_index)

    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline()
        if not header:
            raise ValueError(f"Empty change table: {path}")

        for line_no, line in enumerate(handle, start=2):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 1:
                continue

            fam = parts[0].strip()
            if not fam:
                raise ValueError(f"Missing family id at {path}:{line_no}")

            row = [0] * width
            for i, map_idx in enumerate(col_to_idx, start=1):
                if map_idx is None:
                    continue
                val = parts[i] if i < len(parts) else "0"
                row[map_idx] = _safe_int(val)

            family_ids.append(fam)
            values.append(row)

    return FamilyMatrix(family_ids=family_ids, values=values)


def load_probability_map(
    path: str | Path,
    branch_to_index: Dict[str, int],
    ignored_branch_keys: Iterable[str] | None = None,
) -> Dict[str, List[float]]:
    """Load branch probability table as a mapping family -> vector in branch-index order.

    Missing families are expected and handled by caller.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Branch probability table not found: {path}")

    _, branch_cols = read_cafe_header(path)

    ignored = set(ignored_branch_keys or [])
    col_to_idx: List[int | None] = []
    unknown: List[str] = []
    for col in branch_cols:
        if col in branch_to_index:
            col_to_idx.append(branch_to_index[col])
        elif col in ignored or col == "":
            col_to_idx.append(None)
        else:
            col_to_idx.append(None)
            unknown.append(col)
    if unknown:
        preview = ", ".join(unknown[:6])
        raise ValueError(
            "Branch keys from probability table not found in canonical tree: "
            f"{preview}{' ...' if len(unknown) > 6 else ''}"
        )

    out: Dict[str, List[float]] = {}
    width = len(branch_to_index)

    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline()
        if not header:
            raise ValueError(f"Empty probability table: {path}")

        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            fam = parts[0].strip()
            if not fam:
                continue

            vec = [float("nan")] * width
            for i, map_idx in enumerate(col_to_idx, start=1):
                if map_idx is None:
                    continue
                val = parts[i] if i < len(parts) else "N/A"
                vec[map_idx] = _safe_float(val)
            out[fam] = vec

    return out


def write_tsv(path: str | Path, rows: Sequence[Dict[str, object]], fieldnames: Sequence[str]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_json(path: str | Path, obj: object) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(obj, handle, ensure_ascii=True, indent=2, sort_keys=True)
        handle.write("\n")
