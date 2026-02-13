"""Helpers for obtaining sample test datasets."""

from __future__ import annotations

import shutil
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from urllib.error import HTTPError, URLError


DATASET_FILES = {
    "toy_example": [
        "README.md",
        "species_trait.tsv",
        "cafe_output/Gamma_asr.tre",
        "cafe_output/Gamma_branch_probabilities.tab",
        "cafe_output/Gamma_change.tab",
        "cafe_output/Gamma_family_results.txt",
    ],
    "polar_fish": [
        "README.md",
        "species_trait.tsv",
        "cafe_output/Gamma_asr.tre",
        "cafe_output/Gamma_branch_probabilities.tab",
        "cafe_output/Gamma_category_likelihoods.txt",
        "cafe_output/Gamma_change.tab",
        "cafe_output/Gamma_clade_results.txt",
        "cafe_output/Gamma_count.tab",
        "cafe_output/Gamma_family_likelihoods.txt",
        "cafe_output/Gamma_family_results.txt",
        "cafe_output/Gamma_results.txt",
    ],
}
ROOT_FILES = ["README.md"]
DEFAULT_TEST_DATA_REPOSITORY = "mkrg01/permucn"
DEFAULT_TEST_DATA_REF = "main"


@dataclass(frozen=True)
class TestDataFetchResult:
    out_dir: Path
    datasets: tuple[str, ...]
    source: str


def _local_test_data_root() -> Path:
    # Editable/local usage keeps sample data at the project root.
    return Path(__file__).resolve().parents[1] / "test_data"


def _resolve_datasets(dataset: str) -> list[str]:
    if dataset == "all":
        return list(DATASET_FILES.keys())
    if dataset not in DATASET_FILES:
        valid = ", ".join(["all", *DATASET_FILES.keys()])
        raise ValueError(f"Unknown dataset '{dataset}'. Valid values: {valid}")
    return [dataset]


def _resolve_source(source: str, datasets: list[str]) -> str:
    if source == "github":
        return "github"

    local_root = _local_test_data_root()
    local_ok = local_root.exists() and all((local_root / name).is_dir() for name in datasets)

    if source == "local":
        if not local_ok:
            raise FileNotFoundError(
                "Local test_data directory is unavailable. Use --source github "
                "or run in the project repository."
            )
        return "local"

    if source == "auto":
        return "local" if local_ok else "github"

    raise ValueError(f"Unknown source '{source}'. Valid values: auto, local, github")


def _clear_existing_targets(out_dir: Path, datasets: list[str], force: bool) -> None:
    existing = [name for name in datasets if (out_dir / name).exists()]
    if existing and not force:
        preview = ", ".join(str(out_dir / name) for name in existing)
        raise FileExistsError(f"Output already exists: {preview}. Use --force to overwrite.")

    if force:
        for name in datasets:
            target = out_dir / name
            if target.is_dir():
                shutil.rmtree(target)
            elif target.exists():
                target.unlink()


def _copy_local(datasets: list[str], out_dir: Path, force: bool) -> None:
    source_root = _local_test_data_root()
    _clear_existing_targets(out_dir, datasets, force)

    for dataset in datasets:
        src = source_root / dataset
        if not src.is_dir():
            raise FileNotFoundError(f"Missing local dataset directory: {src}")
        shutil.copytree(src, out_dir / dataset)

    for name in ROOT_FILES:
        src = source_root / name
        dst = out_dir / name
        if src.exists() and (force or not dst.exists()):
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)


def _raw_file_url(repository: str, ref: str, relative_path: str) -> str:
    return f"https://raw.githubusercontent.com/{repository}/{ref}/test_data/{relative_path}"


def _download_file(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        with urllib.request.urlopen(url, timeout=120) as response, destination.open("wb") as handle:
            shutil.copyfileobj(response, handle)
    except HTTPError as exc:
        raise RuntimeError(f"Failed to download {url} (HTTP {exc.code})") from exc
    except URLError as exc:
        raise RuntimeError(f"Failed to download {url}: {exc.reason}") from exc


def _download_from_github(
    datasets: list[str],
    out_dir: Path,
    repository: str,
    ref: str,
    force: bool,
) -> None:
    _clear_existing_targets(out_dir, datasets, force)

    for dataset in datasets:
        for rel in DATASET_FILES[dataset]:
            url = _raw_file_url(repository, ref, f"{dataset}/{rel}")
            _download_file(url, out_dir / dataset / rel)

    for name in ROOT_FILES:
        url = _raw_file_url(repository, ref, name)
        dst = out_dir / name
        if force or not dst.exists():
            _download_file(url, dst)


def fetch_test_data(
    dataset: str,
    out_dir: str | Path,
    *,
    source: str = "auto",
    ref: str = DEFAULT_TEST_DATA_REF,
    repository: str = DEFAULT_TEST_DATA_REPOSITORY,
    force: bool = False,
) -> TestDataFetchResult:
    """Fetch sample datasets either from local checkout or GitHub."""
    datasets = _resolve_datasets(dataset)
    chosen_source = _resolve_source(source, datasets)

    target_dir = Path(out_dir)
    target_dir.mkdir(parents=True, exist_ok=True)

    if chosen_source == "local":
        _copy_local(datasets, target_dir, force)
    else:
        _download_from_github(datasets, target_dir, repository=repository, ref=ref, force=force)

    return TestDataFetchResult(
        out_dir=target_dir.resolve(),
        datasets=tuple(datasets),
        source=chosen_source,
    )
