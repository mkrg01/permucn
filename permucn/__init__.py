"""permucn package."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
import tomllib

__all__ = ["__version__"]


def _read_local_version() -> str:
    """Read version from local pyproject.toml when package metadata is unavailable."""
    pyproject = Path(__file__).resolve().parent.parent / "pyproject.toml"
    data = tomllib.loads(pyproject.read_text(encoding="utf-8"))
    return str(data["project"]["version"])


try:
    __version__ = version("permucn")
except PackageNotFoundError:
    try:
        __version__ = _read_local_version()
    except Exception:
        __version__ = "0+unknown"
