#!/usr/bin/env python3
"""Update conda recipe version and sha256 from a PyPI release."""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

VERSION_PATTERN = re.compile(r'(\{%\s*set\s+version\s*=\s*")[^"]+("\s*%\})')
SHA256_PATTERN = re.compile(r"(^\s*sha256:\s*)([0-9a-fA-F]+)(\s*$)", flags=re.MULTILINE)


def fetch_sdist_sha256(package: str, release_version: str, retries: int, retry_wait: int) -> str:
    url = f"https://pypi.org/pypi/{package}/{release_version}/json"
    last_error: Exception | None = None

    for attempt in range(1, retries + 1):
        try:
            with urlopen(url, timeout=20) as response:
                payload = json.load(response)

            for file_info in payload.get("urls", []):
                if file_info.get("packagetype") == "sdist":
                    digest = file_info.get("digests", {}).get("sha256")
                    if digest:
                        return str(digest)
            raise RuntimeError("No sdist artifact was found in PyPI JSON response")
        except (HTTPError, URLError, TimeoutError, json.JSONDecodeError, RuntimeError) as exc:
            last_error = exc
            if attempt == retries:
                break
            print(
                f"[update_conda_recipe] attempt {attempt}/{retries} failed: {exc} "
                f"(retrying in {retry_wait}s)",
                file=sys.stderr,
            )
            time.sleep(retry_wait)

    raise RuntimeError(
        f"Failed to fetch sdist sha256 for {package}=={release_version}: {last_error}"
    ) from last_error


def update_recipe(recipe_path: Path, release_version: str, sha256: str) -> bool:
    old_text = recipe_path.read_text(encoding="utf-8")
    text = old_text

    text, version_replacements = VERSION_PATTERN.subn(
        lambda m: f"{m.group(1)}{release_version}{m.group(2)}", text, count=1
    )
    if version_replacements != 1:
        raise RuntimeError(
            f"Expected exactly one version assignment in {recipe_path}, found {version_replacements}"
        )

    text, sha_replacements = SHA256_PATTERN.subn(
        lambda m: f"{m.group(1)}{sha256}{m.group(3)}", text, count=1
    )
    if sha_replacements != 1:
        raise RuntimeError(
            f"Expected exactly one sha256 assignment in {recipe_path}, found {sha_replacements}"
        )

    changed = text != old_text
    if changed:
        recipe_path.write_text(text, encoding="utf-8")
    return changed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", required=True, help="Released package version")
    parser.add_argument("--package", default="permucn", help="PyPI package name")
    parser.add_argument(
        "--recipe",
        type=Path,
        default=Path("conda-recipe/meta.yaml"),
        help="Path to conda recipe meta.yaml",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=12,
        help="Maximum number of retries for PyPI JSON fetch",
    )
    parser.add_argument(
        "--retry-wait",
        type=int,
        default=10,
        help="Seconds to wait between retries",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.retries < 1:
        raise ValueError("--retries must be >= 1")
    if args.retry_wait < 0:
        raise ValueError("--retry-wait must be >= 0")
    if not args.recipe.exists():
        raise FileNotFoundError(f"Recipe file not found: {args.recipe}")

    sha256 = fetch_sdist_sha256(args.package, args.version, args.retries, args.retry_wait)
    changed = update_recipe(args.recipe, args.version, sha256)

    status = "updated" if changed else "already up-to-date"
    print(
        f"[update_conda_recipe] {status}: recipe={args.recipe} "
        f"version={args.version} sha256={sha256}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
