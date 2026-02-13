import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from permucn.testdata import fetch_test_data


class TestTestDataFetch(unittest.TestCase):
    def test_auto_uses_github_when_local_data_is_missing(self) -> None:
        td = Path(tempfile.mkdtemp(prefix="permucn_testdata_github_"))
        out_dir = td / "out"
        fake_local_root = td / "missing_local_test_data"
        downloaded: list[tuple[str, Path]] = []

        def fake_download(url: str, destination: Path) -> None:
            downloaded.append((url, destination))
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text("stub\n", encoding="utf-8")

        with patch("permucn.testdata._local_test_data_root", return_value=fake_local_root):
            with patch("permucn.testdata._download_file", side_effect=fake_download):
                result = fetch_test_data(
                    dataset="toy_example",
                    out_dir=out_dir,
                    source="auto",
                    repository="org/repo",
                    ref="feature",
                )

        self.assertEqual(result.source, "github")
        self.assertTrue((out_dir / "toy_example" / "species_trait.tsv").exists())
        urls = [u for u, _ in downloaded]
        self.assertIn(
            "https://raw.githubusercontent.com/org/repo/feature/test_data/toy_example/species_trait.tsv",
            urls,
        )

    def test_local_source_requires_local_data(self) -> None:
        td = Path(tempfile.mkdtemp(prefix="permucn_testdata_local_"))
        out_dir = td / "out"
        fake_local_root = td / "missing_local_test_data"

        with patch("permucn.testdata._local_test_data_root", return_value=fake_local_root):
            with self.assertRaises(FileNotFoundError):
                fetch_test_data(dataset="toy_example", out_dir=out_dir, source="local")


if __name__ == "__main__":
    unittest.main()
