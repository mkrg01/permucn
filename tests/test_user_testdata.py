import json
import tempfile
import unittest
from pathlib import Path

from permucn.cli import main


class TestUserToyData(unittest.TestCase):
    def test_get_test_data_binary_smoke(self) -> None:
        td = Path(tempfile.mkdtemp(prefix="permucn_user_data_"))
        test_data_dir = td / "fetched_test_data"

        rc_fetch = main(
            [
                "get-test-data",
                "--dataset",
                "toy_example",
                "--source",
                "local",
                "--out-dir",
                str(test_data_dir),
            ]
        )
        self.assertEqual(rc_fetch, 0)

        cafe_dir = test_data_dir / "toy_example" / "cafe_output"
        trait_tsv = test_data_dir / "toy_example" / "species_trait.tsv"
        self.assertTrue(cafe_dir.exists())
        self.assertTrue(trait_tsv.exists())

        out_prefix = td / "out" / "toy_binary"

        rc = main(
            [
                "--mode",
                "binary",
                "--cafe-dir",
                str(cafe_dir),
                "--trait-tsv",
                str(trait_tsv),
                "--no-include-trait-loss",
                "--asr-posterior-hi",
                "0.6",
                "--asr-posterior-lo",
                "0.4",
                "--n-perm-initial",
                "5",
                "--n-perm-refine",
                "10",
                "--seed",
                "7",
                "--out-prefix",
                str(out_prefix),
            ]
        )
        self.assertEqual(rc, 0)

        out_tsv = Path(str(out_prefix) + ".family_results.tsv")
        out_json = Path(str(out_prefix) + ".run_metadata.json")
        self.assertTrue(out_tsv.exists())
        self.assertTrue(out_json.exists())

        meta = json.loads(out_json.read_text(encoding="utf-8"))
        self.assertEqual(meta["parameters"]["mode"], "binary")
        self.assertEqual(meta["trait_columns"]["trait_column_source"], "auto")
        self.assertGreater(meta["results"]["n_tested"], 0)

    def test_get_test_data_overwrite_requires_force(self) -> None:
        td = Path(tempfile.mkdtemp(prefix="permucn_user_data_overwrite_"))
        test_data_dir = td / "fetched_test_data"

        first = main(
            [
                "get-test-data",
                "--dataset",
                "toy_example",
                "--source",
                "local",
                "--out-dir",
                str(test_data_dir),
            ]
        )
        self.assertEqual(first, 0)

        second = main(
            [
                "get-test-data",
                "--dataset",
                "toy_example",
                "--source",
                "local",
                "--out-dir",
                str(test_data_dir),
            ]
        )
        self.assertEqual(second, 1)

        third = main(
            [
                "get-test-data",
                "--dataset",
                "toy_example",
                "--source",
                "local",
                "--out-dir",
                str(test_data_dir),
                "--force",
            ]
        )
        self.assertEqual(third, 0)


if __name__ == "__main__":
    unittest.main()
