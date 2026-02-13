import json
import tempfile
import unittest
from pathlib import Path

from permucn.cli import main


class TestUserToyData(unittest.TestCase):
    def test_test_data_binary_smoke(self) -> None:
        root = Path(__file__).resolve().parents[1]
        cafe_dir = root / "test_data" / "toy_example" / "cafe_output"
        trait_tsv = root / "test_data" / "toy_example" / "species_trait.tsv"

        self.assertTrue(cafe_dir.exists())
        self.assertTrue(trait_tsv.exists())

        td = Path(tempfile.mkdtemp(prefix="permucn_user_data_"))
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


if __name__ == "__main__":
    unittest.main()
