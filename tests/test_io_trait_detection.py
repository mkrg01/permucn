import tempfile
import unittest
from pathlib import Path

from permucn.io import load_trait_table


class TestTraitDetection(unittest.TestCase):
    def _write(self, text: str) -> Path:
        td = tempfile.mkdtemp(prefix="permucn_test_")
        path = Path(td) / "trait.tsv"
        path.write_text(text, encoding="utf-8")
        return path

    def test_unique_auto_detect(self) -> None:
        p = self._write(
            "species\tpolar\n"
            "sp1\t0\n"
            "sp2\t1\n"
        )
        tbl = load_trait_table(p)
        self.assertEqual(tbl.species_column, "species")
        self.assertEqual(tbl.trait_column, "polar")
        self.assertEqual(tbl.trait_column_source, "auto")
        self.assertEqual(tbl.species_to_state["sp1"], 0)
        self.assertEqual(tbl.species_to_state["sp2"], 1)

    def test_multiple_binary_columns_fail_fast(self) -> None:
        p = self._write(
            "species\ta\tb\n"
            "sp1\t0\t1\n"
            "sp2\t1\t0\n"
        )
        with self.assertRaises(ValueError) as cm:
            load_trait_table(p)
        self.assertIn("Multiple binary trait columns", str(cm.exception))

    def test_no_binary_columns_fail_fast(self) -> None:
        p = self._write(
            "species\ttrait\n"
            "sp1\t2\n"
            "sp2\t1\n"
        )
        with self.assertRaises(ValueError) as cm:
            load_trait_table(p)
        self.assertIn("No binary trait column", str(cm.exception))

    def test_manual_override(self) -> None:
        p = self._write(
            "species\ta\tb\n"
            "sp1\t0\t1\n"
            "sp2\t1\t0\n"
        )
        tbl = load_trait_table(p, trait_column="a")
        self.assertEqual(tbl.trait_column, "a")
        self.assertEqual(tbl.trait_column_source, "manual")


if __name__ == "__main__":
    unittest.main()
