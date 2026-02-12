import tempfile
import unittest
from pathlib import Path

from permucn.trait_ml import run_trait_asr_ml
from permucn.tree import load_canonical_tree


NEXUS = """#nexus
BEGIN TREES;
  TREE t1 = ((A<0>_0:1.0,B<1>_1:1.0)<2>_0:1.0,(C<3>_1:1.0,D<4>_0:1.0)<5>_1:1.0)<6>_0;
END;
"""


class TestSpeciesMismatch(unittest.TestCase):
    def test_species_mismatch_raises(self) -> None:
        td = tempfile.mkdtemp(prefix="permucn_tree_")
        p = Path(td) / "toy.tre"
        p.write_text(NEXUS, encoding="utf-8")
        tree = load_canonical_tree(p)

        species_to_state = {
            "A": 0,
            "B": 1,
            "C": 0,
            # D missing
            "X": 1,
        }

        with self.assertRaises(ValueError) as cm:
            run_trait_asr_ml(tree, species_to_state)
        self.assertIn("Species mismatch", str(cm.exception))


if __name__ == "__main__":
    unittest.main()
