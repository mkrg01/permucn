import tempfile
import unittest
from pathlib import Path

from permucn.trait_ml import run_trait_asr_ml
from permucn.tree import CanonicalTree, load_canonical_tree


NEXUS = """#nexus
BEGIN TREES;
  TREE t1 = ((A<0>_0:1.0,B<1>_1:1.0)<2>_0:1.0,(C<3>_1:1.0,D<4>_0:1.0)<5>_1:1.0)<6>_0;
END;
"""


class TestTraitMlLogic(unittest.TestCase):
    def _load_tree(self) -> CanonicalTree:
        td = tempfile.mkdtemp(prefix="permucn_trait_ml_")
        p = Path(td) / "toy.tre"
        p.write_text(NEXUS, encoding="utf-8")
        return load_canonical_tree(p)

    @staticmethod
    def _mask_to_keys(tree: CanonicalTree, mask: int) -> set[str]:
        return {
            key
            for idx, key in enumerate(tree.branch_key_by_index)
            if (mask >> idx) & 1
        }

    def test_detects_single_gain_branch(self) -> None:
        tree = self._load_tree()
        species_to_state = {"A": 0, "B": 0, "C": 0, "D": 1}

        asr = run_trait_asr_ml(tree, species_to_state, posterior_lo=0.4, posterior_hi=0.6)

        self.assertSetEqual(self._mask_to_keys(tree, asr.fg_01_mask), {"D<4>"})
        self.assertEqual(asr.n_fg_01, 1)
        self.assertEqual(asr.fg_10_mask, 0)
        self.assertEqual(asr.n_fg_10, 0)

        for node_id, species in enumerate(tree.tip_species_by_node):
            if species is None:
                continue
            self.assertEqual(asr.hard_state_by_node[node_id], species_to_state[species])

    def test_detects_single_loss_branch(self) -> None:
        tree = self._load_tree()
        species_to_state = {"A": 1, "B": 1, "C": 1, "D": 0}

        asr = run_trait_asr_ml(tree, species_to_state, posterior_lo=0.4, posterior_hi=0.6)

        self.assertEqual(asr.fg_01_mask, 0)
        self.assertEqual(asr.n_fg_01, 0)
        self.assertSetEqual(self._mask_to_keys(tree, asr.fg_10_mask), {"D<4>"})
        self.assertEqual(asr.n_fg_10, 1)

    def test_threshold_boundaries_are_inclusive_at_half(self) -> None:
        tree = self._load_tree()
        species_to_state = {"A": 0, "B": 1, "C": 0, "D": 1}

        # p(state=1)==0.5 at internal nodes for this symmetric tip pattern.
        asr_hi_inclusive = run_trait_asr_ml(tree, species_to_state, posterior_lo=0.49, posterior_hi=0.5)
        self.assertSetEqual(self._mask_to_keys(tree, asr_hi_inclusive.fg_10_mask), {"A<0>", "C<3>"})
        self.assertEqual(asr_hi_inclusive.n_fg_10, 2)
        self.assertEqual(asr_hi_inclusive.fg_01_mask, 0)
        for node_id in (2, 5, 6):
            self.assertEqual(asr_hi_inclusive.hard_state_by_node[node_id], 1)

        asr_lo_inclusive = run_trait_asr_ml(tree, species_to_state, posterior_lo=0.5, posterior_hi=0.51)
        self.assertSetEqual(self._mask_to_keys(tree, asr_lo_inclusive.fg_01_mask), {"B<1>", "D<4>"})
        self.assertEqual(asr_lo_inclusive.n_fg_01, 2)
        self.assertEqual(asr_lo_inclusive.fg_10_mask, 0)
        for node_id in (2, 5, 6):
            self.assertEqual(asr_lo_inclusive.hard_state_by_node[node_id], 0)

        asr_ambiguous = run_trait_asr_ml(tree, species_to_state, posterior_lo=0.4999, posterior_hi=0.5001)
        self.assertEqual(asr_ambiguous.fg_01_mask, 0)
        self.assertEqual(asr_ambiguous.fg_10_mask, 0)
        self.assertEqual(asr_ambiguous.n_fg_01, 0)
        self.assertEqual(asr_ambiguous.n_fg_10, 0)
        for node_id in (2, 5, 6):
            self.assertIsNone(asr_ambiguous.hard_state_by_node[node_id])


if __name__ == "__main__":
    unittest.main()
