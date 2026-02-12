import tempfile
import unittest
from pathlib import Path

from permucn.permutation import PermutationGenerator
from permucn.tree import load_canonical_tree


NEXUS = """#nexus
BEGIN TREES;
  TREE t1 = ((A<0>_0:1.0,B<1>_1:1.0)<2>_0:1.0,(C<3>_1:1.0,D<4>_0:1.0)<5>_1:1.0)<6>_0;
END;
"""


class TestPermutationConstraints(unittest.TestCase):
    def _load_tree(self):
        td = tempfile.mkdtemp(prefix="permucn_tree_")
        p = Path(td) / "toy.tre"
        p.write_text(NEXUS, encoding="utf-8")
        return load_canonical_tree(p)

    def _mask(self, tree, *keys):
        m = 0
        for k in keys:
            m |= 1 << tree.branch_index_by_key[k]
        return m

    def test_non_ancestor_within_set(self) -> None:
        tree = self._load_tree()

        # Deliberately choose an observed set containing ancestor/descendant keys
        # so generator must find alternative valid samples with same bin counts.
        obs_01 = self._mask(tree, "<2>", "A<0>")

        gen = PermutationGenerator(
            tree=tree,
            obs_mask_01=obs_01,
            obs_mask_10=0,
            include_trait_loss=False,
            seed=7,
        )
        cache = gen.generate(50)

        for m01 in cache.masks_01:
            idxs = [i for i in range(len(tree.branch_key_by_index)) if (m01 >> i) & 1]
            for i in idxs:
                for j in idxs:
                    if i == j:
                        continue
                    # no ancestor/descendant relationship within same sampled set
                    self.assertEqual((tree.anc_mask_by_branch_index[i] >> j) & 1, 0)
                    self.assertEqual((tree.desc_mask_by_branch_index[i] >> j) & 1, 0)

    def test_loss_branches_descend_from_gain(self) -> None:
        tree = self._load_tree()

        # 0->1 observed on internal <2>, 1->0 observed on descendant tip A<0>
        obs_01 = self._mask(tree, "<2>")
        obs_10 = self._mask(tree, "A<0>")

        gen = PermutationGenerator(
            tree=tree,
            obs_mask_01=obs_01,
            obs_mask_10=obs_10,
            include_trait_loss=True,
            seed=11,
        )
        cache = gen.generate(80)

        for m01, m10 in zip(cache.masks_01, cache.masks_10):
            if m10 == 0:
                continue
            for idx10 in [i for i in range(len(tree.branch_key_by_index)) if (m10 >> i) & 1]:
                ok = False
                for idx01 in [i for i in range(len(tree.branch_key_by_index)) if (m01 >> i) & 1]:
                    if (tree.desc_mask_by_branch_index[idx01] >> idx10) & 1:
                        ok = True
                        break
                self.assertTrue(ok)

    def test_bin_counts_preserved(self) -> None:
        tree = self._load_tree()
        obs_01 = self._mask(tree, "<2>", "<5>")
        obs_10 = self._mask(tree, "A<0>")

        gen = PermutationGenerator(
            tree=tree,
            obs_mask_01=obs_01,
            obs_mask_10=obs_10,
            include_trait_loss=True,
            seed=17,
        )
        cache = gen.generate(60)

        expected_01 = self._bin_counts(tree, obs_01)
        expected_10 = self._bin_counts(tree, obs_10)

        for m01, m10 in zip(cache.masks_01, cache.masks_10):
            self.assertEqual(self._bin_counts(tree, m01), expected_01)
            self.assertEqual(self._bin_counts(tree, m10), expected_10)

    def test_edge_case_zero_gain_nonzero_loss(self) -> None:
        tree = self._load_tree()
        obs_01 = 0
        obs_10 = self._mask(tree, "A<0>", "C<3>")

        gen = PermutationGenerator(
            tree=tree,
            obs_mask_01=obs_01,
            obs_mask_10=obs_10,
            include_trait_loss=True,
            seed=23,
        )
        cache = gen.generate(40)

        self.assertTrue(all(m == 0 for m in cache.masks_01))
        expected_10 = self._bin_counts(tree, obs_10)
        for m10 in cache.masks_10:
            self.assertEqual(self._bin_counts(tree, m10), expected_10)

    def test_parallel_generation_matches_sequential(self) -> None:
        tree = self._load_tree()
        obs_01 = self._mask(tree, "<2>", "<5>")
        obs_10 = self._mask(tree, "A<0>")

        gen = PermutationGenerator(
            tree=tree,
            obs_mask_01=obs_01,
            obs_mask_10=obs_10,
            include_trait_loss=True,
            seed=31,
        )
        cache_seq = gen.generate(50, jobs=1)
        cache_par = gen.generate(50, jobs=2)

        self.assertEqual(cache_seq.masks_01, cache_par.masks_01)
        self.assertEqual(cache_seq.masks_10, cache_par.masks_10)

    @staticmethod
    def _bin_counts(tree, mask):
        out = [0] * 8
        m = mask
        while m:
            lsb = m & -m
            idx = lsb.bit_length() - 1
            out[tree.clade_bin_by_branch_index[idx]] += 1
            m ^= lsb
        return out


if __name__ == "__main__":
    unittest.main()
