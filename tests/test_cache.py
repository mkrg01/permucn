import tempfile
import unittest
from pathlib import Path

from permucn.cache import (
    empty_bundle,
    get_stage_cache,
    is_bundle_compatible,
    make_cache_spec,
    put_stage_cache,
    save_cache_bundle,
    load_cache_bundle,
)
from permucn.permutation import PermutationCache
from permucn.tree import load_canonical_tree


NEXUS = """#nexus
BEGIN TREES;
  TREE t1 = ((A<0>_0:1.0,B<1>_1:1.0)<2>_0:1.0,(C<3>_1:1.0,D<4>_0:1.0)<5>_1:1.0)<6>_0;
END;
"""


class TestCache(unittest.TestCase):
    def _tree(self):
        td = tempfile.mkdtemp(prefix="permucn_cache_tree_")
        p = Path(td) / "toy.tre"
        p.write_text(NEXUS, encoding="utf-8")
        return load_canonical_tree(p)

    def test_cache_roundtrip_and_slice(self) -> None:
        tree = self._tree()
        spec = make_cache_spec(tree, include_trait_loss=True, fg_01_mask=0b1010, fg_10_mask=0b0100)
        bundle = empty_bundle(spec)

        cache = PermutationCache(
            masks_01=[0b1, 0b10, 0b100],
            masks_10=[0b0, 0b1, 0b0],
            total_attempts=12,
            total_restarts=3,
        )
        put_stage_cache(bundle, "initial", cache)

        td = tempfile.mkdtemp(prefix="permucn_cache_io_")
        path = Path(td) / "perm_cache.json.gz"
        save_cache_bundle(path, bundle)

        loaded = load_cache_bundle(path)
        self.assertTrue(is_bundle_compatible(loaded, spec))

        got2 = get_stage_cache(loaded, "initial", 2)
        self.assertIsNotNone(got2)
        self.assertEqual(got2.masks_01, [0b1, 0b10])
        self.assertEqual(got2.masks_10, [0b0, 0b1])

        got4 = get_stage_cache(loaded, "initial", 4)
        self.assertIsNone(got4)


if __name__ == "__main__":
    unittest.main()
