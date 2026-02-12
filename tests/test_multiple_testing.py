import unittest

from permucn.multiple_testing import bh_adjust_with_none


class TestBH(unittest.TestCase):
    def test_bh_with_none(self) -> None:
        p = [0.01, 0.02, None, 0.03]
        q = bh_adjust_with_none(p)
        self.assertEqual(len(q), 4)
        self.assertIsNone(q[2])
        self.assertAlmostEqual(q[0], 0.03)
        self.assertAlmostEqual(q[1], 0.03)
        self.assertAlmostEqual(q[3], 0.03)


if __name__ == "__main__":
    unittest.main()
