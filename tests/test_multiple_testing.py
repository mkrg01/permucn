import unittest
from itertools import product

from permucn.multiple_testing import bh_adjust_with_none, bonferroni_adjust_selected, tarone_screen_min_pvalues


class TestBH(unittest.TestCase):
    def test_bh_with_none(self) -> None:
        p = [0.01, 0.02, None, 0.03]
        q = bh_adjust_with_none(p)
        self.assertEqual(len(q), 4)
        self.assertIsNone(q[2])
        self.assertAlmostEqual(q[0], 0.03)
        self.assertAlmostEqual(q[1], 0.03)
        self.assertAlmostEqual(q[3], 0.03)


class TestTarone(unittest.TestCase):
    def test_tarone_screening(self) -> None:
        out = tarone_screen_min_pvalues([0.001, 0.04, 0.2, None], alpha=0.05)
        self.assertEqual(out["m_total"], 3)
        self.assertEqual(out["bonferroni_denom"], 2)
        self.assertAlmostEqual(out["threshold"], 0.025)
        self.assertEqual(out["testable_mask"], [True, False, False, False])
        self.assertEqual(out["m_testable"], 1)

    def test_bonferroni_adjust_selected(self) -> None:
        adj = bonferroni_adjust_selected(
            pvalues=[0.01, 0.02, None],
            selected_mask=[True, False, True],
            denom=2,
        )
        self.assertEqual(len(adj), 3)
        self.assertAlmostEqual(adj[0], 0.02)
        self.assertIsNone(adj[1])
        self.assertIsNone(adj[2])

    def test_tarone_minimality_and_mask_consistency_exhaustive(self) -> None:
        alpha = 0.05
        values = [None, 0.001, 0.01, 0.03, 0.2]
        eps = 1e-15

        for combo in product(values, repeat=4):
            pmins = list(combo)
            out = tarone_screen_min_pvalues(pmins, alpha=alpha)

            valid = [float(p) for p in pmins if p is not None]
            m_total = len(valid)
            self.assertEqual(out["m_total"], m_total)

            if m_total == 0:
                self.assertEqual(out["bonferroni_denom"], 0)
                self.assertIsNone(out["threshold"])
                self.assertEqual(out["testable_mask"], [False, False, False, False])
                self.assertEqual(out["m_testable"], 0)
                continue

            denom = int(out["bonferroni_denom"])
            threshold = float(out["threshold"])
            mask = list(out["testable_mask"])

            self.assertGreaterEqual(denom, 1)
            self.assertLessEqual(denom, m_total)
            self.assertAlmostEqual(threshold, alpha / denom, places=15)

            def m_k(k: int) -> int:
                cutoff = alpha / k
                return sum(1 for p in valid if p <= cutoff + eps)

            self.assertLessEqual(m_k(denom), denom)
            for k in range(1, denom):
                self.assertGreater(m_k(k), k)

            expected_mask = [
                (p is not None) and (float(p) <= threshold + eps)
                for p in pmins
            ]
            self.assertEqual(mask, expected_mask)
            self.assertEqual(out["m_testable"], sum(1 for v in mask if v))


if __name__ == "__main__":
    unittest.main()
