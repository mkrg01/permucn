import unittest
from fractions import Fraction
from math import comb

from permucn.stats_fisher import fisher_exact_one_sided_from_counts, fisher_min_attainable_pvalue


class TestStatsFisher(unittest.TestCase):
    def test_one_sided_fisher_matches_hypergeom_tail(self) -> None:
        p = fisher_exact_one_sided_from_counts(
            fg_concordant_count=3,
            fg_total=5,
            bg_concordant_count=2,
            bg_total=5,
        )
        self.assertAlmostEqual(p, 0.5, places=12)

    def test_min_attainable_pvalue(self) -> None:
        pmin = fisher_min_attainable_pvalue(
            fg_total=5,
            bg_total=5,
            total_concordant=5,
        )
        self.assertAlmostEqual(pmin, 1.0 / 252.0, places=12)

    def test_invalid_counts_raise(self) -> None:
        with self.assertRaises(ValueError):
            fisher_exact_one_sided_from_counts(
                fg_concordant_count=6,
                fg_total=5,
                bg_concordant_count=1,
                bg_total=5,
            )

    def test_one_sided_fisher_exhaustive_small_tables(self) -> None:
        for fg_total in range(1, 7):
            for bg_total in range(1, 7):
                total = fg_total + bg_total
                for total_concordant in range(0, total + 1):
                    lower = max(0, total_concordant - bg_total)
                    upper = min(fg_total, total_concordant)
                    for fg_concordant in range(lower, upper + 1):
                        bg_concordant = total_concordant - fg_concordant
                        p = fisher_exact_one_sided_from_counts(
                            fg_concordant_count=fg_concordant,
                            fg_total=fg_total,
                            bg_concordant_count=bg_concordant,
                            bg_total=bg_total,
                        )

                        numer = Fraction(0, 1)
                        denom = Fraction(comb(total, fg_total), 1)
                        for x in range(fg_concordant, upper + 1):
                            numer += Fraction(
                                comb(total_concordant, x) * comb(total - total_concordant, fg_total - x),
                                1,
                            )
                        p_ref = float(numer / denom)
                        self.assertAlmostEqual(p, p_ref, places=12)
                        self.assertGreaterEqual(p, -1e-12)
                        self.assertLessEqual(p, 1.0 + 1e-12)

    def test_min_attainable_matches_upper_support_probability(self) -> None:
        for fg_total in range(1, 8):
            for bg_total in range(1, 8):
                total = fg_total + bg_total
                for total_concordant in range(0, total + 1):
                    upper = min(fg_total, total_concordant)
                    pmin = fisher_min_attainable_pvalue(
                        fg_total=fg_total,
                        bg_total=bg_total,
                        total_concordant=total_concordant,
                    )
                    p_ref = fisher_exact_one_sided_from_counts(
                        fg_concordant_count=upper,
                        fg_total=fg_total,
                        bg_concordant_count=total_concordant - upper,
                        bg_total=bg_total,
                    )
                    self.assertAlmostEqual(pmin, p_ref, places=12)


if __name__ == "__main__":
    unittest.main()
