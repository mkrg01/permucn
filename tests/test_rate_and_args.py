import argparse
import unittest

from permucn.cli import _validate_args
from permucn.stats_rate import build_rates


class TestRateAndArgs(unittest.TestCase):
    def test_build_rates_raises_on_non_positive_length(self) -> None:
        with self.assertRaises(ValueError):
            build_rates([1, -1], [1.0, 0.0])

    def test_rate_mode_rejects_cafe_significant_only(self) -> None:
        args = argparse.Namespace(
            asr_posterior_lo=0.1,
            asr_posterior_hi=0.9,
            n_perm_initial=10,
            n_perm_refine=20,
            refine_p_threshold=0.01,
            cafe_alpha=0.05,
            jobs=1,
            mode="rate",
            cafe_significant_only=True,
        )
        with self.assertRaises(ValueError):
            _validate_args(args)

    def test_negative_jobs_rejected(self) -> None:
        args = argparse.Namespace(
            asr_posterior_lo=0.1,
            asr_posterior_hi=0.9,
            n_perm_initial=10,
            n_perm_refine=20,
            refine_p_threshold=0.01,
            cafe_alpha=0.05,
            jobs=-1,
            mode="binary",
            cafe_significant_only=False,
        )
        with self.assertRaises(ValueError):
            _validate_args(args)


if __name__ == "__main__":
    unittest.main()
