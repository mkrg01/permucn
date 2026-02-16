import argparse
import unittest

from permucn.cli import _validate_args, build_parser, build_test_data_parser
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

    def test_negative_pvalue_top_n_rejected(self) -> None:
        args = argparse.Namespace(
            asr_posterior_lo=0.1,
            asr_posterior_hi=0.9,
            n_perm_initial=10,
            n_perm_refine=20,
            refine_p_threshold=0.01,
            cafe_alpha=0.05,
            pvalue_top_n=-1,
            jobs=1,
            mode="binary",
            cafe_significant_only=False,
        )
        with self.assertRaises(ValueError):
            _validate_args(args)

    def test_rate_mode_rejects_fisher_tarone(self) -> None:
        args = argparse.Namespace(
            asr_posterior_lo=0.1,
            asr_posterior_hi=0.9,
            n_perm_initial=10,
            n_perm_refine=20,
            refine_p_threshold=0.01,
            cafe_alpha=0.05,
            jobs=1,
            mode="rate",
            binary_test="fisher-tarone",
            cafe_significant_only=False,
        )
        with self.assertRaises(ValueError):
            _validate_args(args)

    def test_main_help_shows_default_values(self) -> None:
        help_text = build_parser().format_help()
        self.assertIn("--fwer-alpha FWER_ALPHA", help_text)
        self.assertIn("(default: 0.05)", help_text)

    def test_fwer_alpha_default_is_0_05(self) -> None:
        args = build_parser().parse_args([])
        self.assertEqual(args.fwer_alpha, 0.05)

    def test_get_test_data_help_shows_default_values(self) -> None:
        help_text = build_test_data_parser().format_help()
        self.assertIn("--dataset {toy_example,polar_fish,all}", help_text)
        self.assertIn("(default: toy_example)", help_text)


if __name__ == "__main__":
    unittest.main()
