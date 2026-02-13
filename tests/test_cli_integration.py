import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

from permucn.cli import main


TOY_NEXUS = """#nexus
BEGIN TREES;
  TREE t1 = ((A<0>_0:1.0,B<1>_1:1.0)<2>_0:1.0,(C<3>_1:1.0,D<4>_0:1.0)<5>_1:1.0)<6>_0;
END;
"""


class TestCliIntegration(unittest.TestCase):
    def _build_toy_inputs(self, trait_rows: list[list[str]] | None = None) -> tuple[Path, Path, Path]:
        td = Path(tempfile.mkdtemp(prefix="permucn_cli_"))
        cafe_dir = td / "cafe_output"
        cafe_dir.mkdir(parents=True, exist_ok=True)

        # Canonical tree
        (cafe_dir / "Gamma_asr.tre").write_text(TOY_NEXUS, encoding="utf-8")

        # Change table: non-root branches only (tips + internal branches)
        change_header = ["FamilyID", "A<0>", "B<1>", "<2>", "C<3>", "D<4>", "<5>"]
        change_rows = [
            ["fam1", "+1", "0", "+1", "0", "0", "-1"],
            ["fam2", "0", "-1", "0", "+1", "0", "+1"],
            ["fam3", "0", "0", "0", "0", "0", "0"],
        ]
        self._write_tsv(cafe_dir / "Gamma_change.tab", change_header, change_rows)

        prob_header = ["#FamilyID", "A<0>", "B<1>", "<2>", "C<3>", "D<4>", "<5>"]
        prob_rows = [
            ["fam1", "0.01", "0.80", "0.02", "0.90", "0.50", "0.03"],
            ["fam2", "0.80", "0.02", "0.70", "0.01", "0.90", "0.04"],
            ["fam3", "0.90", "0.90", "0.90", "0.90", "0.90", "0.90"],
        ]
        self._write_tsv(cafe_dir / "Gamma_branch_probabilities.tab", prob_header, prob_rows)

        (cafe_dir / "Gamma_family_results.txt").write_text(
            "#FamilyID\tpvalue\tSignificant at 0.05\n"
            "fam1\t0.2\tn\n"
            "fam2\t0.1\tn\n"
            "fam3\t1.0\tn\n",
            encoding="utf-8",
        )

        # Trait TSV
        trait_tsv = td / "species_trait.tsv"
        if trait_rows is None:
            trait_rows = [["A", "0"], ["B", "0"], ["C", "0"], ["D", "1"]]
        self._write_tsv(
            trait_tsv,
            ["species", "polar"],
            trait_rows,
        )

        return td, cafe_dir, trait_tsv

    @staticmethod
    def _write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
        lines = ["\t".join(header)]
        lines.extend("\t".join(r) for r in rows)
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    def test_binary_run_end_to_end(self) -> None:
        td, cafe_dir, trait_tsv = self._build_toy_inputs()
        out_prefix = td / "out" / "binary"

        rc = main(
            [
                "--mode",
                "binary",
                "--cafe-dir",
                str(cafe_dir),
                "--trait-tsv",
                str(trait_tsv),
                "--no-include-trait-loss",
                "--asr-posterior-hi",
                "0.6",
                "--asr-posterior-lo",
                "0.4",
                "--n-perm-initial",
                "5",
                "--n-perm-refine",
                "10",
                "--refine-p-threshold",
                "0.01",
                "--seed",
                "7",
                "--out-prefix",
                str(out_prefix),
            ]
        )
        self.assertEqual(rc, 0)

        out_tsv = Path(str(out_prefix) + ".family_results.tsv")
        out_json = Path(str(out_prefix) + ".run_metadata.json")
        self.assertTrue(out_tsv.exists())
        self.assertTrue(out_json.exists())

        header = out_tsv.read_text(encoding="utf-8").splitlines()[0].split("\t")
        self.assertIn("family_id", header)
        self.assertIn("p_empirical", header)
        self.assertIn("fg_concordant_count", header)

        meta = json.loads(out_json.read_text(encoding="utf-8"))
        self.assertEqual(meta["parameters"]["mode"], "binary")
        self.assertEqual(meta["trait_columns"]["trait_column_source"], "auto")
        self.assertGreater(meta["results"]["n_tested"], 0)

    def test_progress_logs_are_emitted(self) -> None:
        td, cafe_dir, trait_tsv = self._build_toy_inputs()
        out_prefix = td / "out" / "progress"
        stderr = io.StringIO()

        with contextlib.redirect_stderr(stderr):
            rc = main(
                [
                    "--mode",
                    "binary",
                    "--cafe-dir",
                    str(cafe_dir),
                    "--trait-tsv",
                    str(trait_tsv),
                    "--no-include-trait-loss",
                    "--n-perm-initial",
                    "5",
                    "--n-perm-refine",
                    "10",
                    "--seed",
                    "7",
                    "--out-prefix",
                    str(out_prefix),
                ]
            )
        self.assertEqual(rc, 0)
        logs = stderr.getvalue()
        self.assertIn("[permucn] [1/8]", logs)
        self.assertIn("[permucn] [5/8]", logs)
        self.assertIn("[permucn] [8/8] Run complete; outputs were written successfully", logs)

    def test_rate_run_end_to_end(self) -> None:
        td, cafe_dir, trait_tsv = self._build_toy_inputs()
        out_prefix = td / "out" / "rate"

        rc = main(
            [
                "--mode",
                "rate",
                "--cafe-dir",
                str(cafe_dir),
                "--trait-tsv",
                str(trait_tsv),
                "--no-include-trait-loss",
                "--asr-posterior-hi",
                "0.6",
                "--asr-posterior-lo",
                "0.4",
                "--n-perm-initial",
                "5",
                "--n-perm-refine",
                "10",
                "--refine-p-threshold",
                "0.01",
                "--seed",
                "7",
                "--out-prefix",
                str(out_prefix),
            ]
        )
        self.assertEqual(rc, 0)

        out_tsv = Path(str(out_prefix) + ".family_results.tsv")
        out_json = Path(str(out_prefix) + ".run_metadata.json")
        self.assertTrue(out_tsv.exists())
        self.assertTrue(out_json.exists())

        header = out_tsv.read_text(encoding="utf-8").splitlines()[0].split("\t")
        self.assertIn("family_id", header)
        self.assertIn("stat_obs", header)
        self.assertIn("fg_mean_signed_rate", header)

        meta = json.loads(out_json.read_text(encoding="utf-8"))
        self.assertEqual(meta["parameters"]["mode"], "rate")
        self.assertGreater(meta["results"]["n_tested"], 0)

    def test_binary_cafe_significant_only(self) -> None:
        td, cafe_dir, trait_tsv = self._build_toy_inputs()
        out_prefix = td / "out" / "binary_sig"

        rc = main(
            [
                "--mode",
                "binary",
                "--cafe-dir",
                str(cafe_dir),
                "--trait-tsv",
                str(trait_tsv),
                "--cafe-significant-only",
                "--no-include-trait-loss",
                "--asr-posterior-hi",
                "0.6",
                "--asr-posterior-lo",
                "0.4",
                "--n-perm-initial",
                "5",
                "--n-perm-refine",
                "10",
                "--refine-p-threshold",
                "0.01",
                "--seed",
                "7",
                "--out-prefix",
                str(out_prefix),
            ]
        )
        self.assertEqual(rc, 0)

        out_tsv = Path(str(out_prefix) + ".family_results.tsv")
        self.assertTrue(out_tsv.exists())

    def test_no_valid_foreground_status(self) -> None:
        td, cafe_dir, trait_tsv = self._build_toy_inputs(
            trait_rows=[["A", "0"], ["B", "0"], ["C", "0"], ["D", "0"]]
        )
        out_prefix = td / "out" / "no_fg"

        rc = main(
            [
                "--mode",
                "binary",
                "--cafe-dir",
                str(cafe_dir),
                "--trait-tsv",
                str(trait_tsv),
                "--n-perm-initial",
                "5",
                "--n-perm-refine",
                "10",
                "--seed",
                "7",
                "--out-prefix",
                str(out_prefix),
            ]
        )
        self.assertEqual(rc, 0)

        out_tsv = Path(str(out_prefix) + ".family_results.tsv")
        lines = out_tsv.read_text(encoding="utf-8").splitlines()
        self.assertGreaterEqual(len(lines), 2)
        status_idx = lines[0].split("\t").index("status")
        for row in lines[1:]:
            self.assertEqual(row.split("\t")[status_idx], "no_valid_foreground")

        meta = json.loads(Path(str(out_prefix) + ".run_metadata.json").read_text(encoding="utf-8"))
        self.assertEqual(meta["results"]["n_tested"], 0)
        self.assertEqual(meta["permutation"]["cache"]["initial_source"], "skipped_no_foreground")
        self.assertEqual(meta["permutation"]["cache"]["refine_source"], "skipped_no_foreground")
        self.assertEqual(meta["permutation"]["initial"]["n_perm"], 0)

    def test_warning_when_potential_transition_branches_are_skipped_by_thresholding(self) -> None:
        td, cafe_dir, trait_tsv = self._build_toy_inputs(
            trait_rows=[["A", "0"], ["B", "1"], ["C", "0"], ["D", "1"]]
        )
        out_prefix = td / "out" / "asr_warning"
        stderr = io.StringIO()

        with contextlib.redirect_stderr(stderr):
            rc = main(
                [
                    "--mode",
                    "binary",
                    "--cafe-dir",
                    str(cafe_dir),
                    "--trait-tsv",
                    str(trait_tsv),
                    "--asr-posterior-hi",
                    "0.6",
                    "--asr-posterior-lo",
                    "0.4",
                    "--n-perm-initial",
                    "5",
                    "--n-perm-refine",
                    "10",
                    "--seed",
                    "7",
                    "--out-prefix",
                    str(out_prefix),
                ]
            )
        self.assertEqual(rc, 0)
        logs = stderr.getvalue()
        self.assertIn("[permucn][warning]", logs)
        self.assertIn("ASR posterior thresholding skipped potential phenotype-transition branches", logs)

    def test_reproducible_outputs_with_seed(self) -> None:
        td, cafe_dir, trait_tsv = self._build_toy_inputs()
        out_prefix1 = td / "out" / "rep1"
        out_prefix2 = td / "out" / "rep2"

        argv = [
            "--mode",
            "binary",
            "--cafe-dir",
            str(cafe_dir),
            "--trait-tsv",
            str(trait_tsv),
            "--no-include-trait-loss",
            "--asr-posterior-hi",
            "0.6",
            "--asr-posterior-lo",
            "0.4",
            "--n-perm-initial",
            "7",
            "--n-perm-refine",
            "11",
            "--refine-p-threshold",
            "0.5",
            "--seed",
            "123",
        ]

        rc1 = main(argv + ["--out-prefix", str(out_prefix1)])
        rc2 = main(argv + ["--out-prefix", str(out_prefix2)])
        self.assertEqual(rc1, 0)
        self.assertEqual(rc2, 0)

        tsv1 = Path(str(out_prefix1) + ".family_results.tsv").read_text(encoding="utf-8")
        tsv2 = Path(str(out_prefix2) + ".family_results.tsv").read_text(encoding="utf-8")
        self.assertEqual(tsv1, tsv2)

    def test_permutation_cache_reuse(self) -> None:
        td, cafe_dir, trait_tsv = self._build_toy_inputs()
        cache_path = td / "perm_cache.json.gz"
        out_prefix1 = td / "out" / "cache1"
        out_prefix2 = td / "out" / "cache2"

        argv = [
            "--mode",
            "binary",
            "--cafe-dir",
            str(cafe_dir),
            "--trait-tsv",
            str(trait_tsv),
            "--no-include-trait-loss",
            "--asr-posterior-hi",
            "0.6",
            "--asr-posterior-lo",
            "0.4",
            "--n-perm-initial",
            "7",
            "--n-perm-refine",
            "7",
            "--seed",
            "101",
            "--perm-cache",
            str(cache_path),
        ]

        rc1 = main(argv + ["--out-prefix", str(out_prefix1)])
        self.assertEqual(rc1, 0)
        self.assertTrue(cache_path.exists())

        meta1 = json.loads(Path(str(out_prefix1) + ".run_metadata.json").read_text(encoding="utf-8"))
        self.assertEqual(meta1["permutation"]["cache"]["initial_source"], "generated")

        rc2 = main(argv + ["--out-prefix", str(out_prefix2)])
        self.assertEqual(rc2, 0)
        meta2 = json.loads(Path(str(out_prefix2) + ".run_metadata.json").read_text(encoding="utf-8"))
        self.assertTrue(meta2["permutation"]["cache"]["cache_loaded"])
        self.assertEqual(meta2["permutation"]["cache"]["initial_source"], "cache")

    def test_visual_summary_outputs_exist(self) -> None:
        td, cafe_dir, trait_tsv = self._build_toy_inputs()
        out_prefix = td / "out" / "viz"

        rc = main(
            [
                "--mode",
                "binary",
                "--cafe-dir",
                str(cafe_dir),
                "--trait-tsv",
                str(trait_tsv),
                "--no-include-trait-loss",
                "--asr-posterior-hi",
                "0.6",
                "--asr-posterior-lo",
                "0.4",
                "--n-perm-initial",
                "5",
                "--n-perm-refine",
                "10",
                "--pvalue-top-n",
                "2",
                "--seed",
                "9",
                "--out-prefix",
                str(out_prefix),
            ]
        )
        self.assertEqual(rc, 0)

        top_hits = Path(str(out_prefix) + ".top_hits.tsv")
        top_pvalues = Path(str(out_prefix) + ".top_pvalues.tsv")
        hist_tsv = Path(str(out_prefix) + ".pvalue_hist.tsv")
        qq_tsv = Path(str(out_prefix) + ".qq.tsv")
        self.assertTrue(top_hits.exists())
        self.assertTrue(top_pvalues.exists())
        self.assertTrue(hist_tsv.exists())
        self.assertTrue(qq_tsv.exists())

        top_lines = top_pvalues.read_text(encoding="utf-8").splitlines()
        self.assertGreaterEqual(len(top_lines), 1)
        self.assertLessEqual(len(top_lines), 3)  # header + top 2 rows
        if len(top_lines) > 2:
            header = top_lines[0].split("\t")
            p_idx = header.index("p_empirical")
            pvals = [float(row.split("\t")[p_idx]) for row in top_lines[1:]]
            self.assertEqual(pvals, sorted(pvals))

        meta = json.loads(Path(str(out_prefix) + ".run_metadata.json").read_text(encoding="utf-8"))
        viz = meta["results"]["visual_outputs"]
        self.assertEqual(viz["top_hits_tsv"], str(top_hits))
        self.assertEqual(viz["top_pvalues_tsv"], str(top_pvalues))
        self.assertEqual(viz["pvalue_hist_tsv"], str(hist_tsv))
        self.assertEqual(viz["qq_tsv"], str(qq_tsv))
        self.assertIn("pvalue_hist_pdf", viz)
        self.assertIn("qq_pdf", viz)
        self.assertNotIn("pvalue_hist_png", viz)
        self.assertNotIn("qq_png", viz)

    def test_parallel_matches_sequential(self) -> None:
        td, cafe_dir, trait_tsv = self._build_toy_inputs()
        out_prefix1 = td / "out" / "seq"
        out_prefix2 = td / "out" / "par"

        base = [
            "--mode",
            "binary",
            "--cafe-dir",
            str(cafe_dir),
            "--trait-tsv",
            str(trait_tsv),
            "--no-include-trait-loss",
            "--asr-posterior-hi",
            "0.6",
            "--asr-posterior-lo",
            "0.4",
            "--n-perm-initial",
            "7",
            "--n-perm-refine",
            "11",
            "--refine-p-threshold",
            "0.5",
            "--seed",
            "42",
        ]

        rc1 = main(base + ["--jobs", "1", "--out-prefix", str(out_prefix1)])
        rc2 = main(base + ["--jobs", "2", "--out-prefix", str(out_prefix2)])
        self.assertEqual(rc1, 0)
        self.assertEqual(rc2, 0)

        tsv1 = Path(str(out_prefix1) + ".family_results.tsv").read_text(encoding="utf-8")
        tsv2 = Path(str(out_prefix2) + ".family_results.tsv").read_text(encoding="utf-8")
        self.assertEqual(tsv1, tsv2)


if __name__ == "__main__":
    unittest.main()
