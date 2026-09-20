"""Independent scientific/numerical invariants for the synthetic QSP module."""
import hashlib
import contextlib
import csv
import io
import json
from pathlib import Path
import sys
import unittest
from unittest import mock
from tempfile import TemporaryDirectory

import numpy as np
from PIL import Image
from scipy.integrate import quad, solve_ivp
from matplotlib.figure import Figure

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from task7_qsp import run_task7_qsp_tumor_immune_pkpd_synergy as qsp


class QSPTests(unittest.TestCase):
    def test_numerical_invariants(self):
        checks = qsp.numeric_checks(qsp.DEFAULT)
        self.assertGreaterEqual(len(checks), 10)
        self.assertTrue(all(c["passed"] for c in checks))

    def test_invalid_units_and_settings_rejected(self):
        for field, value in [("small_volume_l", 0), ("transit_tau_day", -1),
                             ("oral_bioavailability", 1.1), ("small_dose_mg", -50),
                             ("mab_dose_mg_kg", float("nan")), ("n_patients", 2.5),
                             ("infusion_duration_day", 14), ("progression_mass_ratio", 1),
                             ("permutation_replicates", 4.5)]:
            with self.subTest(field=field), self.assertRaises(ValueError):
                qsp.validate_config(dict(qsp.DEFAULT, **{field: value}))

    def test_mean_cv_parameterization(self):
        draw = qsp.lognormal_mean_cv(np.random.default_rng(99), 15, 0.3, 200000)
        self.assertAlmostEqual(draw.mean(), 15, delta=0.04)
        self.assertAlmostEqual(draw.std() / draw.mean(), 0.3, delta=0.002)

    def test_oral_auc_dose_over_clearance(self):
        cfg = dict(qsp.DEFAULT, horizon_day=1, small_interval_day=100)
        area = quad(lambda t: float(qsp.oral_concentration(t, 50, [15], cfg)[0]), 0, 150, epsabs=1e-9)[0]
        self.assertAlmostEqual(area, cfg["oral_bioavailability"] * 50 / 15, places=7)

    def test_infusion_total_auc(self):
        cfg = dict(qsp.DEFAULT, horizon_day=1)
        end = cfg["infusion_duration_day"]
        area = quad(lambda t: float(qsp.antibody_concentration(t, 10, cfg)), 0, end)[0]
        area += quad(lambda t: float(qsp.antibody_concentration(t, 10, cfg)), end, 600)[0]
        self.assertAlmostEqual(area, 10 * cfg["body_weight_kg"] / cfg["mab_clearance_l_day"], places=6)

    def test_transit_chain_erlang_three_not_tau_total(self):
        cfg = dict(qsp.DEFAULT, immune_kill_ctl_day_inv=0, antigen_stim_ctl_day=0, baseline_ctl=0)
        cohort = qsp.make_cohort(cfg, 1, nominal=True)
        tau = cfg["transit_tau_day"]
        times = np.linspace(0, 30, 301)
        solution = solve_ivp(qsp.rhs, (0, 30), [0, 10, 0, 0, 0],
                             args=(cohort, 0, 0, cfg), t_eval=times, rtol=1e-10, atol=1e-12)
        expected = 10 * np.exp(-times / tau) * (1 + times / tau + 0.5 * (times / tau) ** 2)
        np.testing.assert_allclose(solution.y[1:4].sum(axis=0), expected, rtol=1e-7, atol=1e-9)
        self.assertGreater(solution.y[1:4, 15].sum(), 9)  # At one tau, most mass remains.

    def test_paired_cohort_reproducibility(self):
        a, b = qsp.make_cohort(qsp.DEFAULT), qsp.make_cohort(qsp.DEFAULT)
        for name in a:
            np.testing.assert_array_equal(a[name], b[name])

    def test_valid_custom_baseline_and_zero_baseline(self):
        for baseline in [0.0, 0.5, 3.0]:
            cfg = dict(qsp.DEFAULT, baseline_ctl=baseline)
            qsp.validate_config(cfg)
            self.assertTrue(all(c["passed"] for c in qsp.numeric_checks(cfg)))
            cohort = qsp.make_cohort(cfg)
            self.assertTrue(np.all(np.isfinite(cohort["ctl_baseline"])))
            if baseline == 0:
                np.testing.assert_array_equal(cohort["ctl_baseline"], 0)

    def test_loewe_exact_internal_knot_and_axes(self):
        inverse, status = qsp.invert_monotherapy([0, 1, 2], [0, 0.4, 0.8], 0.4)
        self.assertEqual(inverse, 1)
        ci, status = qsp.loewe_index(0, 1, 0.4, [0, 1], [0, 0.1], [0, 1, 2], [0, 0.4, 0.8])
        self.assertEqual(ci, 1)
        inverse, status = qsp.invert_monotherapy([0, 1, 2], [0, 0.5, 0.4], 0.3)
        self.assertIsNone(inverse)

    def test_cox_does_not_emit_zero_hr_for_no_combo_events(self):
        t = np.arange(1.0, 6.0)
        result = qsp.paired_survival_statistics((t, np.ones(5, dtype=bool)),
                 (np.full(5, 60.0), np.zeros(5, dtype=bool)), qsp.DEFAULT)
        self.assertIsNone(result["cox_hr_combo_vs_mono"])
        self.assertEqual(result["cox_status"], "not_estimable")

    def test_cox_identical_paired_observations_hr_one(self):
        t = np.arange(1.0, 6.0)
        ep = (t, np.ones(5, dtype=bool))
        result = qsp.paired_survival_statistics(ep, ep, qsp.DEFAULT)
        self.assertAlmostEqual(result["cox_hr_combo_vs_mono"], 1, places=8)
        self.assertEqual(result["paired_permutation_p"], 1)

    def test_custom_ten_people_thirty_days_outputs_and_figure_labels(self):
        cfg = dict(qsp.DEFAULT, n_patients=10, horizon_day=30.0)
        captions = []
        original_savefig = Figure.savefig
        def record_and_save(figure, *args, **kwargs):
            captions.extend(ax.get_title() for ax in figure.axes)
            return original_savefig(figure, *args, **kwargs)
        with TemporaryDirectory(prefix="qsp_custom_") as folder:
            out = Path(folder) / "run"
            with mock.patch.object(Figure, "savefig", new=record_and_save), contextlib.redirect_stdout(io.StringIO()):
                qsp.run(cfg, out)
            summary = json.loads((out / "results_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["n_unique_virtual_people"], 10)
            self.assertEqual(summary["endpoint_day"], 30)
            self.assertEqual(summary["trajectory_rows"], 10 * 4 * 121)
            self.assertNotIn("day60", json.dumps(summary))
            with (out / "loewe_monotherapy_support.csv").open(encoding="utf-8") as handle:
                support = list(csv.DictReader(handle))
            self.assertTrue(all(float(row["endpoint_day"]) == 30 for row in support))
            self.assertIn("endpoint_effect", support[0])
            with (out / "progression_endpoints.csv").open(encoding="utf-8") as handle:
                endpoints = list(csv.DictReader(handle))
            self.assertEqual(len(endpoints), 40)
            self.assertTrue(all(float(row["duration_day"]) <= 30 for row in endpoints))
            self.assertEqual(len(list((out / "figures_task7").glob("*.png"))), 4)
        self.assertIn("Synthetic QSP cohort: mean tumor burden ± SEM (10 paired people)", captions)
        self.assertIn("Day-30 nominal-person Bliss landscape", captions)
        self.assertIn("10 paired virtual people; administrative censoring at day 30", captions)

    def test_delivered_scientific_artifacts_and_dpi(self):
        folder = ROOT / "task7_qsp"
        manifest = json.loads((folder / "manifest.json").read_text(encoding="utf-8"))
        for relative, digest in manifest["files_sha256"].items():
            with self.subTest(file=relative):
                self.assertEqual(hashlib.sha256((folder / relative).read_bytes()).hexdigest(), digest)
        self.assertEqual(hashlib.sha256(Path(qsp.__file__).read_bytes()).hexdigest(), manifest["script_sha256"])
        figures = list((folder / "figures_task7").glob("*.png"))
        self.assertEqual(len(figures), 4)
        for file in figures:
            with Image.open(file) as image:
                self.assertTrue(all(abs(dpi - 300) < 0.02 for dpi in image.info["dpi"]))
        summary = json.loads((folder / "results_summary.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["n_unique_virtual_people"], 50)
        self.assertFalse(summary["clinical_validation"])
        self.assertTrue(all(c["passed"] for c in json.loads((folder / "verification.json").read_text(encoding="utf-8"))["checks"]))


if __name__ == "__main__":
    unittest.main()
