"""Independent physical and regression checks for the Task 10 RNA model."""
import copy
import importlib.util
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT/"projects/task10_rna_splicing/run_task10_rna_targeted_small_molecule_dynamics.py"
spec = importlib.util.spec_from_file_location("task10", PATH)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class Task10Tests(unittest.TestCase):
    def test_partition_against_exhaustive_structures(self):
        self.assertTrue(m.numerical_checks(m.defaults())["partition_enumeration_pass"])

    def test_probability_normalization_and_entropy(self):
        e = m.secondary_ensemble(m.defaults()["sequence"], 37)
        np.testing.assert_allclose(e["p"], e["p"].T, atol=1e-14)
        np.testing.assert_allclose(e["p"].sum(axis=1)+e["unpaired"], 1., atol=1e-12)
        self.assertTrue(np.all((e["p"] >= 0) & (e["p"] <= 1)))
        self.assertTrue(np.all(e["entropy"] >= 0))

    def test_unpaired_state_is_normalized_and_zero_entropy(self):
        e = m.secondary_ensemble("AAAAAAAA", 37)
        np.testing.assert_array_equal(e["unpaired"], np.ones(8))
        np.testing.assert_array_equal(e["entropy"], np.zeros(8))

    def test_pdb_identity_models_and_modified_base(self):
        for pid in ("6HMI", "6HMO"):
            models = m.read_pdb(PATH.parent/"inputs"/(pid+".pdb"))
            self.assertEqual(len(models), 20)
            self.assertEqual(len({(a["chain"], a["resid"]) for a in models[0] if a["resname"] in m.EDGES}), 22)
            self.assertTrue(any(a["resname"] == "PSU" for a in models[0]))
            self.assertTrue(any(a["chain"] == "B" and a["resid"] == 14 and a["resname"] == "A" for a in models[0]))

    def test_debye_units_and_screening(self):
        self.assertAlmostEqual(m.debye_length_angstrom(.15, 310.15), 8.01106499, places=6)
        p = np.array([[0., 0., 0.]])
        point = np.array([10., 0., 0.])
        low = m.phosphate_potential_mV(point, p, .01, 310.15)
        high = m.phosphate_potential_mV(point, p, 1., 310.15)
        self.assertLess(low, high)
        self.assertLess(high, 0.)

    def test_open_only_kd_and_full_cycle(self):
        x = m.defaults()["ligands"][0]
        p = m.ligand_parameters(x, 310.15)
        self.assertAlmostEqual(p["intrinsic_kd_M"]/p["competent_fraction"], p["open_only_kdeff_M"])
        self.assertAlmostEqual(p["cycle_residual_kcal"], 0., places=13)
        self.assertLess(p["full_kdeff_M"], p["open_only_kdeff_M"])

    def test_reciprocity_u1_and_ligand(self):
        p = m.state_weights(1e-6, 1e-7, .01, 1e-8, 1e-3, 1e-6, 100)
        self.assertAlmostEqual(p[7]*p[1]/p[3]/p[5], 100.)
        self.assertAlmostEqual(p.sum(), 1.)

    def test_finite_ligand_depletion_and_u1_mass_balance(self):
        c = m.defaults()
        s = m.make_sites(c, c["ligands"][0])
        for dose in (0., 1e-12, 1e-8, 1e-4):
            e = m.solve_equilibrium(dose, c["u1_total_M"], s)
            self.assertLess(abs(e["ligand_mass_error_M"]), 1e-16)
            self.assertLess(abs(e["u1_mass_error_M"]), 1e-16)
            self.assertLessEqual(e["free_l_M"], dose)

    def test_zero_u1_and_no_ligand_limits(self):
        c = m.defaults()
        s = m.make_sites(c, c["ligands"][0])
        e = m.solve_equilibrium(1e-3, 0., s)
        self.assertEqual(e["p"][0][4:].sum(), 0.)
        p = m.state_weights(0., 1e-7, .1, 1e-8, 1e-3, 2e-6, 100.)
        self.assertAlmostEqual(p[4:].sum(), 1/21)

    def test_analytic_saturation_matches_large_ligand_limit(self):
        c = m.defaults()
        s = m.make_sites(c, c["ligands"][0])
        exact = m.saturation_probabilities(c["u1_total_M"], s)
        finite = m.solve_equilibrium(1., c["u1_total_M"], s)["p"]
        self.assertLess(abs(exact[0][4:].sum()-finite[0][4:].sum()), 2e-6)

    def test_inactive_mutant_remains_uncoupled(self):
        c = m.defaults()
        s = m.make_sites(c, c["ligands"][2], mutant=True)
        p0 = m.solve_equilibrium(0., c["u1_total_M"], s)["p"][0][4:].sum()
        p1 = m.solve_equilibrium(.01, c["u1_total_M"], s)["p"][0][4:].sum()
        self.assertAlmostEqual(p0, p1, places=13)

    def test_selectivity_bands_do_not_bridge_failed_points(self):
        rows = [{"total_l_M": i, "inclusion_proxy_percent": p, "offtarget_inclusion_proxy_percent": 0} for i, p in enumerate([90, 0, 90, 90])]
        bands = m.contiguous_grid_bands(rows, 0)
        self.assertEqual(len(bands), 2)
        self.assertFalse(bands[0]["right_censored_at_grid_max"])
        self.assertTrue(bands[1]["right_censored_at_grid_max"])

    def test_nonmonotonic_endpoint_retains_crossing(self):
        c = m.defaults()
        c["dose_points"] = 41
        _, endpoints = m.splicing_scan(c)
        e = [e for e in endpoints if e["ligand"] == c["ligands"][3]["name"] and e["site"] == "mutant scenario"][0]
        self.assertEqual(e["response_direction_on_grid"], "nonmonotonic")
        self.assertIsNone(e["ec50_total_M"])
        self.assertGreaterEqual(len(e["half_dynamic_range_crossings"]), 1)

    def test_invalid_configuration_and_output_protection(self):
        c = m.defaults()
        c["salt_M"] = -1
        with self.assertRaises(ValueError):
            m.validate_config(c)
        with tempfile.TemporaryDirectory() as t:
            p = Path(t)
            (p/"keep.txt").write_text("user content")
            with self.assertRaises(FileExistsError):
                m.run(m.defaults(), p)
            self.assertEqual((p/"keep.txt").read_text(), "user content")


if __name__ == "__main__":
    unittest.main()
