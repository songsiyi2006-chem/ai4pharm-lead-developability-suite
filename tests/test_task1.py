"""Scientific invariants and an offline end-to-end artifact check."""
import copy
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from PIL import Image
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import run_task1_mpo_admet_developability as task


class NumericalTests(unittest.TestCase):
    def test_original_mpo_best_and_worst(self):
        self.assertEqual(task.cns_mpo(3, 2, 360, 40, .5, 8)["cns_mpo_approx"], 6)
        self.assertEqual(task.cns_mpo(5, 4, 500, 120, 3.5, 10)["cns_mpo_approx"], 0)
        self.assertEqual(task.cns_mpo(3, 2, 360, 90, .5, None)["cns_mpo_approx"], 6)

    def test_interpolation_and_continuity(self):
        terms = task.cns_mpo(4, 3, 430, 105, 2, 9)
        self.assertAlmostEqual(terms["cns_mpo_approx"], 3)
        for x, expected in ((20, 0), (30, .5), (40, 1), (90, 1), (105, .5), (120, 0)):
            self.assertAlmostEqual(task.window_desirability(x, 20, 40, 90, 120), expected)
        for boundary in (20, 40, 90, 120):
            left = task.window_desirability(boundary - 1e-8, 20, 40, 90, 120)
            right = task.window_desirability(boundary + 1e-8, 20, 40, 90, 120)
            self.assertLess(abs(left-right), 1e-7)

    def test_custom_logp_variant_is_distinct(self):
        result = task.cns_mpo(-1, 0, 300, 60, 0, None)
        self.assertEqual(result["cns_mpo_approx"], 6)
        self.assertEqual(result["custom_mpo_logp_2_4"], 5)

    def test_esol_coefficients_and_units(self):
        self.assertAlmostEqual(task.esol_log_s(3, 400, 4, .5), -4.316)
        self.assertAlmostEqual(task.molar_to_ug_ml(-3, 100), 100)
        delta = task.esol_log_s(3, 400, 5, .5) - task.esol_log_s(3, 400, 4, .5)
        self.assertAlmostEqual(delta, .066)

    def test_stability_and_invalid_inputs(self):
        self.assertEqual(task.sigmoid(1000), 1)
        self.assertEqual(task.sigmoid(-1000), 0)
        for value in (float("nan"), float("inf"), -float("inf")):
            with self.assertRaises(ValueError):
                task.low_desirability(value, 3, 5)
        with self.assertRaises(ValueError):
            task.low_desirability(2, 3, 3)


class ChemistryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.panel = task.load_panel()

    def test_panel_identity_and_frozen_copy(self):
        self.assertEqual(self.panel, json.loads((ROOT / "data/reference_panel.json").read_text(encoding="utf-8")))
        self.assertEqual(len(self.panel), 30)
        for row in self.panel:
            mol = task.parse_molecule(row)
            self.assertEqual(Chem.MolToInchiKey(mol), row["inchikey"])
            self.assertIn("pubchem.ncbi.nlm.nih.gov", row["structure_source"])

    def test_invalid_structures_and_identity_are_rejected(self):
        for smiles in ("not_a_smiles", "", "CC.O", "C[N+](C)(C)C"):
            with self.assertRaises(ValueError):
                task.parse_molecule({"name": "invalid", "smiles": smiles})
        row = dict(self.panel[0], inchikey="WRONG")
        with self.assertRaises(ValueError):
            task.parse_molecule(row)

    def test_duplicate_panel_is_rejected(self):
        records = copy.deepcopy(self.panel)
        records[1] = records[0]
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "panel.json"
            path.write_text(json.dumps(records), encoding="utf-8")
            with self.assertRaises(ValueError):
                task.load_panel(path)

    def test_ionization_direction_and_amide_exclusion(self):
        neutral = task.ionization(Chem.MolFromSmiles("CCC"), 2)
        base = task.ionization(Chem.MolFromSmiles("CCN"), 2)
        acid = task.ionization(Chem.MolFromSmiles("CC(=O)O"), 2)
        amide = task.ionization(Chem.MolFromSmiles("CC(=O)N"), 2)
        self.assertEqual(neutral["logd74_approx"], 2)
        self.assertLess(base["logd74_approx"], 2)
        self.assertLess(acid["logd74_approx"], 2)
        self.assertIsNone(amide["basic_pka_assumed"])
        override = task.ionization(Chem.MolFromSmiles("CCN"), 2,
                                  {"source": "synthetic test", "basic_pka": 7.4})
        self.assertAlmostEqual(override["neutral_fraction_approx"], .5)
        self.assertAlmostEqual(override["logd74_approx"], 2-math.log10(2))

    def test_override_validation(self):
        valid = {"Aspirin": {"source": "test", "basic_pka": None, "acidic_pka": 4}}
        self.assertEqual(task.validate_overrides(valid, {"Aspirin"}), valid)
        for value in ({"Aspirin": {}}, {"Wrong": {"source": "test"}},
                      {"Aspirin": {"source": "test", "logd74": None}},
                      {"Aspirin": {"source": "test", "basic_pka": float("nan")}},
                      {"Aspirin": {"source": "test", "typo": 5}}):
            with self.assertRaises(ValueError):
                task.validate_overrides(value, {"Aspirin"})

    def test_clinical_labels_cannot_leak_into_scores(self):
        row = self.panel[0]
        first = task.calculate_record(row, conformers=0)
        altered = task.calculate_record(dict(row, archetype="Toxic Dropouts", clinical_note="changed"), conformers=0)
        for column in ("cns_mpo_approx", "herg_liability_proxy", "oral_developability_score_proxy"):
            self.assertEqual(first[column], altered[column])
        self.assertIsNone(first["absolute_chameleonic_hb_index"])
        self.assertIsNone(first["chi_geometry_range_proxy"])
        self.assertEqual(first["geometry_status"], "disabled")

    def test_real_3d_geometry_reproducibility(self):
        first = task.geometry_calculation("CCO", 3, 42)
        second = task.geometry_calculation("CCO", 3, 42)
        self.assertEqual(first, second)
        self.assertGreater(first["conformers_converged"], 0)
        self.assertGreater(first["volume_3d_a3"], 0)
        self.assertEqual(first["imhb_vectors_mean"], 0)
        self.assertEqual(first["unpaired_donor_vectors_mean"], 1)
        self.assertIsNone(first["absolute_chameleonic_hb_index"])

    def test_isolated_worker_timeout_is_explicit(self):
        # Tiny time budget exercises termination without waiting for expensive
        # macrocycle embedding. The public CLI has a >=1 second lower bound.
        result = task.bounded_geometry(self.panel[20]["smiles"], 100, 42, .001)
        self.assertEqual(result["geometry_status"], "timeout")
        self.assertIsNone(result["volume_3d_a3"])

    def test_negative_embedding_sentinel_does_not_reach_mmff(self):
        with patch.object(task.AllChem, "EmbedMultipleConfs", return_value=[-1]), \
             patch.object(task.AllChem, "MMFFOptimizeMoleculeConfs") as optimize:
            result = task.geometry_calculation("CCO", 3, 42)
        self.assertEqual(result["geometry_status"], "embedding_failed")
        self.assertEqual(result["conformers_embedded"], 0)
        optimize.assert_not_called()


class EndToEndTests(unittest.TestCase):
    def test_offline_cli_and_publication_artifacts(self):
        with tempfile.TemporaryDirectory() as d:
            output = Path(d) / "output with spaces"
            run = subprocess.run([sys.executable, str(ROOT / "run_task1_mpo_admet_developability.py"),
                                  "--conformers", "0", "--output-dir", str(output)],
                                 capture_output=True, text=True, timeout=180)
            self.assertEqual(run.returncode, 0, run.stderr)
            records = json.loads((output / "results_task1/developability_results.json").read_text(encoding="utf-8"))
            self.assertEqual(len(records), 30)
            for row in records:
                self.assertTrue(0 <= row["cns_mpo_approx"] <= 6)
                self.assertTrue(0 <= row["herg_liability_proxy"] <= 1)
                self.assertTrue(0 <= row["hia_percent_proxy"] <= 100)
                self.assertTrue(0 <= row["oral_developability_score_proxy"] <= 100)
                self.assertGreater(row["solubility_ug_ml_esol"], 0)
                self.assertGreater(row["caco2_papp_proxy_1e_minus6_cm_s"], 0)
                self.assertIsNone(row["absolute_chameleonic_hb_index"])
                self.assertLessEqual(row["mpo_pka_sensitivity_min"], row["cns_mpo_approx"])
                self.assertGreaterEqual(row["mpo_pka_sensitivity_max"], row["cns_mpo_approx"])
                total = sum(row[k] for k in row if k.startswith("mpo_d_"))
                self.assertAlmostEqual(total, row["cns_mpo_approx"])
            manifest = json.loads((output / "results_task1/run_manifest.json").read_text(encoding="utf-8"))
            for name, digest in manifest["artifacts"].items():
                self.assertEqual(hashlib.sha256((output/name).read_bytes()).hexdigest(), digest)
            pngs = list((output / "figures_task1").glob("*.png"))
            self.assertEqual(len(pngs), 3)
            for path in pngs:
                with Image.open(path) as image:
                    self.assertGreaterEqual(image.width, 3000)
                    self.assertTrue(all(abs(dpi-300)<.1 for dpi in image.info["dpi"]))
                    self.assertGreater(np.asarray(image).std(), 10)

    def test_bad_cli_does_not_succeed(self):
        for arguments in (("--conformers", "-1"), ("--conformer-timeout", "nan"),
                          ("--strict-3d", "--conformers", "0")):
            with self.subTest(arguments=arguments):
                run = subprocess.run([sys.executable, str(ROOT / "run_task1_mpo_admet_developability.py"), *arguments],
                                     capture_output=True, timeout=30)
                self.assertNotEqual(run.returncode, 0)


if __name__ == "__main__":
    unittest.main()
