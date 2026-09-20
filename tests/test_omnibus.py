"""Integration safety and plotting-data invariants, without expensive reruns."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import run_ai4pharm_omnibus_suite as suite
from omnibus.figures_1_5 import grid, vals
from omnibus.catalog import PROJECTS, DRIVERS


class OmnibusTests(unittest.TestCase):
    def test_output_protects_repository_ancestors_and_artifacts(self):
        for path in [ROOT, ROOT.parent, ROOT/"projects/new", ROOT/"docs/new", ROOT/"data_omnibus/new"]:
            with self.subTest(path=path), self.assertRaises(ValueError):
                suite.prepare_output(path)

    def test_output_refuses_existing_even_empty_directory(self):
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaises(FileExistsError):
                suite.prepare_output(Path(folder))

    def test_new_output_preserves_siblings(self):
        with tempfile.TemporaryDirectory() as folder:
            base = Path(folder)
            marker = base/"evidence.txt"
            marker.write_text("unchanged")
            created = suite.prepare_output(base/"new")
            self.assertTrue(created.is_dir())
            self.assertEqual(marker.read_text(), "unchanged")

    def test_routing_uses_original_drivers_without_git_flags(self):
        for task in range(1, 11):
            cmd, target = suite.child_command(task, ROOT/"work/example", python="python")
            self.assertEqual(Path(cmd[1]), ROOT/"projects"/PROJECTS[task-1]/DRIVERS[task-1])
            self.assertNotIn("--git-sync", cmd)
            self.assertNotIn("--git-publish", cmd)
            self.assertIn("--output-dir" if task <= 3 else "--out", cmd)
            self.assertTrue(target.is_relative_to(ROOT/"work/example"))
            self.assertEqual(target.name, "results" if task in (5, 10) else PROJECTS[task-1])

    def test_vienna_library_is_explicit_and_task10_only(self):
        for task in (1, 9, 10):
            cmd, _ = suite.child_command(task, ROOT/"work/example", rna_library=ROOT/"work/library")
            self.assertEqual("--rna-library" in cmd, task == 10)

    def test_failed_child_stops_chain_and_keeps_record(self):
        with tempfile.TemporaryDirectory() as folder:
            with patch.object(suite, "run_child", return_value=7) as run:
                with self.assertRaisesRegex(RuntimeError, "Task 1 failed"):
                    suite.recompute(Path(folder), 60)
            self.assertEqual(run.call_count, 1)
            record = json.loads((Path(folder)/"data_omnibus/execution.json").read_text())[0]
            self.assertEqual(record["returncode"], 7)
            self.assertEqual(record["status"], "failed")
            self.assertTrue((Path(folder)/record["log"]).is_file())

    def test_timeout_stops_chain_without_success_status(self):
        with tempfile.TemporaryDirectory() as folder:
            with patch.object(suite, "run_child", side_effect=subprocess.TimeoutExpired("python", 1)):
                with self.assertRaisesRegex(RuntimeError, "Task 1 timeout"):
                    suite.recompute(Path(folder), 1)
            record = json.loads((Path(folder)/"data_omnibus/execution.json").read_text())[0]
            self.assertEqual(record["status"], "timeout")
            self.assertIsNone(record["returncode"])

    def test_timeout_terminates_descendants(self):
        with patch.object(suite.subprocess, "Popen") as popen, patch.object(suite, "stop_process_tree") as stop:
            process = popen.return_value.__enter__.return_value
            process.wait.side_effect = subprocess.TimeoutExpired("python", 1)
            stop.return_value = {"method": "test tree termination"}
            with self.assertRaises(subprocess.TimeoutExpired) as caught:
                suite.run_child(["python"], None, 1)
            stop.assert_called_once_with(process)
            self.assertEqual(caught.exception.tree_termination["method"], "test tree termination")

    def test_missing_grid_cells_never_become_zero(self):
        data = [{"x": 1, "y": 1, "z": 2}, {"x": 2, "y": 2, "z": 3}]
        with self.assertRaises(ValueError):
            grid(data, "x", "y", "z")

    def test_duplicate_grid_cells_are_rejected(self):
        with self.assertRaises(ValueError):
            grid([{"x": 1, "y": 1, "z": 2}]*2, "x", "y", "z")

    def test_missing_numerical_plot_values_are_rejected(self):
        for value in ("nan", "inf", ""):
            with self.subTest(value=value), self.assertRaises(ValueError):
                vals([{"x": value}], "x")

    def test_input_hash_audit_detects_corruption_and_path_escape(self):
        spec = importlib.util.spec_from_file_location("omnibus_validator", ROOT/"tools/validate_omnibus.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as folder:
            base = Path(folder)
            data = base/"input.csv"
            data.write_text("x\n1\n")
            expected = module.digest(data)
            data.write_text("x\n2\n")
            errors = []
            module.verify_map(base, {"input.csv": expected, "../escape.csv": expected}, "input", errors)
            self.assertEqual(len(errors), 2)
            self.assertIn("checksum mismatch", errors[0])
            self.assertIn("escaping path", errors[1])


if __name__ == "__main__":
    unittest.main()
