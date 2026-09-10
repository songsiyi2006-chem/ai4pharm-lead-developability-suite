"""Run the PBPK invariants and check versioned artifacts in the suite CI."""
import hashlib
import json
from pathlib import Path
import sys
import unittest

from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from task4_pbpk import run_task4_pbpk_pharmacokinetics_dose_prediction as pbpk


class PBPKTests(unittest.TestCase):
    def test_numerical_invariants(self):
        checks = pbpk.self_tests()
        self.assertGreaterEqual(len(checks), 25)
        for check in checks:
            with self.subTest(check=check['test']):
                self.assertTrue(check['passed'])

    def test_delivered_artifacts_match_manifest(self):
        folder = ROOT / 'task4_pbpk'
        manifest = json.loads((folder / 'manifest.json').read_text(encoding='utf-8'))
        for relative, digest in manifest['files_sha256'].items():
            with self.subTest(file=relative):
                path = folder / relative
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), digest)
        script = folder / 'run_task4_pbpk_pharmacokinetics_dose_prediction.py'
        self.assertEqual(hashlib.sha256(script.read_bytes()).hexdigest(), manifest['script_sha256'])
        figures = list((folder / 'figures_task4').glob('*.png'))
        self.assertEqual(len(figures), 4)
        for path in figures:
            with Image.open(path) as image:
                self.assertTrue(all(abs(dpi - 300) < 0.02 for dpi in image.info['dpi']))


if __name__ == '__main__':
    unittest.main()
