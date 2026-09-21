"""Independent invariants and controls for the publishable Tasks 12–14.

Task 11 is intentionally withheld from this release; its implementation and
tests are not part of the published tree.
"""
import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
def load(task):
    spec=importlib.util.spec_from_file_location(task,ROOT/'projects'/task/'driver.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
CHROM=load('task12_chromatin');BBB=load('task13_bbb_tfr');BITE=load('task14_bite')

class ChromatinTests(unittest.TestCase):
    def test_binding_mass_balances(self):
        for dose in [.001,1,100,1e4]:
            y=CHROM.equilibrium(dose)
            self.assertAlmostEqual(y[[2,4,5,6,7,8]].sum(),dose,places=7)
            self.assertAlmostEqual(y[[3,6,7,8]].sum(),20,places=7)
            self.assertAlmostEqual(y[[0,4,7]].sum(),30,places=7)
            self.assertAlmostEqual(y[[1,5,8]].sum(),70,places=7)
    def test_thermodynamic_cycle(self):
        y=CHROM.equilibrium(100,penalty=2)
        self.assertAlmostEqual(y[7]/y[8],10/(10*np.exp(-2/CHROM.CONFIG['RT_kcal_per_mol']))*y[0]/y[1],places=8)
    def test_no_drug_matches_analytic_pool_relaxation(self):
        s=CHROM.simulate(0)
        rate=.14+.06+.025;steady=.14*100/rate
        expected=steady+(70-steady)*np.exp(-rate*s.t)
        np.testing.assert_allclose(s.y[1],expected,rtol=1e-7,atol=2e-6)
    def test_stoichiometric_conservation_under_arbitrary_state(self):
        y=np.arange(1.,12.);dy=CHROM.rhs(0,y)
        self.assertAlmostEqual(sum(dy[[2,4,5,6,7,8]]),0,places=9)
        self.assertAlmostEqual(sum(dy[[3,6,7,8]]),0,places=9)
        self.assertAlmostEqual(sum(dy[[0,1,4,5,7,8]])+dy[9]-dy[10],0,places=9)

class BBBTests(unittest.TestCase):
    def test_antibody_and_receptor_ledgers(self):
        self.assertLess(max(BBB.balances(BBB.simulate())),1e-7)
    def test_no_receptor_matches_analytic_plasma_clearance(self):
        s=BBB.simulate(receptor=0.,rtol=2e-10)
        np.testing.assert_allclose(s.y[0],300*np.exp(-.03*s.t),rtol=1e-8)
        self.assertEqual(float(s.y[6].max()),0.)
    def test_no_internalization_blocks_brain_entry(self):
        s=BBB.simulate(internalization=0.);self.assertEqual(float(s.y[6].max()),0.)
    def test_vascular_contamination_is_exact_volume_ratio(self):
        _,_,kp,contaminated=BBB.exposure(BBB.simulate())
        self.assertAlmostEqual(contaminated-kp,.015/.3,places=12)

class BiTETests(unittest.TestCase):
    def test_mass_balance_and_terminal_symmetry(self):
        y=BITE.equilibrium(25,8.,16.,100.,1.,3.)
        self.assertAlmostEqual(y[[0,3,5]].sum(),8.,places=8)
        self.assertAlmostEqual(y[[1,4,5]].sum(),16.,places=8)
        self.assertAlmostEqual(y[[2,3,4,5]].sum(),25.,places=8)
        z=BITE.equilibrium(25,16.,8.,1.,100.,3.)
        self.assertAlmostEqual(y[5],z[5],places=8)
    def test_binary_only_matches_quadratic_solution(self):
        p,e,k=20.,8.,3.
        expected=2*p*e/(p+e+k+np.sqrt((p+e+k)**2-4*p*e))
        y=BITE.equilibrium(p,e,0.,k,1.)
        self.assertAlmostEqual(y[3],expected,places=9);self.assertEqual(y[5],0)
    def test_molecule_conversion_cannot_exceed_available_receptors(self):
        for dose in [1.,10.,100.,1e5]:
            self.assertLessEqual(BITE.synapses(dose),5000.+1e-8)
    def test_cell_loss_ledger_and_missing_effector_control(self):
        s=BITE.simulate();np.testing.assert_allclose(s.y.sum(axis=0),1,atol=1e-12)
        self.assertEqual(BITE.simulate(ratio=0.).y[1,-1],0.)

class OutputSafetyTests(unittest.TestCase):
    def test_existing_directory_is_rejected_before_any_write(self):
        with tempfile.TemporaryDirectory() as d:
            sentinel=Path(d)/'keep.txt';sentinel.write_text('preserve',encoding='utf-8')
            for module in [CHROM,BBB,BITE]:
                with self.assertRaises(FileExistsError):module.run(Path(d))
            self.assertEqual(sentinel.read_text(encoding='utf-8'),'preserve')

if __name__=='__main__':unittest.main()
