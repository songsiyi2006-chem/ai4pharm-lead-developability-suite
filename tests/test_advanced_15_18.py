"""Independent physical invariants and phenotype boundaries, Tasks 15, 16 and 18."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
def load(name,folder):
    spec=importlib.util.spec_from_file_location(name,ROOT/'projects'/folder/'driver.py')
    m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
a=load('test_adv15','task15_receptor_signaling')
b=load('test_adv16','task16_transporters')
d=load('test_adv18','task18_pgx')

class ProofreadingPhysicalTests(unittest.TestCase):
    def test_receptor_flux_conservation_arbitrary_states(self):
        y=np.arange(1.,14.);y/=y.sum()
        self.assertAlmostEqual(a.receptor_rhs(y,[73.,14.]).sum(),0.,places=11)

    def test_equal_antigens_are_indistinguishable(self):
        free,blocks=a.steady_proofreading(koff=[36.,36.])
        np.testing.assert_allclose(blocks[0],blocks[1],rtol=1e-14)
        self.assertAlmostEqual(free+blocks.sum(),1.,places=14)

    def test_terminal_conditional_probability(self):
        _,blocks=a.steady_proofreading()
        for i,off in enumerate([36.,360.]):
            self.assertAlmostEqual(blocks[i,-1]/blocks[i].sum(),(720/(720+off))**5,places=14)

    def test_ligand_clearance_and_signaling_are_distinct(self):
        y=np.zeros(19);y[0]=1.;y[15]=2.
        # At the same instantaneous cytokine state, blocking the receptor can
        # reduce ligand removal; it must not be modeled as forced IL6 depletion.
        self.assertGreater(a.derivative(12.,y,100.)[15],a.derivative(12.,y,0.)[15])
        self.assertGreater(a.blockade(12.,100.),.99-1e-4)

class TransportConservationTests(unittest.TestCase):
    def test_competitive_inhibition_changes_km_not_vmax(self):
        self.assertLess(b.michaelis(.1,16.,.1,100.),b.michaelis(.1,16.,.1))
        self.assertAlmostEqual(b.michaelis(1e12,16.,.1,100.),16.,places=8)

    def test_augmented_amount_ledger(self):
        m=b.TransportPBPK(10.);y=np.arange(1.,15.)
        dy=m.derivative(0.,y)
        self.assertAlmostEqual(b.mass_ledger(dy),0.,places=9)

    def test_nonnegative_boundary_inward_flux(self):
        m=b.TransportPBPK(0.)
        for idx in list(range(10))+[11,12,13]:
            y=np.ones(14);y[idx]=0.
            self.assertGreaterEqual(m.derivative(0.,y)[idx],-1e-10)

    def test_uses_existing_task4_class(self):
        self.assertIsInstance(b.TransportPBPK().base,b.task4.PBPK)
        self.assertEqual(b.TASK4_PATH.resolve(),(ROOT/'projects/task04_pbpk/run_task4_pbpk_pharmacokinetics_dose_prediction.py').resolve())

    def test_no_saturable_flux_from_empty_donor(self):
        self.assertEqual(b.michaelis(0.,1.,1.),0.)

class PharmacogenomicBoundaryTests(unittest.TestCase):
    def test_cyp2d6_current_cpic_boundaries(self):
        self.assertEqual([d.cyp2d6_phenotype(x) for x in [0.,.25,1.,1.25,2.,2.25,2.5]],
                         ['PM','IM','IM','NM','NM','NM','UM'])

    def test_phased_copy_number_and_unknown_alleles(self):
        self.assertEqual(d.cyp2d6_diplotype(['*1x2','*1'])['activity_score'],3.)
        self.assertEqual(d.cyp2d6_diplotype(['*10','*10'])['activity_score'],.5)
        self.assertEqual(d.cyp2d6_diplotype(['*999','*1'])['phenotype'],'indeterminate')

    def test_cyp2c19_uses_its_own_translation(self):
        self.assertEqual(d.cyp2c19_diplotype(['*1','*17']),'RM')
        self.assertEqual(d.cyp2c19_diplotype(['*2','*17']),'IM')
        self.assertEqual(d.cyp2c19_diplotype(['*2','*9']),'likely_PM')
        self.assertEqual(d.cyp2c19_diplotype(['*12','*1']),'indeterminate')

    def test_molar_equivalent_mass_conservation(self):
        for factor in [0.,.35,1.,1.8]:
            np.testing.assert_allclose(d.pk_matrix(factor)[:4].sum(axis=0),0.,atol=1e-15)

    def test_parent_and_prodrug_metabolite_have_opposite_direction(self):
        _,pm,_=d.simulate(0.);_,nm,_=d.simulate(1.)
        self.assertGreater(pm[4,-1],nm[4,-1]);self.assertEqual(pm[5,-1],0.);self.assertGreater(nm[5,-1],0.)

    def test_inhibition_does_not_modify_genetic_classification(self):
        classification=d.cyp2d6_diplotype(['*1','*1'])
        _,y,_=d.simulate(1.,100.)
        self.assertEqual(classification['phenotype'],'NM');self.assertGreater(y[4,-1],0.)

class OutputProtectionTests(unittest.TestCase):
    def test_existing_directory_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            for module in [a,b,d]:
                with self.assertRaises(FileExistsError):module.run(Path(tmp))

if __name__=='__main__':unittest.main()
