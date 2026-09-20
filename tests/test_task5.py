"""Mechanistic invariants and independent analytic limits for Task 5."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np
from scipy.integrate import solve_ivp

P=Path(__file__).resolve().parents[1]/'projects/task05_cyp_ddi/run_task5_cyp_ddi_mechanism_based_inhibition.py'
spec=importlib.util.spec_from_file_location('task5_test_target',P)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class Task5Tests(unittest.TestCase):
    def test_units(self):
        self.assertAlmostEqual(float(m.mg_l_to_uM(1.,500.)),2.)

    def test_bad_kinetics(self):
        d=m.panel()[0];d['ki_uM'][0]=0
        with self.assertRaises(ValueError): m.validate_compound(d)

    def test_no_inhibitor_limit_and_correct_direction(self):
        self.assertAlmostEqual(m.official_static_aucr(1.,1.),1.,places=14)
        self.assertGreater(m.official_static_aucr(.1,.2),m.official_static_aucr(.5,.6))

    def test_static_complete_blockade_limit(self):
        self.assertAlmostEqual(m.official_static_aucr(0.,0.,.9,.5),20.)

    def test_enzyme_constant_exposure_analytic(self):
        kd=.019;k=.2;ss=kd/(kd+k)
        t=np.linspace(0,36,501)
        exact=ss+(1-ss)*np.exp(-(kd+k)*t)
        sol=solve_ivp(lambda t,y:kd*(1-y)-k*y,(0,36),[1.],t_eval=t,rtol=1e-10,atol=1e-12)
        np.testing.assert_allclose(sol.y[0],exact,rtol=2e-9,atol=1e-10)

    def test_task4_uninhibited_rhs_identical(self):
        p=m.probes()[0];model=m.pbpk_model(p);y=np.arange(11,dtype=float)/11
        np.testing.assert_allclose(m.victim_rhs(model,p,y,np.ones(3),np.ones(3)),model.derivative(y),rtol=1e-13,atol=1e-13)

    def test_modified_rhs_mass_conserved(self):
        for p in m.probes():
            model=m.pbpk_model(p);y=np.arange(1,12,dtype=float)
            derivative=m.victim_rhs(model,p,y,np.array([.01,.2,.3]),np.array([.05,.4,.8]))
            self.assertAlmostEqual(float(derivative[:10].sum()),0.,places=9)

    def test_constant_pbpk_analytic_auc(self):
        p=m.probes()[0];model=m.pbpk_model(p);ah=np.array([.1,.3,.7]);ag=np.array([.2,.5,.8])
        y=np.zeros(11);y[m.t4.GUT]=p['dose_mg']
        s=solve_ivp(lambda t,y:m.victim_rhs(model,p,y,ah,ag),(0,1500),y,method='BDF',rtol=1e-9,atol=1e-12)
        expected=m.generalized_static(model,p,ah,ag)
        self.assertLess(s.y[:7,-1].sum()/p['dose_mg'],1e-8)
        self.assertAlmostEqual(s.y[m.t4.AUC,-1]/expected,1.,places=7)

    def test_recovery_analytic_and_censoring(self):
        class Mock:
            def __init__(self,kd): self.kd=kd
            def at(self,t):
                v=1-.8*np.exp(-self.kd*(np.asarray(t)-m.LAST_DOSE))
                return np.broadcast_to(v,(17,)+np.shape(v))
        expected=np.log(8)/.019/24
        self.assertAlmostEqual(m.recovery(Mock(.019),'hepatic',0),expected,places=8)
        self.assertIsNone(m.recovery(Mock(.0001),'hepatic',0))

    def test_classification_boundaries(self):
        self.assertEqual([m.severity(r) for r in [1.,1.25,2.,5.]],['below_weak_threshold','weak_range','moderate_range','strong_range'])

    def test_m12_factor_five_not_old_factor_fifty(self):
        r1,r2=m.basic_ratios(.1,5.,.3,1.,.02)
        self.assertAlmostEqual(r1,1.02)
        self.assertAlmostEqual(r2,6.)


if __name__=='__main__': unittest.main()
