"""Conservation, equilibrium and boundary tests for Task 6."""
import copy
import importlib.util
import math
from pathlib import Path
import unittest

import numpy as np

SCRIPT = Path(__file__).resolve().parents[1]/"projects"/"task06_asd"/"run_task6_asd_formulation_supersaturation_kinetics.py"
spec = importlib.util.spec_from_file_location("task6_model", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class Task6PhysicsTests(unittest.TestCase):
    def test_mass_volume_roundtrip_and_density_effect(self):
        w=np.array([0,.1,.5,.9,1])
        phi=m.mass_to_volume_fraction(w,1.3,1.1)
        np.testing.assert_allclose(m.volume_to_mass_fraction(phi,1.3,1.1),w,atol=1e-15)
        self.assertLess(phi[2],w[2])

    def test_chi_unit_conversion(self):
        drug={"molar_volume_cm3_mol":100,"hansen_MPa_half":[2,0,0]}
        p={"hansen_MPa_half":[0,0,0]}
        self.assertAlmostEqual(m.chi_hansen(drug,p,300),100e-6*4e6/(m.R*300))

    def test_subcritical_has_no_phase_boundary(self):
        self.assertFalse(m.phase_boundaries(.6,100)["two_phase"])
        self.assertFalse(m.phase_boundaries(.605,100)["two_phase"])

    def test_critical_composition_orientation(self):
        p=m.phase_boundaries(.606,100)
        self.assertGreater(p["spinodal"][0],.8)
        self.assertLess(p["spinodal"][1],.97)

    def test_binodal_matches_common_tangent_and_convex_envelope(self):
        ch=.75; n=100
        b=m.phase_boundaries(ch,n)
        low,high=b["binodal"]
        f1,f2=m.free_energy(np.array([low,high]),ch,n)
        slope=(f2-f1)/(high-low)
        def derivative(x):
            return math.log(x)+1-(math.log1p(-x)+1)/n+ch*(1-2*x)
        self.assertAlmostEqual(derivative(low),slope,places=7)
        self.assertAlmostEqual(derivative(high),slope,places=7)
        x=np.linspace(low,high,501)
        self.assertGreaterEqual(float(np.min(m.free_energy(x,ch,n)-(f1+slope*(x-low)))),-1e-10)
        for x in b["spinodal"]:
            self.assertAlmostEqual(1/x+1/(n*(1-x))-2*ch,0,places=10)

    def test_GT_kelvin_endpoints_and_rh_sorption(self):
        self.assertEqual(m.gordon_taylor(0,300,400,2),400)
        self.assertEqual(m.gordon_taylor(1,300,400,2),300)
        d,p,s=m.DEFAULT["drug"],m.DEFAULT["polymers"][0],m.DEFAULT["storage"]
        dry=m.storage_result(d,p,.2,298.15,0,s)
        wet=m.storage_result(d,p,.2,298.15,.75,s)
        self.assertLess(wet["wet_Tg_K"],dry["wet_Tg_K"])
        self.assertLess(wet["assumed_water_mass_fraction"],.05)

    def test_vft_pole_is_undefined_not_shelf_life(self):
        row=m.storage_result(m.DEFAULT["drug"],m.DEFAULT["polymers"][3],.1,298.15,.60,m.DEFAULT["storage"])
        self.assertFalse(row["VFT_domain_valid"])
        self.assertIsNone(row["log10_induction_h_hypothetical"])

    def test_cnt_dimensions_and_limit(self):
        v=300e-6/m.NA; gamma=.004; temp=310.15
        b=m.cnt_barrier(2,gamma,v,temp)
        barrier_J=16*math.pi*gamma**3*v*v/(3*(m.KB*temp*math.log(2))**2)
        self.assertAlmostEqual(b,barrier_J/(m.KB*temp))
        self.assertTrue(math.isinf(m.cnt_barrier(1,gamma,v,temp)))
        self.assertAlmostEqual(m.cnt_barrier(2,2*gamma,v,temp),8*b)

    def test_closed_vessel_conserves_finite_dose(self):
        for form in ("crystalline","amorphous","ASD"):
            r=m.simulation(m.DEFAULT,form)
            np.testing.assert_allclose(r["states"][:4].sum(axis=0),100,atol=1e-7)
            self.assertGreaterEqual(float(r["states"][:4].min()),-1e-7)
            self.assertLessEqual(float(r["states"][1].max()),100+1e-7)
            self.assertEqual(float(r["states"][3,-1]),0)

    def test_crystal_cannot_supersaturate_and_absorption_monotone(self):
        r=m.simulation(m.DEFAULT,"crystalline",absorb=True)
        self.assertLessEqual(float(r["concentration"].max()),.01+1e-9)
        self.assertGreaterEqual(float(np.diff(r["states"][3]).min()),-1e-8)
        self.assertGreater(float(r["states"][3,-1]),0)

    def test_volume_factor_dilutes_initial_concentration_rate(self):
        cfg=copy.deepcopy(m.DEFAULT)
        cfg["kinetics"].update(hours=1e-5,samples=2)
        a=m.simulation(cfg,"amorphous")
        cfg["kinetics"]["volume_ml"]*=2
        b=m.simulation(cfg,"amorphous")
        self.assertAlmostEqual(a["concentration"][-1]/b["concentration"][-1],2,places=4)

    def test_absorption_proxy_is_integral_and_loading_changes_activity(self):
        r=m.simulation(m.DEFAULT,"ASD",absorb=True)
        expected=m.DEFAULT["kinetics"]["absorption_h_inv"]*900*r["summary"]["AUC_mg_h_ml"]
        self.assertAlmostEqual(expected,r["summary"]["absorbed_mg"],delta=.002)
        dilute=m.simulation(m.DEFAULT,"ASD",loading=.1)
        concentrated=m.simulation(m.DEFAULT,"ASD",loading=.5)
        self.assertLess(dilute["summary"]["source_solubility_mg_ml"],concentrated["summary"]["source_solubility_mg_ml"])

    def test_invalid_inputs_rejected(self):
        with self.assertRaises(ValueError):
            m.mass_to_volume_fraction(1.1,1,1)
        with self.assertRaises(ValueError):
            m.cnt_barrier(2,-.1,1,300)
        with self.assertRaises(ValueError):
            m.simulation(m.DEFAULT,"ASD",loading=0)

    def test_configuration_rejects_ignored_segments_and_unphysical_values(self):
        self.assertEqual(m.validate_config(copy.deepcopy(m.DEFAULT)),m.DEFAULT)
        for section,key,value in [("drug","N_segments",2),("kinetics","volume_ml",-10),("drug","molar_volume_cm3_mol",1),("storage","conditions",[[298.15,60],[313.15,.75]])]:
            cfg=copy.deepcopy(m.DEFAULT)
            cfg[section][key]=value
            with self.assertRaises(ValueError):
                m.validate_config(cfg)


if __name__ == "__main__":
    unittest.main()
