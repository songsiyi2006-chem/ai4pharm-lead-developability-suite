"""ADC numerical tests target physical limits and independent conservation laws."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import unittest
import numpy as np
from PIL import Image
from scipy.sparse.linalg import expm_multiply

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from projects.task08_adc import run_task8_adc_dar_cleavage_bystander_dynamics as adc


class ADCTests(unittest.TestCase):
    def test_analytic_and_conservation_invariants(self):
        checks=adc.numeric_checks(adc.DEFAULT)
        self.assertGreaterEqual(len(checks),13)
        self.assertTrue(all(x["passed"] for x in checks),checks)

    def test_invalid_configuration(self):
        for key,value in [("radial_cells",2.5),("extracellular_fraction",1),("core_radius_um",201),("time_step_h",0),("kon_per_nm_h",float("nan")),("ssa_replicates",1),("downstream_payload_equivalents",3),("output_interval_h",0.37)]:
            with self.subTest(key=key),self.assertRaises(ValueError):
                adc.validate_config(dict(adc.DEFAULT,**{key:value}))
        with self.assertRaises(ValueError): adc.validate_config(dict(adc.DEFAULT,unexpected=1))

    def test_reduction_probability_and_capacity(self):
        c=dict(adc.DEFAULT,reduction_rate_h=np.log(2),reduction_time_h=1)
        y=adc.initial_conjugation(c,5)
        capacities=np.array([m for m,n in adc.STATES])
        self.assertAlmostEqual(float(y[:-1].sum()),1)
        self.assertAlmostEqual(float(y[:-1]@capacities),4)
        self.assertEqual(y[-1],5)

    def test_finite_linker_bounds_average_dar(self):
        c=dict(adc.DEFAULT,conjugation_time_h=30)
        for equivalents in (0.1,2,12):
            r=adc.conjugate(c,equivalents)
            self.assertLessEqual(r["mean"][-1],min(equivalents,8*adc.reduction_prob(c))+1e-8)
            self.assertGreaterEqual(r["free_linker"].min(),-1e-9)

    def test_gillespie_reproducible_and_stoichiometric(self):
        c=dict(adc.DEFAULT,ssa_replicates=2,ssa_antibodies=100)
        a=adc.gillespie(c,2);b=adc.gillespie(c,2)
        np.testing.assert_array_equal(a,b)
        np.testing.assert_allclose(a.sum(axis=1),1)
        self.assertTrue(np.all(a@np.arange(9)<=2+1e-12))

    def test_hic_all_nine_species_and_area_recovery(self):
        fractions=np.arange(1,10,dtype=float);fractions/=fractions.sum()
        h=adc.hic(adc.DEFAULT,fractions)
        self.assertEqual(h["components"].shape[1],9)
        np.testing.assert_allclose(np.trapezoid(h["components"],h["time"],axis=0),fractions,atol=1e-12)
        np.testing.assert_allclose(h["fitted_fractions"],fractions,atol=1e-12)
        self.assertTrue(0<h["dar4_purity"]<1)

    def test_no_binding_no_cellular_payload(self):
        c=dict(adc.DEFAULT,kon_per_nm_h=0,horizon_h=2)
        r=adc.cell_model(c,4,.2)
        np.testing.assert_allclose(r["values"][[2,3,4,6,7,8,9]],0,atol=1e-9)

    def test_finite_volume_closed_uniform_decay(self):
        c=dict(adc.DEFAULT,diffusion_um2_h=120,extracellular_clearance_h=.1,payload_metabolism_h=.1)
        grid=adc.radial_operator(c,.3,20,absorbing=False)
        initial=np.r_[grid["ve"],grid["vi"]]
        answer=expm_multiply(grid["matrix"]*3,initial)
        np.testing.assert_allclose(answer,initial*np.exp(-.3),rtol=2e-12,atol=1e-8)

    def test_radial_operator_metzler_and_boundary_loss(self):
        g=adc.radial_operator(adc.DEFAULT,.2,20)
        a=g["matrix"].toarray();off=a.copy();np.fill_diagonal(off,0)
        self.assertGreaterEqual(off.min(),0)
        sums=a.sum(axis=0)
        np.testing.assert_allclose(sums[:19],-adc.DEFAULT["extracellular_clearance_h"],atol=1e-12)
        self.assertAlmostEqual(sums[19],-adc.DEFAULT["extracellular_clearance_h"]-g["boundary_rate"],places=11)

    def test_zero_permeability_no_bystander_and_custom_horizon(self):
        c=dict(adc.DEFAULT,horizon_h=2,radial_cells=20)
        source=adc.cell_model(c,4,0)
        ti=adc.tissue_model(c,source,0)
        self.assertEqual(ti["horizon_h"],2)
        self.assertEqual(ti["bystander_survival_at_horizon"],1)
        self.assertEqual(ti["ce"].max(),0)
        self.assertIsNone(ti["threshold_radius_um"])

    def test_archived_evidence_and_300dpi_figures(self):
        base=ROOT/"projects/task08_adc"
        summary=json.loads((base/"summary.json").read_text(encoding="utf-8"))
        self.assertTrue(summary["verification_passed"])
        self.assertEqual(summary["computation_counts"]["radial_data_rows"],23200)
        for file in (base/"figures_task8").glob("*.png"):
            with Image.open(file) as im:
                self.assertGreaterEqual(min(im.info["dpi"]),299)
        self.assertEqual(len(list((base/"figures_task8").glob("*.png"))),4)

    def test_project_manifest(self):
        base=ROOT/"projects/task08_adc"
        manifest=json.loads((base/"manifest.json").read_text(encoding="utf-8"))
        for name,digest in manifest["files_sha256"].items():
            self.assertEqual(hashlib.sha256((base/name).read_bytes()).hexdigest(),digest,name)
        self.assertEqual(hashlib.sha256(Path(adc.__file__).read_bytes()).hexdigest(),manifest["script_sha256"])


if __name__=="__main__": unittest.main()
