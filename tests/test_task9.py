"""Analytic invariants for Task 9; no scientific validation implied."""
import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch
import contextlib
import io
import tempfile
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/'projects/task09_cryoem_allostery/run_task9_cryoem_cryptic_pocket_allostery.py'
spec=importlib.util.spec_from_file_location('task9_under_test',SCRIPT)
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

class TestTask9(unittest.TestCase):
    def test_generator_detailed_balance(self):
        t,p=m.generator(m.DEFAULT)
        np.testing.assert_allclose(t.sum(1),1)
        np.testing.assert_allclose(p[:,None]*t,p[None,:]*t.T,atol=1e-14)
        np.testing.assert_allclose(p@t,p,atol=1e-14)

    def test_estimator_balance_and_observed_counts(self):
        states=np.array([0,1,2,3,4,4,3,2,1,0])
        t,p,c=m.estimate_msm(states)
        self.assertEqual(c.sum(),9)
        np.testing.assert_allclose(p[:,None]*t,p[None,:]*t.T,atol=1e-14)

    def test_mfpt_two_state_analytic(self):
        t=np.array([[.8,.2],[.1,.9]])
        self.assertAlmostEqual(m.mfpt(t,1,10)[0],50)

    def test_mfpt_bellman_equation(self):
        t,_=m.generator(m.DEFAULT); a=m.mfpt(t,4,10)
        np.testing.assert_allclose(a[:4],10+t[:4]@a)
        self.assertEqual(a[4],0)

    def test_sasa_isolated_sphere(self):
        area=m.lining_sasa(np.zeros((1,3)),np.array([1.7]),[0],1.4,96)
        self.assertAlmostEqual(area,4*np.pi*3.1**2,places=9)

    def test_surface_distance_accounts_for_radii(self):
        x=np.array([[0.,0,0],[3.,0,0]])
        pts=np.array([[1.,0,0],[2.5,0,0]])
        radii=np.array([1.,2.])
        got=m.surface_clearance(m.spatial.cKDTree(x),pts,radii)
        expected=(m.spatial.distance.cdist(pts,x)-radii).min(1)
        np.testing.assert_allclose(got,expected)

    def test_anm_tetrahedron_pseudoinverse_and_prs(self):
        x=np.array([[0.,0,0],[1.,0,0],[0.,1,0],[0.,0,1]])
        p=m.prs_calculation(x,3,[64,256,4096],123)
        self.assertEqual(p['zero_modes'],6)
        self.assertLess(p['pseudoinverse_residual'],1e-12)
        self.assertLess(p['convergence'][-1]['relative_frobenius_error'],.05)
        np.testing.assert_allclose(p['hessian'],p['hessian'].T)
        np.testing.assert_allclose(np.diag(p['normalized']),1)

    def test_hessian_translation_zero(self):
        x=np.array([[0.,0,0],[1.,0,0],[0.,1,0],[0.,0,1]])
        h,_=m.anm_hessian(x,3)
        for direction in np.eye(3): np.testing.assert_allclose(h@np.tile(direction,4),0,atol=1e-14)

    def test_pdb_matching_and_ligand_excluded(self):
        keys,x,y,ca,res,elements,center=m.align_anchors()
        self.assertGreater(len(ca),150)
        self.assertEqual(x.shape,y.shape)
        self.assertNotIn('3GZ',res)
        self.assertTrue(np.isfinite(center).all())
        self.assertEqual(len(set(keys)),len(keys))

    def test_synthetic_state_reproducibility(self):
        t,p=m.generator(m.DEFAULT)
        a=m.simulate_states(t,100,np.random.default_rng(42),p)
        b=m.simulate_states(t,100,np.random.default_rng(42),p)
        np.testing.assert_array_equal(a,b)

    def test_config_validation(self):
        m.validate_config(m.DEFAULT)
        for k,v in [('lag_ns',0),('n_conformers',10),('grid_spacing_A',.01),('prs_force_counts',[1,2,3])]:
            c=dict(m.DEFAULT); c[k]=v
            with self.assertRaises(ValueError): m.validate_config(c)

    def test_empty_cavity_is_finite_and_zero(self):
        c=dict(m.DEFAULT); c['roi_halfwidth_A']=1.; c['grid_spacing_A']=1.
        p,lining=m.pocket_metrics(np.zeros((1,3)),np.array(['C']),['ALA'],np.zeros(3),c)
        self.assertEqual(p['volume_A3'],0)
        self.assertEqual(p['bounded_uncalibrated_score'],0)
        self.assertEqual(len(lining),0)

    def test_cli_refuses_archived_output_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp, patch('sys.argv',['task9','--out',tmp]), contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as caught: m.main()
            self.assertEqual(caught.exception.code,2)

    def test_boundary_cavity_has_unknown_full_volume(self):
        keys,x,y,ca,res,elements,center=m.align_anchors()
        p,_=m.pocket_metrics(y,elements,res,center,m.DEFAULT,with_sasa=False)
        self.assertTrue(p['roi_boundary_contact'])
        self.assertIsNone(p['full_component_volume_A3'])
        self.assertEqual(p['threshold_classification'],'indeterminate_roi_truncated')

if __name__=='__main__': unittest.main()
