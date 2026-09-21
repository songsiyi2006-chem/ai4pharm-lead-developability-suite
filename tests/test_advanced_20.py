"""Chemical validity, objective ordering and rigid geometry invariants."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np
from rdkit import Chem
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('denovo',ROOT/'projects/task20_denovo/driver.py')
d=importlib.util.module_from_spec(spec);spec.loader.exec_module(d)

class DeNovoChecks(unittest.TestCase):
    def test_all_fragment_graphs_sanitize(self):
        for a in range(len(d.CORES)):
            for b in range(len(d.LEFT)):
                for c in range(len(d.RIGHT)):
                    mol=d.assemble((a,b,c))
                    self.assertFalse(any(x.GetAtomicNum()==0 for x in mol.GetAtoms()))
                    self.assertEqual(len(Chem.GetMolFrags(mol)),1)
    def test_invalid_anchor_is_rejected(self):
        with self.assertRaises(ValueError):d.join_fragment(Chem.MolFromSmiles('CC'),Chem.MolFromSmiles('[1*]O'),1)
    def test_rigid_transform_preserves_geometry_score(self):
        mol=Chem.MolFromSmiles('CO');xyz=np.array([[0.,0,0],[1.4,0,0]])
        p={'atoms':np.array([[0.,4,0],[2.,4,0]]),'radii':np.array([1.7,1.52]),
           'acceptor_targets':np.array([[0.,3,0]]),'donor_targets':np.array([[2.,3,0]])}
        rotation=d.Rotation.from_euler('xyz',[.4,.7,.2]).as_matrix();offset=np.array([12.,-3,8])
        transformed={k:(v@rotation.T+offset if k!='radii' else v) for k,v in p.items()}
        before=d.pose_score(xyz,mol,p);after=d.pose_score(xyz@rotation.T+offset,mol,transformed)
        np.testing.assert_allclose(before,after,atol=1e-12)
    def test_steric_clash_penalizes_exact_overlap(self):
        mol=Chem.MolFromSmiles('C');p={'atoms':np.array([[0.,0,0]]),'radii':np.array([1.7]),
                'acceptor_targets':np.zeros((1,3)),'donor_targets':np.zeros((1,3))}
        clash=d.pose_score(np.zeros((1,3)),mol,p)[0][0]
        contact=d.pose_score(np.array([[4.,0,0]]),mol,p)[0][0]
        self.assertLess(clash,contact)
    def test_nondominated_sort_preserves_tradeoff(self):
        def row(i,g,m,q,s):return {'smiles':str(i),'pocket_geometry_score':g,'cns_mpo_proxy':m,'qed':q,'sa_score':s}
        records=[row(0,1,1,.5,3),row(1,2,2,.7,2),row(2,3,.5,.4,4)]
        self.assertEqual(d.pareto_fronts(records),[[1,2],[0]])
        self.assertEqual({r['smiles'] for r in d.select_population(records,2)},{'1','2'})
    def test_pocket_landmarks_have_angstrom_scale(self):
        p=d.pocket(ROOT/'projects/task20_denovo/inputs/6YB7.pdb')
        self.assertGreater(p['protein_atoms'],2000);self.assertGreater(p['pocket_atoms'],100)
        self.assertTrue(np.all(np.linalg.norm(p['atoms']-p['center'],axis=1)<13))
    def test_sascore_and_mpo_finite_on_druglike_graph(self):
        x=d.descriptors(Chem.MolFromSmiles('CC(=O)Nc1ccc(O)cc1'))
        self.assertTrue(1<=x['sa_score']<=10);self.assertTrue(0<=x['cns_mpo_proxy']<=6)

if __name__=='__main__':unittest.main()
