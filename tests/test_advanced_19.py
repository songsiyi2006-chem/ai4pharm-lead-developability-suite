"""Physical and evidence-boundary tests; no long molecular dynamics in CI."""
import gzip
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
PROJECT=ROOT/'projects/task19_metadynamics'
spec=importlib.util.spec_from_file_location('advanced19',PROJECT/'driver.py')
t19=importlib.util.module_from_spec(spec);spec.loader.exec_module(t19)
try:
    import openmm as mm
    from openmm import unit
except ImportError:
    mm=None


class StructureAndEvidence(unittest.TestCase):
    def test_actual_structure_identity_and_internal_gap(self):
        atoms=t19.structure_atoms(PROJECT/'inputs/5T35.pdb')
        self.assertEqual(len([a for a in atoms if a['chain']=='D' and a['resname']=='759']),69)
        prep=json.loads((PROJECT/'inputs/prepared/preparation.json').read_text())
        self.assertEqual(len(prep['internal_missing_residues_modeled']['(2, 32)']),10)
        self.assertGreater(prep['particles'],100000)

    def test_orientation_rigid_motion_invariance(self):
        xyz=np.array([[1.,0.,0.],[0.,0.,0.],[0.,1.,0.],[.5,1.,1.]])
        theta=.47
        rot=np.array([[math.cos(theta),-math.sin(theta),0.],[math.sin(theta),math.cos(theta),0.],[0.,0.,1.]])
        self.assertAlmostEqual(t19.orientation_angle(xyz),t19.orientation_angle(xyz@rot.T+12.),places=12)

    def test_undefined_orientation_rejected(self):
        with self.assertRaises(ValueError):t19.orientation_angle([[0,0,0],[1,0,0],[2,0,0],[3,0,0]])

    def test_output_overwrite_guard(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'important.txt').write_text('keep')
            with self.assertRaises(FileExistsError):t19.run(p)
            self.assertEqual((p/'important.txt').read_text(),'keep')

    def test_cached_parameters_are_hash_verified(self):
        cache=PROJECT/'inputs/prepared';m=json.loads((cache/'manifest.json').read_text())
        for name,digest in m['files_sha256'].items():self.assertEqual(t19.sha(cache/name),digest,name)


@unittest.skipIf(mm is None,'Optional OpenMM engine is not installed')
class AtomisticMechanics(unittest.TestCase):
    def test_merge_preserves_full_cross_coulomb_energy(self):
        a,b=mm.System(),mm.System()
        for s,q in ((a,1.),(b,-.4)):
            s.addParticle(12.);f=mm.NonbondedForce();f.setNonbondedMethod(mm.NonbondedForce.NoCutoff);f.addParticle(q,.3,0.);s.addForce(f)
        result=t19.merge_systems(a,b)
        integ=mm.VerletIntegrator(.001)
        context=mm.Context(result,integ,mm.Platform.getPlatformByName('Reference'))
        context.setPositions([[0,0,0],[.8,0,0]])
        e=context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        self.assertAlmostEqual(e,138.935456*(-.4)/.8,places=4)
        self.assertEqual(result.getNumParticles(),2)

    def test_real_centroid_force_matches_manual_value_and_gradient(self):
        c=dict(t19.CONFIG);c['metadynamics_steps']=50
        system=mm.System()
        for m in [12.,14.,12.,14.]:system.addParticle(m)
        bias=t19.make_bias_force(c,{'D':[0],'A':[3]},[[0],[1],[2],[3]])
        system.addForce(bias)
        context=mm.Context(system,mm.VerletIntegrator(.001),mm.Platform.getPlatformByName('Reference'))
        xyz=np.array([[1.,0.,0.],[0.,0.,0.],[0.,1.,0.],[.5,1.,1.]])
        d=np.linalg.norm(xyz[0]-xyz[3]);phi=t19.orientation_angle(xyz)
        context.setParameter('h0',.25);context.setParameter('r0',d+.02);context.setParameter('p0',phi+.04)
        context.setPositions(xyz)
        state=context.getState(getEnergy=True,getForces=True)
        observed=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        expected=.25*math.exp(-.5*(.02/.05)**2-.5*(.04/.15)**2)
        self.assertAlmostEqual(observed,expected,places=10)
        force=state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)[0,0]
        step=1e-5
        plus=xyz.copy();plus[0,0]+=step;context.setPositions(plus)
        ep=context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        minus=xyz.copy();minus[0,0]-=step;context.setPositions(minus)
        em=context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        self.assertAlmostEqual(force,-(ep-em)/(2*step),places=6)
        context.setPositions(xyz);context.setParameter('p0',phi+.04+2*math.pi)
        periodic=context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        self.assertAlmostEqual(periodic,observed,places=10)

    def test_bias_configuration_gates(self):
        for change in ({'bias_factor':1.},{'hill_frequency_steps':0},{'metadynamics_steps':999}):
            c=dict(t19.CONFIG);c.update(change)
            with self.assertRaises(ValueError):t19.make_bias_force(c,{'D':[0],'A':[3]},[[0],[1],[2],[3]])

    def test_bias_archive_retains_recorded_nonzero_hills(self):
        c=dict(t19.CONFIG);c['metadynamics_steps']=50
        s=mm.System()
        for _ in range(4):s.addParticle(12.)
        s.addForce(t19.make_bias_force(c,{'D':[0],'A':[3]},[[0],[1],[2],[3]]))
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'biased_system.xml.gz').write_bytes(gzip.compress(mm.XmlSerializer.serialize(s).encode(),mtime=0))
            t19.table(p/'hill_history.csv',[dict(hill=1,height_kJ_mol=.23,center_distance_nm=4.,center_angle_rad=-1.)])
            result=t19.sync_bias_archive(p)
            restored=mm.XmlSerializer.deserialize(gzip.decompress((p/'biased_system.xml.gz').read_bytes()).decode())
            force=restored.getForce(0)
            self.assertEqual([force.getGlobalParameterDefaultValue(i) for i in range(3)],[.23,4.,-1.])
            self.assertEqual(result['default_values_synchronized'],3)
            self.assertEqual(t19.sync_bias_archive(p)['default_values_synchronized'],0)

    def test_checkpoint_continues_positions_and_random_state(self):
        from openmm import app
        def new_simulation():
            top=app.Topology();chain=top.addChain();res=top.addResidue('TEST',chain)
            sys=mm.System()
            for _ in range(2):top.addAtom('C',app.element.carbon,res);sys.addParticle(12.)
            integrator=mm.LangevinMiddleIntegrator(300.,1.,.001);integrator.setRandomNumberSeed(57)
            return app.Simulation(top,sys,integrator,mm.Platform.getPlatformByName('Reference'))
        first=new_simulation();first.context.setPositions([[0,0,0],[.5,0,0]])
        first.context.setVelocitiesToTemperature(300.,61);first.step(4)
        with tempfile.TemporaryDirectory() as temporary:
            source=Path(temporary)/'source';target=Path(temporary)/'resumed';source.mkdir();target.mkdir()
            t19.save_checkpoint(source,first,t19.CONFIG,'test')
            first.step(3)
            a=first.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            second=new_simulation()
            t19.load_checkpoint(target,second,dict(t19.CONFIG,restart_from=str(source)))
            self.assertEqual(second.currentStep,4)
            self.assertTrue((target/'checkpoint_metadata.json').exists())
            second.step(3)
            b=second.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            np.testing.assert_allclose(a,b,rtol=0,atol=1e-12)

    def test_corrupted_or_incompatible_checkpoint_is_rejected(self):
        from openmm import app
        top=app.Topology();chain=top.addChain();res=top.addResidue('TEST',chain);top.addAtom('C',app.element.carbon,res)
        system=mm.System();system.addParticle(12.)
        sim=app.Simulation(top,system,mm.VerletIntegrator(.001),mm.Platform.getPlatformByName('Reference'))
        sim.context.setPositions([[0,0,0]])
        with tempfile.TemporaryDirectory() as temporary:
            source=Path(temporary)/'source';target=Path(temporary)/'resumed';source.mkdir();target.mkdir()
            t19.save_checkpoint(source,sim,t19.CONFIG,'test')
            with self.assertRaises(ValueError):
                t19.load_checkpoint(target,sim,dict(t19.CONFIG,restart_from=str(source),temperature_K=310.))
            with (source/'latest.checkpoint.gz').open('ab') as handle:handle.write(b'corrupted')
            with self.assertRaises(ValueError):t19.load_checkpoint(target,sim,dict(t19.CONFIG,restart_from=str(source)))


if __name__=='__main__':unittest.main()
