# Task 19 — Real ternary-complex atomistics and a bounded WTMetaD pilot

This project uses experimental **5T35 BRD4–MZ1–VHL/Elongin B/Elongin C** coordinates, complete classical atomistic parameters, explicit TIP3P solvent and a two-variable well-tempered metadynamics engine. The deliverable is a local CPU **engine pilot**, with a strict gate against calling its short bias history a converged PMF, binding free energy, cooperativity factor or productive-degradation basin.

The archived calculation completed: **150,838 atoms, 0.2 ps initial thermalization + 1.0 ps biased dynamics, 20 real hills, one CPU thread, 16.90 min wall time**. Positions and recorded energies stayed finite. The PMF is unconverged and cooperativity remains `null`.

| Locate | Contents |
| --- | --- |
| [driver.py](driver.py) | Structure audit, preparation, system merge, actual CPU dynamics, analytic WTMetaD hills, figure and manifest |
| [SOURCES.md](SOURCES.md) | Primary literature, force-field and source-file provenance |
| [REPORT_ZH.md](REPORT_ZH.md) / [REPORT_EN.md](REPORT_EN.md) | Methods, equations, evidence boundaries and production requirements |
| [inputs](inputs) | 5T35/CCD originals and compressed, hashed prepared all-atom system |
| [outputs/summary.json](outputs/summary.json) | Actual execution status, duration, particle count and deliberately null free-energy endpoints |
| [outputs/hill_history.csv](outputs/hill_history.csv) | Time, center, deposited height and prior bias for every real hill |
| [outputs/bias_archive_audit.json](outputs/bias_archive_audit.json) | Audit/synchronization of serialized force defaults against actual recorded hills |
| [outputs/executed_driver.py](outputs/executed_driver.py) | Source snapshot of the actual numerical execution |
| [outputs/atomistic_cv_trajectory.csv](outputs/atomistic_cv_trajectory.csv) | Real MD distance/orientation and energy samples |
| [outputs/figures/fig19_ternary_metadynamics_pmf_landscape.png](outputs/figures/fig19_ternary_metadynamics_pmf_landscape.png) | 300 DPI structure/trajectory/bias/evidence figure; filename follows the request, panel explicitly states that the bias is not a converged PMF |

From the repository root, use an interpreter with OpenMM installed:

```powershell
$env:OPENMM_CPU_THREADS='1'
$env:OMP_NUM_THREADS='1'
python projects/task19_metadynamics/driver.py --out work/task19_reproduce
```

The optional pip dependencies for this module are declared in [requirements.txt](requirements.txt). The archived run uses OpenMM 8.6.0. Fresh-system reruns use the OpenMM 8.6 series; binary checkpoint continuation requires the matching version and platform. Preparation from the original crystal additionally uses RDKit, OpenFF Toolkit/Interchange with Sage 2.3 and Ash 1.0, and PDBFixer; those are not required to integrate the published parameterized system.

The output directory must be new or empty. Running this local pilot requires no HPC. It does not train a model or approximate proteins as beads. The root advanced-suite entry point delegates to this module. A single CPU thread and a bounded number of steps control resource use; stochastic/CPU/library differences mean MD trajectories are not guaranteed byte-identical across machines.

The verified parameter cache saves repeating preparation. It includes the full physical system, including solvent and force terms; it is not a cached synthetic PMF. A full rebuild can be executed with `prepare_system(Path(...), CONFIG, ignore_cache=True)` after creating an empty destination and loading this module. The installed preparation environment is recorded in the reports. No dependency is installed by the driver.

The driver writes numerical results and a verified binary checkpoint after every hill, and checkpoints after minimization and initial thermalization. An interrupted compatible local run can be continued into a **new** output directory without repeating completed dynamics:

```powershell
python projects/task19_metadynamics/driver.py --restart-from work/task19_interrupted --out work/task19_resumed
```

Configuration, prepared-system identity, step counter and checkpoint-file hashes must agree. Binary checkpoints require the same OpenMM version and compatible CPU platform; they are not a cross-platform trajectory format. The latest checkpoint metadata is published after its numerical files. A partial checkpoint transaction fails verification instead of silently continuing from inconsistent data.
