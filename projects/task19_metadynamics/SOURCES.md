# Task 19 sources and preparation provenance

All sources below were inspected on 2026-09-21. The calculation uses archived structure files; it does not fetch a replacement on rerun.

| Source | Role and boundary |
| --- | --- |
| [RCSB 5T35](https://www.rcsb.org/structure/5T35) and [PDB archive coordinates](https://files.rcsb.org/download/5T35.pdb) | Experimental X-ray structure at 2.70 angstrom. Author chains A/B/C/D are BRD4 BD2, Elongin B, Elongin C and VHL. The calculation uses the first assembly and MZ1, CCD 759, author chain D. These are crystal coordinates, not a molecular-dynamics trajectory. |
| [Gadd et al., 2017](https://doi.org/10.1038/nchembio.2329) | Original ternary-complex structural, mutagenesis and binding experiments. The existence of published cooperativity measurements does not make the present simulation a calculation of those measurements. |
| [CCD 759](https://www.rcsb.org/ligand/759), [CCD dictionary](https://files.rcsb.org/ligands/download/759.cif), [ideal SDF](https://files.rcsb.org/ligands/download/759_ideal.sdf) | Bond orders and stereochemical identity of MZ1. Crystal and CCD canonical isomeric SMILES must match before parameterization. |
| [OpenFF force-field release registry](https://github.com/openforcefield/openff-forcefields) | The installed `openff-2.3.0.offxml` is Sage 2.3 and uses NAGL charges. The standard file itself declares the Ash 1.0 model and hash. No substitute coarse-grained ligand potential is used. |
| [Ash / AM1-BCC GNN 1.0](https://docs.openforcefield.org/projects/nagl-models/en/latest/models/openff-gnn-am1bcc-1.0.0/) | Learned partial-charge model. Charges resemble the AM1-BCC ELF10 target; they are not a fresh semiempirical or ab-initio charge calculation for MZ1. |
| [OpenMM simulation guide](https://docs.openmm.org/latest/userguide/application/02_running_sims.html) | Amber protein force fields, explicit water, PME, constraints, Langevin dynamics and custom force support. |
| [Barducci, Bussi and Parrinello, 2008](https://doi.org/10.1103/PhysRevLett.100.020603) | Well-tempered metadynamics and its asymptotic bias/free-energy relation. The relation does not guarantee convergence in a picosecond pilot. |

The source-file SHA256 values are in `outputs/structure_audit.json`. The full prepared system and coordinates are compressed in `inputs/prepared/`; that directory's manifest verifies their exact bytes. Preparation is executable through `prepare_system(out, CONFIG, ignore_cache=True)` and requires the already installed OpenFF, RDKit, PDBFixer and OpenMM environments. Ordinary reruns reuse these audited atomistic parameters and need OpenMM, NumPy, SciPy and Matplotlib.

The original PDB SEQRES records are retained during preparation. Unresolved terminal tails are omitted, and the internal Elongin C gap is built by PDBFixer. Its modeled geometry is not experimental data. All hydrogenation uses the stated pH 7 default variants; no protein pKa prediction or protonation ensemble was performed.

The docking/structural contacts identify geometric candidates for future mutation work only. No mutation is reported as thermodynamic resistance without separate mutation free-energy calculations and experiments.
