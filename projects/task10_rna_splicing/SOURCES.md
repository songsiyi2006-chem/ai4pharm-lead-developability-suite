# Sources, parameters, and provenance

Sources checked on 2026-09-20. Public structure files are redistributed as PDB archive data with attribution; see [wwPDB data access](https://www.wwpdb.org/about/usage). Article text and published figures have not been copied. All four figures are computed here.

| Source | What it supports | What it does not supply here |
|---|---|---|
| [Campagne et al., 2019, DOI 10.1038/s41589-019-0384-5](https://doi.org/10.1038/s41589-019-0384-5) | SMN-C5, bulged adenine, RNA duplex and splicing correction mechanism | The four hypothetical scaffold free-energy components |
| [RCSB 6HMI](https://www.rcsb.org/structure/6HMI) / [raw PDB](https://files.rcsb.org/download/6HMI.pdb) | Apo RNA coordinates, 20 submitted NMR conformers, 22 residues | Boltzmann populations or dynamic time ordering |
| [RCSB 6HMO](https://www.rcsb.org/structure/6HMO) / [raw PDB](https://files.rcsb.org/download/6HMO.pdb) | SMN-C5-bound RNA coordinates, 20 submitted conformers | Clinical risdiplam exposure or affinity of our model compounds |
| [ViennaRNA partition-function documentation](https://viennarna.readthedocs.io/en/latest/partfunc/global.html) and [Python API](https://viennarna.readthedocs.io/en/latest/api_python.html) | Genuine partition and BPP functions | Tertiary flipping barriers or spliceosome activity |
| [McCaskill, 1990, DOI 10.1002/bip.360290621](https://doi.org/10.1002/bip.360290621) | Equilibrium secondary-structure partition-function formalism | A measurement on this engineered sequence |
| [White et al., 2024, DOI 10.1038/s41467-024-53124-5](https://www.nature.com/articles/s41467-024-53124-5) | Reconstituted U1 snRNP, U1-C and branaplam-dependent recognition, industry collaboration | Universal kinetic order for all splice modulators |
| [Specificity, synergy, and mechanisms of splice-modifying drugs, 2024](https://www.nature.com/articles/s41467-024-46090-5) | Sequence/context-specific drug responses and mechanistic distinctions | Clinical safety of a single generic off-target pool |

## Archived structures

| File | SHA256 |
|---|---|
| `inputs/6HMI.pdb` | `36941bba2a15d0bd4aa75d32bb1417d193e2a940bb9220e307edf4f4100b2f7c` |
| `inputs/6HMO.pdb` | `9ef5a62593c3356424f386c066b4332fe7e745f83407281be49f740bf7eb17aa` |

Chain A sequence: `AUACΨΨACCUG`; chain B: `GGAGUAAGUCU`. Author residue numbering is A1–A11 and B12–B22; bulged adenine is B14 (third nucleotide of chain B). The 3D site definitions use author residue IDs, including A3/B20, A4/B19, A5/B18, A6/B17, A7/B16, A8/B15, A9/B13. The secondary-structure model concatenates A (Ψ→U), UUCG, and B into one 26-nt chain. It is an engineered intramolecular control, not an experimentally validated SMN2 hairpin.

## Default assumptions and units

All non-coordinate scenario parameters are in [example_config.json](inputs/example_config.json). RNA folding uses ViennaRNA 2.7.2 default Turner parameter set, `dangles=2`, explicitly set temperature and monovalent salt 0.15 M. The NMR acquisition conditions are not replaced by this calculation salt.

| Parameter | Default | Status |
|---|---:|---|
| Temperature | 37 C | Calculation condition |
| U1 total; target RNA; off-target RNA | 100; 2; 30 nM | Hypothetical finite-pool concentrations |
| U1 KD target/off-target | 2 / 4 uM | Assumed |
| Opening penalties, four scaffolds | 5.5 / 4.0 / 5.5 / 6.5 kcal/mol | Distinct ligand-competent substates, not extracted from BPP |
| Stack + H-bond + electrostatic + desolvation + entropy | Per-ligand fields in config | Assumed decomposition; no docking/FEP/DFT energies |
| WT splice stabilization | -3.5 / -2 / 0 / -3 kcal/mol | Assumed, shared reciprocal coupling factor |
| Mutant scenario | +1.5 kcal opening; intrinsic KD ×10; U1 KD ×3; coupling no stronger than -1 kcal/mol | Independent counterfactual assumptions, not computed from the A18C sequence mutation |
| Off-target scenario | Opening penalty reduced by 1.5 (floor 0); ligand-specific KD multiplier/coupling | Assumed aggregate RNA pool |
| Four scaffold names | Planar, polycationic, inactive, flexible | Qualitative models, not four synthesized/assayed molecules |

The dielectric is fixed at 78.5, each phosphate P is a -1 point charge, and screening uses a monovalent Debye–Hückel expression. Waters, counterion condensation, Mg2+, atomic partial charges and dielectric boundaries are omitted. Values exceeding thermal voltage are outside a quantitatively reliable linear-response regime. These are signed electrostatic descriptors, not binding energies or Poisson–Boltzmann solutions.
