# Task 20: Three-dimensional pocket-conditioned generation

The local CPU run used the public unliganded SARS-CoV-2 Mpro structure 6YB7. Fragment graphs were assembled, ETKDG/MMFF conformers were generated, 40 random rigid poses were scored and translations were locally refined. NSGA-II evolved 100 candidates for 20 generations using a Pareto balance of pocket-geometry contact proxy, a CNS-MPO proxy, QED and RDKit SA; an equal unique-graph budget was evaluated in random order as a control.

## Executed results

The evolutionary run evaluated 1,034 unique graphs; the random control expanded the cache to 1,597 records. The best geometry score changed from 19.21995 to 20.66983, while the random-control best was 20.65233. The small difference does not establish an evolutionary advantage. All 100 final candidates passed the model's SA<3.5 proxy threshold. AutoDock Vina independently docked 12 candidates with a rigid receptor. Two independent X77 redocking controls in 6W63 gave unaligned heavy-atom RMSDs of 1.242 Å and 1.055 Å.

## Equations and evidence boundary

The contact term uses a Gaussian function of atomic van der Waals gaps and an overlap penalty normalized by ligand heavy atoms; the geometry score is dimensionless and ranking-only. The CNS-MPO proxy uses structural heuristic pKa values 0/5/9 and has no experimental or trained calibration. SA is the Ertl fragment-complexity proxy, not a synthesis-route review. Vina scores were not converted to Kd or nM.

6YB7 is biologically dimeric, while generation used a rigid chain A. Crystal waters, receptor flexibility, microstate ensembles, long-timescale sampling and experimental binding measurements are absent. Novelty, synthetic feasibility, cellular efficacy and clinical utility remain unknown. Acceptance requires multi-conformer/multi-microstate docking, replicated seeds, synthesis and binding/function assays.

Inputs, sources and SHA256 records are in [README](README.md), [sources.md](sources.md) and [outputs/summary.json](outputs/summary.json). Figure: [fig20_denovo_pareto_lead_optimization.png](outputs/figures/fig20_denovo_pareto_lead_optimization.png).
