# Task 15 source and assumption register

Sources checked on 2026-09-21. Literature establishes mechanisms; no patient dataset or fitted cytokine parameter was supplied.

| Source | Supported use | Not inferred |
|---|---|---|
| [McKeithan 1995, original kinetic proofreading model](https://pmc.ncbi.nlm.nih.gov/articles/PMC41844/) | Sequential receptor modification with ligand dissociation resetting the cascade | The five stages or numerical rates of a particular engineered receptor construct |
| [Britain et al. 2022, optogenetic proofreading experiment](https://pmc.ncbi.nlm.nih.gov/articles/PMC9536835/) | Experimental discrimination depends on ligand dwell time and signaling level | Direct parameter transfer to patient-specific signaling |
| [Norelli et al. 2018, primary mouse-model study](https://pubmed.ncbi.nlm.nih.gov/29808007/) | Monocyte-derived cytokines contribute to the modeled inflammatory feedback mechanisms | Human cytokine homeostasis or a patient treatment dose |

| Model input | Value / unit | Evidence |
|---|---|---|
| Receptor total, cascade stages | 1 normalized receptor pool; 5 stages | Assumed |
| Ligand pseudo-first-order association | 90, 90 h^-1 | Assumed reservoir concentrations times association constants |
| Phosphorylation; target/self dissociation | 720; 36/360 h^-1 | Assumed |
| Target antigen decline | Starts 24 h; rate 0.08 h^-1 | Assumed external forcing; self antigen fixed |
| T-cell activation, decay | 0.22, 0.06 h^-1; trigger half-response 0.15 receptor fraction | Assumed |
| Macrophage activation through IFN / IL6R; decay | 0.35 / 0.20; 0.10 h^-1 | Assumed |
| IL6 production from T / M | 0.18 / 1.10 normalized concentration h^-1 | Assumed |
| IL6 independent / receptor-linked clearance | 0.12 / 0.08 h^-1 | Assumed |
| TNF production from T / M; clearance | 0.55 / 0.60; 0.30 h^-1 | Assumed normalized units |
| IFN production from T / M; clearance | 0.90 / 0.10; 0.25 h^-1 | Assumed normalized units |
| IL1 production from T / M; clearance | 0.08 / 0.80; 0.20 h^-1 | Assumed normalized units |
| Cytokine activation half-response constants | 1 normalized concentration | Assumed |
| Antagonist concentration / KD | 0, 1, 10, 100; begins 12 h; half-life 48 h | Assumed, not tocilizumab dosage |

No measured numerical biological parameter is claimed. Blood pressure, organ failure and clinical intervention recommendations are not model outputs.
