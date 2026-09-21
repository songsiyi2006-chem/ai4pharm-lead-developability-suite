# Task 18 source and assumption register

Sources checked on 2026-09-21. The supplied prompt's CYP2D6 thresholds are corrected; the same scoring scheme is not applied to CYP2C19.

| Source | Implemented use | Boundary |
|---|---|---|
| [CPIC CYP2D6 / 5-HT3 antagonist guideline, 2026 update, Table 1](https://files.cpicpgx.org/data/guideline/publication/ondansetron/2026/41979467.pdf) | CYP2D6 AS=0 PM; 0<AS<1.25 IM; 1.25≤AS≤2.25 NM; AS>2.25 UM. Selected alleles include *10 and *41 at 0.25. | Allele subset, not a comprehensive clinical genotyping service |
| [CPIC/DPWG consensus translation](https://pmc.ncbi.nlm.nih.gov/articles/PMC6951851/) | Activity-score translation is distinct from substrate-specific pharmacokinetics | No universal clearance multiplier is supplied by AS |
| [CPIC CYP2C19 / clopidogrel guideline, 2022 update, Table 1](https://files.cpicpgx.org/data/guideline/publication/clopidogrel/2022/35034351.pdf) | Independent allele-function categories; *1/*17 rapid, *17/*17 ultrarapid, *2/*17 intermediate; decreased-function “likely” categories | Drug recommendations in the source are not repurposed as generic dosing advice |
| [PharmVar CYP2D6 gene review](https://pmc.ncbi.nlm.nih.gov/articles/PMC6925641/) | Star allele and copy-number interpretation requires structural-variant context | Unknown alleles return indeterminate |

| Computational input | Value | Status |
|---|---|---|
| Example diplotypes | Four CYP2D6 and five CYP2C19 examples | Source-backed category examples; no patient DNA |
| PM/IM/NM/RM/UM CYP clearance factors | 0 / 0.35 / 1 / 1.4 / 1.8 | Assumed for a hypothetical substrate; independent of activity-score arithmetic |
| NM enzyme / other clearance | 8 / 2 L h^-1 | Assumed |
| Parent / metabolite volume | 50 / 30 L | Assumed |
| Metabolite clearance and molar yield | 4 L h^-1, 0.8 | Assumed |
| Absorption, bioavailability | 1 h^-1, 0.85 | Assumed |
| Input schedule | 1000 nmol every 24 h, seven inputs | Scenario definition, not a clinical prescription |
| Competitive inhibitor I/Ki | 0, 1, 100 | Assumed exposure; no actual perpetrator dose |
| NM enzyme contribution sensitivity | 0.3, 0.6, 0.8, 0.95 of 10 L h^-1 | Assumption sensitivity |

Outputs labeled “exposure normalization” are algebraic multipliers within the hypothetical linear model. They are not precision dosing protocols. No therapeutic range, adverse-event risk estimate or patient recommendation is calculated.
