# Task 18 — Correct phenotype translation and conditional parent/metabolite PK

This module separates evidence-based genotype labels from hypothetical pharmacokinetic factors. The supplied old CYP2D6 thresholds are corrected, and CYP2C19 is translated independently. It calculates exposure changes and mathematical normalization factors; clinical dosing recommendations remain null. See [sources](sources.md), [summary](outputs/summary.json), [checks](outputs/verification.json) and [figure](outputs/figures/fig18_pgx_cyp2d6_phenotype_kinetics.png).

## Genotype translation

The current CPIC CYP2D6 scheme classifies AS=0 as poor, 0<AS<1.25 as intermediate, 1.25≤AS≤2.25 as normal, and AS>2.25 as ultrarapid. For phased copy-number examples, each allele activity is multiplied by its copy number before summing. The supported allele subset includes *1/*2 at 1, *4/*5 at 0, *10/*41 at 0.25, and *17 at 0.5. Unknown alleles return indeterminate; this is not a comprehensive clinical genotyping system.

CYP2C19 uses allele functions rather than the CYP2D6 activity-score thresholds. In particular *1/*17 is rapid, *17/*17 ultrarapid, and *2/*17 intermediate because an increased-function allele is not assumed to compensate for a no-function allele. Decreased-function combinations return the appropriate “likely” categories. The five common CYP2C19 categories appear in the data even though the overview figure emphasizes four CYP2D6 cohorts.

## PK derivation

The modeled molecule is hypothetical. D is the oral depot, P the parent amount, M a metabolite amount expressed in molar equivalents, and L a cumulative eliminated-equivalent ledger. Absorption ka, bioavailability F, metabolic yield y, parent volume Vp and metabolite volume Vm give:

$$\dot D=-k_aD,$$
$$\dot P=Fk_aD-(CL_{other}+CL_{CYP})P/V_P,$$
$$\dot M=yCL_{CYP}P/V_P-CL_MM/V_M,$$
$$\dot L=(1-F)k_aD+[CL_{other}+(1-y)CL_{CYP}]P/V_P+CL_MM/V_M.$$

Adding all four equations gives zero. Across dose events the total rises by exactly 1000 nmol; molecular masses need not be equal because this is a parent-equivalent molar ledger. Two further states integrate parent and metabolite concentrations in nM·h.

The assumed enzyme pathway is CLCYP=10·fm·g/(1+I/Ki) L/h, while CLother=10(1−fm). Default fm=0.8; g is an explicitly hypothetical phenotype-dependent factor, not activity score itself. The factors PM/IM/NM/RM/UM are 0/0.35/1/1.4/1.8. This distinction prevents turning an ordinal phenotype into a universal empirical clearance law. Constant I/Ki changes modeled enzyme function without changing inherited genotype.

The linear six-state system is propagated by its matrix exponential between seven 24-hour dose events. Dose-boundary rows retain the post-dose depot; the final 168-hour row precedes any new dose. Parent and metabolite AUC comparisons therefore use an identical observation window and input schedule.

## Executed results

Nine diplotypes, each with three inhibitor conditions, produce **27 scenarios and 18,171 trajectory rows**. For uninhibited CYP2D6 PM/IM/NM/UM, parent AUC0–168 is **2700.647/1217.852/594.118/362.775 nM·h**, whereas active-metabolite AUC is **0/663.207/934.816/1032.604 nM·h**. A poor metabolizer can accumulate active parent yet fail to form an active metabolite. Exposure direction depends on which molecular species mediates the response.

An 84-point scan varies enzyme pathway fraction and activity factor. Exposure-normalization factors are AUCreference/AUCscenario under the model's linearity. For a PM with zero bioactivation, no finite input multiplier can match the reference active-metabolite AUC: the corresponding value is null. Even when a finite factor exists, matching one analyte does not establish matched efficacy or safety.

Checks cover CPIC boundary cases, phased duplication, unknown alleles, CYP2C19 rapid/likely categories, mass conservation, nonnegative solutions and timestep-independent matrix propagation. The molar ledger error is below 3e-11 nmol. Deployment as a clinical dosing tool would require an actual drug, verified diplotype/CNV calling, measured genotype-specific PK, active species identification, interacting medications, disease effects and a drug-specific guideline. None is substituted by the synthetic model.
