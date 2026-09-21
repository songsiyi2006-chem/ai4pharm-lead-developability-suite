# AI4Pharm advanced modules: delivery, biological context and molecular design

**Tasks 11–20 · Integrated technical report · 21 September 2026**  
[中文版](AI4PHARM_ADVANCED_TREATISE_ZH.md) · [Tasks 1–10](AI4PHARM_DECADE_TREATISE_EN.md) · [Repository navigation](README.md)

This report integrates the equations, executed local calculations, parameter provenance and acceptance criteria of eight delivered advanced modules; Tasks 11 and 17 are skipped under the revised scope. The evidence is heterogeneous: public structural coordinates and source-backed pharmacogenetic categories coexist with hypothetical kinetics, finite-pool simulations, a short atomistic pilot and uncalibrated molecular-design scores. Integration makes these differences inspectable; it does not turn them into a clinically calibrated model.

The governing research requirements are reproducible inputs, correct units, explicit conservation, source attribution, numerical checks and conclusions supported by the observed outputs. A requested optimum, mechanism or free-energy result is retained as an unmet objective when the calculation does not establish it. No supercomputer is assumed. Parallel agents can develop and independently review bounded modules while CPU jobs retain explicit resource limits.

## Navigation

[11 LNP](#task11) · [12 Chromatin degradation](#task12) · [13 BBB–TfR](#task13) · [14 BiTE](#task14) · [15 Receptor signaling](#task15) · [16 Transporters](#task16) · [17 Pulmonary delivery](#task17) · [18 Pharmacogenetics](#task18) · [19 Atomistic sampling](#task19) · [20 Molecular generation](#task20) · [Integrated interpretation](#integration)

<a id="task11"></a>
## 11. LNP module — skipped in this release

Task 11 is excluded from the current delivery under the user's revised scope. This report contains no Task 11 model, calculations or figure. The task number is retained for navigation. [Status](projects/task11_lnp/README.md)
<a id="task12"></a>
## 12. Chromatin-context competition and target degradation

### Binding thermodynamics and the kinetic network

Unliganded free-context and chromatin-context targets are $F,C$, with degrader $P$ and ligase $E$. Binary states are $PF,PC,EP$; ternary states are $EFP,ECP$. Binding-only equilibrium satisfies

$$PF=pf/K_T,\quad PC=pc/K_T,\quad EP=ep/K_E,$$
$$EFP=\alpha_Fepf/(K_EK_T),\quad ECP=\alpha_Cepc/(K_EK_T),\quad\alpha_C=\alpha_Fe^{-\Delta G_{occ}/RT}.$$

Finite totals include every species containing the relevant component. Eliminating $f=F_{tot}/[1+p/K_T+\alpha_Fep/(K_EK_T)]$, and analogously $c$, reduces equilibrium to bounded roots for free $p,e$. The nominal $K_E=50$ nM, $K_T=100$ nM, $\alpha_F=10$ and assumed occlusion penalty 2 kcal/mol give $\alpha_C=0.38968$.

The degrading system is integrated kinetically rather than assumed to remain at equilibrium. For each reversible reaction $r$, $\dot x=\sum_r\nu_r(v_{r,+}-v_{r,-})$. Binary off-rates equal $k_{on}K_D$; the appropriate ternary off-rate is divided by $\alpha_F$ or $\alpha_C$. Both paths, $EP+F\rightleftharpoons EFP$ and $PF+E\rightleftharpoons EFP$, therefore have the same equilibrium product. The same construction applies to chromatin-context species.

Only unliganded targets exchange contexts, $J_{chrom}=k_{cap}F-k_{rel}C$. All target-containing species turn over at $k_0$, releasing surviving ligase and degrader; synthesis enters $F$ at $s=k_0T_0$. Ternary-specific degradation adds $k_{cat}(EFP+ECP)$. Summation cancels binding and context exchange:

$$\dot T=s-k_0T-k_{cat}(EFP+ECP),\qquad T+D-S=T_0.$$

$D,S$ are cumulative degraded and synthesized target. Degrader and ligase totals remain finite and conserved. With no drug, $T=T_0$ and

$$C(t)=C_\infty+[C(0)-C_\infty]e^{-\lambda t},\quad\lambda=k_{cap}+k_{rel}+k_0,\quad C_\infty=k_{cap}T_0/\lambda.$$

This independent baseline explains why a changed chromatin fraction need not imply drug-induced degradation.

### Executed result and development interface

At 100 nM degrader and 48 h, total target is **51.2557 nM**, comprising **22.7889 nM free-context** and **28.4668 nM chromatin-context** pools. Each pool includes its bound species. Untreated target redistributes from 30:70 to **37.7776:62.2224 nM** without losing total mass. There are **152 dynamic scenarios and 152 binding equilibria**, spanning 38 doses and four penalties; the largest grid balance error is approximately **1.29×10⁻⁹ nM**.

Numerical biological parameters are hypothetical. There is no nucleosome coordinate calculation, measured steric penalty, finite chromatin-site occupancy or explicit multimeric complex disassembly. Experimental selective BET degradation motivates the question but does not provide these parameter values. [Zengerle et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC4548256/) Industrially useful next measurements are matched soluble/chromatin fractionation, target synthesis/turnover, ligase abundance, washout and accessibility perturbations. The energy penalty is a testable input, not a computed nucleosome uncoupling free energy.

[Parameters](projects/task12_chromatin/outputs/config.json) · [Results](projects/task12_chromatin/outputs/summary.json) · [Sources](projects/task12_chromatin/sources.md) · [Full report](projects/task12_chromatin/REPORT_EN.md)

![Task 12: finite-ligase binding and context-dependent target depletion](figures_advanced/fig12_epigenetic_chromatin_depletion.png)

<a id="task13"></a>
## 13. TfR transport across the BBB: affinity, release and assay definition

### Finite antibody and receptor amounts

Amounts are nmol and volumes L, so amount/volume is nM. States are plasma antibody $A_p$, unoccupied surface receptor $R_s$, surface complex $C_s$, endosomal receptor $R_e$, endosomal complex $C_e$, free endosomal antibody $A_e$ and free brain antibody $A_b$. Net association fluxes are

$$J_s=k_{on}[(A_p/V_p)R_s-K_{D,s}C_s],\qquad J_e=k_{on}[(A_e/V_e)R_e-K_{D,e}(t)C_e].$$

The common pH closure $pH_e=5.8+1.6e^{-t/(1h)}$ modifies $K_{D,e}=K_{D,s}10^{q(7.4-pH_e)}$, with assumed slope $q=0.7$. This is a specified acid-sensitive binding hypothesis, not measured histidine-mediated release.

Define complex/receptor internalization $I_C=k_iC_s$, $I_R=k_{ir}R_s$, receptor recycling $Q=k_rR_e$, complex/free-cargo degradation $D_C=k_{dc}C_e$, $D_F=k_{df}A_e$, brain/plasma release $T_b=k_bA_e$, $T_p=k_pA_e$, and clearances $L_p=k_{cl,p}A_p$, $L_b=k_{cl,b}A_b$:

$$\dot A_p=-J_s+T_p-L_p,\quad\dot R_s=-J_s-I_R+Q,\quad\dot C_s=J_s-I_C,$$
$$\dot R_e=I_R-Q-J_e+D_C,\quad\dot C_e=I_C+J_e-D_C,$$
$$\dot A_e=-J_e-D_F-T_b-T_p,\quad\dot A_b=T_b-L_b.$$

The corresponding cumulative losses receive $L_p$, $D_C+D_F$ and $L_b$. Their sum with antibody inventories equals initial antibody, and $R_s+C_s+R_e+C_e=R_{tot}$. Degrading cargo returns receptor to the endosomal pool. Thus receptor downregulation is intentionally absent; it must be added before interpreting systems where trafficking also destroys TfR.

### Exposure and executed outcome

The observable is explicitly finite-time and extravascular:

$$K_{p,0-72h}=\frac{\int_0^{72h}A_b/V_b\,dt}{\int_0^{72h}A_p/V_p\,dt}.$$

Nominal ratio is **0.00227242**, brain AUC **6.69954 nM·h**. Because there is no brain target-binding pool, this is not a general $K_{p,uu}$ prediction for target-binding antibodies. Adding assumed vascular contamination $C_pV_{vascular}/V_b$ increases the apparent ratio to **0.0522724**, exactly **0.05** higher. This illustrates a measurement-definition error, not improved parenchymal delivery.

The **328 ODE scenarios** are one-factor slices through affinity, pH sensitivity, dose and receptor amount, not their full Cartesian product. The nominal affinity optimum is **13.3352 nM**, outside the requested 50–500 nM range. Zero receptor or complex internalization gives zero delivery; the independent no-receptor solution is $A_p=A_p(0)e^{-0.03t}$.

All numerical trafficking inputs are assumptions. Reduced-affinity BBB uptake has experimental precedent in a specific mouse system, but no universal optimum follows. [Yu et al.](https://pubmed.ncbi.nlm.nih.gov/21613623/) [TfR trafficking study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3920563/) Translation requires species-specific receptor abundance, endogenous transferrin competition, valency/avidity, receptor turnover, capillary-depleted exposure and target engagement. The model's main value is exposing release/retention and measurement trade-offs.

[Parameters](projects/task13_bbb_tfr/outputs/config.json) · [Results](projects/task13_bbb_tfr/outputs/summary.json) · [Sources](projects/task13_bbb_tfr/sources.md) · [Full report](projects/task13_bbb_tfr/REPORT_EN.md)

![Task 13: affinity-dependent trafficking and extravascular exposure](figures_advanced/fig13_bbb_tfr_transcytosis_profile.png)

<a id="task14"></a>
## 14. BiTE crosslinking and target-cell loss

### Equilibrium stoichiometry and cellular units

Free CD3, tumor antigen and BiTE concentrations are $e,t,b$ in nM. Binary states are $EB=eb/K_E$, $TB=tb/K_T$, and molecular trimer $X=\alpha etb/(K_EK_T)$. The finite balances are

$$E_0=e+EB+X,\qquad T_0=t+TB+X,\qquad B_0=b+EB+TB+X.$$

Nested bounded roots enforce these balances without equating free and total BiTE. Exchanging receptor labels preserves the trimer; removing one receptor recovers a binary quadratic limit. These are independent numerical checks of the topology.

An assumed contact volume $v=10^{-12}$ L per initial target links accessible receptor copies to concentration. If $n=N/N_0$, effector/initial-target ratio is $r$ and copy counts are $n_{CD3},n_{TAA}$,

$$E_0=\frac{r n_{CD3}}{N_Av}10^9,\qquad T_0(n)=\frac{n n_{TAA}}{N_Av}10^9\quad[\mathrm{nM}],$$
$$S(n)=X(n)10^{-9}N_Av/n.$$

$S$ counts molecular trimers per surviving target, not cell–cell immunological synapses. An actual synapse population would need spatial contacts, interacting cell counts and assembly kinetics. With declared $k_{max}=0.12$ h⁻¹, $K_S=500$ complexes/cell and $h=2$,

$$\dot n=-k_{max}\frac{S(n)^h}{K_S^h+S(n)^h}n,\qquad\dot d=-\dot n,\qquad n+d=1.$$

Effector count and externally maintained total BiTE are fixed. Lost target receptors release ligand into that maintained pool; antibody PK and irreversible drug loss are outside this model.

### Executed outcome and development interface

At CD3 $K_D=100$ nM and TAA $K_D=1$ nM, the sampled maximum is **26.1016 nM total BiTE**, with **572.752 trimers/target**. At 100 µM, trimer abundance falls to **0.144815%** of peak. Maintained 10 nM exposure leaves **3.82331%** of initial targets at 48 h. The run contains **981 binding equilibria and 270 lysis ODEs**. Zero drug, effector or crosslink cooperativity gives zero lysis.

The high-dose hook is a consequence of this reversible finite-receptor model. It is not a universal clinical downturn above 1 µM. All affinities, copy numbers, contact volume and lysis parameters are assumptions; the model does not include cytokine toxicity, serial killing kinetics, exhaustion, antigen shedding or cell motility. Primary mechanistic models support asking these questions without supplying the present rates. [Betts et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6531394/) [Cell-level synapse population dynamics](https://elifesciences.org/articles/83659) The next industrial interface is matched receptor quantification, binding kinetics, cell-pair imaging and dose/time killing assays.

[Parameters](projects/task14_bite/outputs/config.json) · [Results](projects/task14_bite/outputs/summary.json) · [Sources](projects/task14_bite/sources.md) · [Full report](projects/task14_bite/REPORT_EN.md)

![Task 14: molecular crosslinking hook and conditional target-cell loss](figures_advanced/fig14_bite_synapse_crosslinking_curve.png)

<a id="task15"></a>
## 15. Finite-receptor proofreading and cytokine signaling

### Dwell-time discrimination with a conserved receptor pool

For antigen class $a$, a reservoir creates pseudo-first-order binding $k_aR$. Five phosphorylation transitions proceed at $k_p$; all bound states dissociate at $k_{off,a}$ and return receptor to $R$:

$$\dot C_{a,0}=k_aR-(k_p+k_{off,a})C_{a,0},$$
$$\dot C_{a,j}=k_pC_{a,j-1}-(k_p+k_{off,a})C_{a,j},\quad j=1,\ldots,4,$$
$$\dot C_{a,5}=k_pC_{a,4}-k_{off,a}C_{a,5},\qquad\dot R=\sum_a k_{off,a}\sum_{j=0}^5C_{a,j}-R\sum_a k_a.$$

Summation gives $R+\sum_{a,j}C_{a,j}=1$. Ligand is an external reservoir, with target antigen declining after 24 h. At steady state, $R=[1+\sum_a k_a/k_{off,a}]^{-1}$; the total bound pool for antigen $a$ is $k_aR/k_{off,a}$. Each of five successful modification steps has probability $k_p/(k_p+k_{off,a})$, so the terminal fraction among bound receptors is its fifth power. This analytic limit tests the implemented kinetic proofreading topology. [McKeithan](https://pmc.ncbi.nlm.nih.gov/articles/PMC41844/)

Let $T,M$ be active T-cell/macrophage fractions and $I_6,N,F,I_1$ normalized cytokines. Set $q=(C_{target,5}+C_{self,5})/(0.15+C_{target,5}+C_{self,5})$, antagonist occupancy $B=X/(1+X)$ and IL6R signal $s=I_6(1-B)/(1+I_6)$. The implemented network, with time in hours, is

$$\dot T=0.22q(1-T)-0.06T,\quad\dot M=[0.35F/(1+F)+0.20s](1-M)-0.10M,$$
$$\dot I_6=0.18T+1.10M-[0.12+0.08(1-B)]I_6,$$
$$\dot N=0.55T+0.60M-0.30N,\quad\dot F=0.90T+0.10M-0.25F,\quad\dot I_1=0.08T+0.80M-0.20I_1.$$

The active-fraction derivatives point inward at zero and one. All coefficients and concentration units are assumptions. Antagonist forcing is zero until 12 h and then $X_0e^{-\ln2(t-12)/48}$; the solver segments that intervention boundary.

### Results and meaningful endpoints

The assumed target/self terminal-state ratio is **59.499**, not measured construct specificity. Four 96 h scenarios contain **3,844 rows**, with a **65-point** self-antigen off-rate scan. For $X_0=0,1,10,100$, IL6R-signal AUC is **75.201, 51.158, 17.106, 6.781 h**, whereas IL6 peaks are **5.108, 5.944, 7.186, 7.587 normalized units**. Blocking receptor-linked removal allows lower signaling and higher cytokine concentration simultaneously; it does not abolish IL6 production.

There is no blood pressure, endothelial leak, shock, organ injury or clinical toxicity probability in this model. `homeostasis_restoration_established` remains null. Dimensionless $C/K_D$ is not a tocilizumab dose. Literature on cytokines in specific experimental systems supports mechanistic distinctions, not human rescue claims. [Norelli et al.](https://pubmed.ncbi.nlm.nih.gov/29808007/) Validation needs construct-specific binding/phosphorylation, antigen density, cytokine units, macrophage assays and an independently checked antagonist exposure model. This module provides an assay-design and signaling-hypothesis interface.

[Parameters](projects/task15_receptor_signaling/outputs/config.json) · [Results](projects/task15_receptor_signaling/outputs/summary.json) · [Sources](projects/task15_receptor_signaling/sources.md) · [Full report](projects/task15_receptor_signaling/REPORT_EN.md)

![Task 15: receptor proofreading, cytokine concentration and receptor signal](figures_advanced/fig15_receptor_proofreading_cytokine_balance.png)

<a id="task16"></a>
## 16. Transporter-limited cellular access inside the existing PBPK circulation

### Actual reuse and avoidance of double counting

This module imports Task 4's `PBPK.derivative`, retaining its series lung/systemic circulation and glomerular filtration. It first reverses the old hepatic metabolic flux in both the liver and loss derivatives. The new intracellular pathway therefore replaces that mechanism instead of adding a second parallel elimination sink. Seven original amount states gain finite hepatocyte, renal-cell and enterocyte states; inherited partition factors remain approximations.

For unbound donor concentration $C_{u,ext}$, hepatic carrier uptake plus passive exchange is

$$J_H=\sum_{j\in\{1B1,1B3\}}\frac{V_{max,j}C_{u,ext}}{K_{m,j}(1+I/K_i)+C_{u,ext}}+PS_H(C_{u,ext}-C_{u,H}).$$

$C_{u,H}=f_{u,H}A_H/V_H$; saturable terms have mass/time units and passive $PS$ has volume/time units. Competitive inhibition increases apparent $K_m$ while leaving limiting $V_{max}$ unchanged. Passive exchange can reverse sign, so it is not a one-way sink. The hepatocyte equation is

$$\dot A_H=J_H-CL_{met}C_{u,H}-\frac{V_{max,bile}C_{u,H}}{K_{m,bile}+C_{u,H}}.$$

Uptake is subtracted from the donor liver space. Metabolism and bile go to distinct ledgers. For renal cells, $\dot A_R=J_R-J_{urine}$ uses separate renal uptake/passive terms and saturable P-gp/BCRP efflux. Existing filtration remains exactly once; renal uptake is not assigned to hepatic OATP.

The old gut-depot→liver shortcut is removed. With luminal amount $A_g$, enterocyte amount $A_E$, luminal-return efflux $J_{eff}$, basolateral rate $k_b$ and transit rate $k_f$,

$$\dot A_E=k_aA_g-J_{eff}-k_bA_E,\qquad\dot A_g=-k_aA_g+J_{eff}-k_fA_g.$$

Basolateral export enters liver; transit enters fecal loss. Every exchange appears with opposite signs in donor and receiver. Summing gives $\sum A_i+L_{met}+L_{urine}+L_{feces}=Dose$. The AUC state has concentration×time units and is excluded from this mass sum. The exact implementation is [TransportPBPK.derivative](projects/task16_transporters/driver.py).

### Executed exposure and limitations

Eight detailed IV/oral scenarios compare $I/K_i=0,1,10,100$ at 100 mg; **24** additional dose/inhibition combinations inspect saturation. IV AUC0–96 values are **7.266, 10.250, 16.399, 18.687 mg·h/L**, corresponding to AUCR **1, 1.411, 2.257, 2.572**. Oral AUCR is **1, 1.495, 2.533, 2.919**. At IV $I/K_i=100$, **55.72 mg remains at 96 h**. These are truncated AUC ratios, not complete exposure ratios.

All transporter kinetics and added cell volumes are assumptions; a constant imposed inhibitor ratio has no perpetrator PK. Adding effective cellular pools does not establish anatomically validated volume decomposition. Numerical checks address stoichiometry, empty-donor limits, competition and refinement. Industrial use requires measured transporter $V_{max}/K_m$, intracellular binding, abundance/scaling and external PK data. Transporter-mediated interactions belong in development assessment, but the guideline does not validate this scenario. [FDA/ICH M12](https://www.fda.gov/regulatory-information/search-fda-guidance-documents/m12-drug-interaction-studies)

[Parameters](projects/task16_transporters/outputs/config.json) · [Results](projects/task16_transporters/outputs/summary.json) · [Sources](projects/task16_transporters/sources.md) · [Full report](projects/task16_transporters/REPORT_EN.md)

![Task 16: saturable uptake, efflux and finite-time PK](figures_advanced/fig16_transporter_oatp_pgp_kinetics.png)

<a id="task17"></a>
## 17. Pulmonary module — skipped in this release

Task 17 is excluded from the current delivery under the user's revised scope. This report contains no Task 17 model, calculations or figure. The task number is retained for navigation. [Status](projects/task17_pulmonary/README.md)
<a id="task18"></a>
## 18. Pharmacogenetic translation and parent–metabolite exposure

### Evidence-based categories are not universal clearance laws

The source-backed CYP2D6 activity-score translation is PM at AS=0, IM at $0<AS<1.25$, NM at $1.25\le AS\le2.25$, and UM at AS>2.25. Phased allele copy numbers multiply allele activity before summation; the supported subset assigns *10 and *41 a value of 0.25. Unknown alleles return indeterminate. This is a limited translator, not a complete sequencing/CNV interpretation service. [CPIC 2026 update, Table 1](https://files.cpicpgx.org/data/guideline/publication/ondansetron/2026/41979467.pdf)

CYP2C19 uses its own allele-function scheme. *1/*17 is rapid, *17/*17 ultrarapid, and *2/*17 intermediate; decreased-function combinations retain the appropriate “likely” category. Applying CYP2D6 cutoffs to CYP2C19 would be incorrect. No drug-specific recommendations from the classification sources are converted into generic dosage advice. [CPIC CYP2C19 update](https://files.cpicpgx.org/data/guideline/publication/clopidogrel/2022/35034351.pdf)

### Parent-equivalent conservation and exact propagation

For depot $D$, parent $P$, metabolite $M$ and cumulative eliminated parent-equivalent amount $L$, in nmol,

$$\dot D=-k_aD,\qquad\dot P=Fk_aD-(CL_{other}+CL_{CYP})P/V_P,$$
$$\dot M=yCL_{CYP}P/V_P-CL_MM/V_M,$$
$$\dot L=(1-F)k_aD+[CL_{other}+(1-y)CL_{CYP}]P/V_P+CL_MM/V_M.$$

Their derivatives sum to zero; each administration adds 1000 nmol. Parent and metabolite molecular weights need not agree because the ledger counts parent-equivalent molar units. Two more states integrate their concentrations. The linear six-state system propagates by matrix exponentials between seven 24 h dose events; the 168 h observation precedes a new dose.

The assumed pathway is $CL_{CYP}=10f_mg/(1+I/K_i)$ L/h, with $CL_{other}=10(1-f_m)$ and nominal $f_m=0.8$. The phenotype factors PM/IM/NM/RM/UM = **0/0.35/1/1.4/1.8** are hypothetical substrate parameters, independent of activity-score arithmetic. Inhibition alters enzyme function without rewriting genotype.

### Executed findings and next evidence

Nine diplotypes under three inhibitor conditions give **27 scenarios and 18,171 trajectory rows**; an **84-point** scan varies pathway fraction and activity. Uninhibited CYP2D6 PM/IM/NM/UM parent AUC0–168 is **2700.647/1217.852/594.118/362.775 nM·h**, while active-metabolite AUC is **0/663.207/934.816/1032.604 nM·h**. Parent accumulation and failed bioactivation can therefore occur together; phenotype alone does not establish a generic direction of clinical toxicity or efficacy.

An exposure-normalization multiplier $AUC_{reference}/AUC_{scenario}$ follows from model linearity. If metabolite formation is zero, no finite multiplier matches a nonzero reference metabolite exposure, so the result is null. Even a finite multiplier does not imply equal efficacy or safety. Clinical doses and mortality probability remain null.

Industrial translation needs a real drug and active species, verified diplotype/CNV calling, substrate-specific genotype PK, interactions, disease effects and the relevant drug-specific guideline. Numerical matrix accuracy and correct phenotype labels do not validate the assumed pharmacokinetic factors.

[Parameters](projects/task18_pgx/outputs/config.json) · [Results](projects/task18_pgx/outputs/summary.json) · [Sources](projects/task18_pgx/sources.md) · [Full report](projects/task18_pgx/REPORT_EN.md)

![Task 18: source-backed categories with hypothetical parent/metabolite exposure](figures_advanced/fig18_pgx_cyp2d6_phenotype_kinetics.png)

## 19. Atomistic ternary-interface sampling

The real 5T35 BRD4–MZ1–VHL/Elongin B/C system completed a 150,838-particle CPU pilot: 0.2 ps unbiassed thermalization, 1.0 ps biased dynamics and 20 hills. Coordinates and bias values are finite and audited, but the short trajectory does not converge a PMF; alpha, DeltaDeltaG and productive degradation geometries remain null. [Task 19 report](projects/task19_metadynamics/REPORT_EN.md) · [figure](figures_advanced/fig19_ternary_metadynamics_pmf_landscape.png)

## 20. Three-dimensional pocket-conditioned generation

The 6YB7 Mpro pocket run completed 100 candidates over 20 NSGA-II generations, an equal-budget random-order control, 12 Vina subset docks and two X77 redocking controls. Geometry, QED, MPO and SA are model-ranking quantities; redocking RMSDs were 1.242 Å and 1.055 Å, but no score was converted to nM and no novelty, synthetic-feasibility or clinical-efficacy claim is made. [Task 20 report](projects/task20_denovo/REPORT_EN.md) · [figure](figures_advanced/fig20_denovo_pareto_lead_optimization.png)
