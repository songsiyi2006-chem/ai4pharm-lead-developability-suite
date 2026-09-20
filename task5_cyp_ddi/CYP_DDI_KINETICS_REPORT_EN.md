# CYP DDI kinetics, mechanism-based inactivation, and enzyme recovery

This module implements six named scenarios, three CYP isoforms, 28 BID administrations over fourteen treatment days, two oral probes, and ten days of recovery after the last perpetrator dose. It directly imports and calls the repository's actual Task 4 PBPK derivative. The principal compound PK and inhibition inputs are **constructed research assumptions**, so the outputs are not clinical predictions for the named drugs and do not establish safety, contraindications, or patient washout periods.

## 1. Evidence and completed computations

| Evidence level | Included work |
|---|---|
| Primary-source/regulatory context | FDA/ICH M12 inhibition models, FDA CYP examples, experiments on heme and apoprotein modification |
| Published model estimates | Quinney et al., Table 2: clarithromycin KI 5.3 µM; hepatic kinact 0.4 h⁻¹; intestinal kinact 4 h⁻¹ |
| User-specified inputs | Midazolam 2 mg and metoprolol 50 mg orally; three isoforms and fourteen-day BID protocol; illustrative hepatic CYP3A4 kdeg 0.019 h⁻¹ |
| Assumptions | Six complete PK/Ki/KI/kinact scenarios, other turnover constants, all tissue Kp=1, intestinal concentration proxy, probe clearance fractions |
| Executed calculations | 24 matched dynamic AUCR comparisons, 18 isoform screens, 36 recovery entries, three turnover scenarios, synthetic kinetic-parameter recovery |
| Missing evidence | Measured binding and inhibition, clinical PK/DDI calibration, identifiability, external validation, population and exposure–response uncertainty |

The three literature numbers are the authors' model estimates, not newly obtained in-vitro measurements. `literature_parameter_response.csv` contains 246 calculated transfer-function rows: two sites, three assumed constant concentrations, and 41 time points. It does not reinterpret the source KI as a measured unbound constant or insert these estimates into the principal PBPK scenarios. It is not a replication of the paper's clinical fitting or experimental curves. [Quinney et al., 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2812061/)

The six scenarios are Ketoconazole, Clarithromycin, Ritonavir, Fluconazole, Quinidine, and Amoxicillin. The first three illustrate reversible, time-dependent, and combined inhibition; fluconazole and quinidine illustrate isoform preferences. The amoxicillin scenario uses zero kinact and an intentionally huge Ki as a near-null control; it does not exclude real interactions involving that drug. Qualitative roles use the [FDA CYP examples](https://www.fda.gov/drugs/drug-interactions-labeling/drug-development-and-drug-interactions-table-substrates-inhibitors-and-inducers), without deriving any numerical kinetic or PK parameter from those roles. Metoprolol is the user-requested CYP2D6 probe scenario; it is not asserted to be the guideline's selected sensitive index substrate.

## 2. Chemical mechanisms and what the ODE can identify

Reversible inhibition follows available inhibitor concentration, whereas time-dependent activity loss can involve metabolic activation, tight coordination, or covalent modification. Heme-iron ligation, porphyrin destruction, and covalent apoprotein adduction are chemically different events. A declining activity curve alone cannot distinguish them: kinact in this model is an effective loss rate, not an atomistic mechanism assignment.

One reconstituted-system ritonavir study reported heme loss/heme-related modification; another used radiolabel and mass spectrometry to support modification of CYP3A4 Lys257. Their assay contexts and observations should be retained separately rather than promoting one mechanism to a universal human explanation. [Lin et al., 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3781371/), [Rock et al., 2014](https://pubmed.ncbi.nlm.nih.gov/25274602/)

Evidence for MBI requires appropriately controlled NADPH, time/concentration, drug-removal/recovery, and metabolite or protein/heme studies. Enzyme instability, inhibitor depletion, and nonspecific binding can confound interpretation. Figure 1 generates noise-free exponential activity curves from known parameters and refits them; this checks numerical parameter recovery, not experimental identifiability under realistic noise.

## 3. Static-model correction and regulatory interpretation

The prompt's proposed net-AUCR expression can decrease when its inhibition ratios increase. That direction is inconsistent with the inhibition-only model implemented here, so it was corrected. Basic screening ratios remain distinct from predicted exposure ratios:

\[
R_1=1+C_{max,u}/K_i,\qquad
R_2=1+\frac{k_{inact}\,5C_{max,u}}{k_{deg}(K_I+5C_{max,u})}.
\]

The implementation uses M12's factor 5. The 1.02 and 1.25 basic hepatic cutoffs and oral CYP3A `(Dose/0.250 L)/Ki` cutoff of 10 indicate when further evaluation is needed. They are not AUCR values. Here, Cmax is the simulated day-14 cycle peak sampled every 0.25 h, without a claim that periodic steady state has been established. The requested three-isoform subset does not replace a complete regulatory panel. [FDA/ICH M12, §2.1.2](https://www.fda.gov/media/161199/download)

For inhibition without induction, define residual activity `a=[kdeg/(kdeg+kobs)]/(1+I/Ki)`. The separately reported simplified M12 limit is:

\[
AUCR_{M12}=\frac{1}{a_g(1-F_g)+F_g}\frac{1}{a_hf_m+(1-f_m)}.
\]

Here fm denotes the affected hepatic-clearance fraction under the relevant clearance assumptions. `generalized_constant_mean_static_aucr` instead derives from this module's well-stirred liver, intestinal first pass, and renal clearance. With `s_h=1−Σw_i+Σw_i a_hi`, the weights refer to **intrinsic hepatic clearance**, not measured systemic fm:

\[
x=f_{u,p}CL_{int}s_h/R_b,\quad F_h=Q_h/(Q_h+x),\quad
CL_{h,p}=R_bQ_hx/(Q_h+x),
\]
\[
F'_g=[1+(1/F_g-1)s_g]^{-1},\quad
AUC=Dose\,F_aF'_gF_h/(CL_{h,p}+CL_{renal}).
\]

A test independently compares this expression with the integrated constant-inhibition ODE. AUCR returns to 1 without inhibition; the fully blocked simplified limit is 20 for fm=0.9 and Fg=0.5. This also detects the sign/direction problem in the prompt.

The outputs call the 1.25/2/5 bands `weak_range`, `moderate_range`, and `strong_range`. These names describe an illustrative model value. Clinical contraindications require the affected drug's exposure–response relationship, therapeutic window, and product-specific evidence; neither a green cell nor AUCR≥5 determines that decision here.

## 4. Dynamic coupling, units, and mass balance

Perpetrator PK uses the actual eleven-state Task 4 implementation: seven amount states (blood, liver, GI luminal depot, kidney, brain, lung, rest) plus hepatic elimination, renal elimination, loss before systemic entry, and AUC integrals. Pulmonary and systemic circulation are in series. Six additional dimensionless enzyme states cover hepatic/intestinal CYP3A4, CYP2D6, and CYP2C9, each normalized to baseline 1:

\[
\dot E_i=k_{deg,i}(1-E_i)-k_{inact,i}\frac{I_u(t)}{K_{I,i}+I_u(t)}E_i.
\]

The probes are simulated independently on the same perpetrator/enzyme trajectory, without probe–probe interactions; this is not a validated clinical cocktail. The probe derivative calls Task 4 first and then modifies hepatic metabolism and intestinal first-pass fluxes together with their matching loss integrals. Task 4's `FEC` bookkeeping becomes combined unabsorbed plus intestinal presystemic loss, rather than exclusively fecal excretion.

Amounts/doses are mg, time is h, concentrations are mg/L, and `µM=mg/L×1000/MW(g/mol)`. Unbound hepatic concentration derives from hepatic amount, Kp, and fu,p. The intestinal proxy is `I_gut,u=I_plasma,u+fu_gut Fa ka A_depot/Qgut`, converted to µM. Task 4 has no perfused enterocyte compartment, and none is falsely claimed here. All tissue Kp overrides are 1 as explicit assumptions.

The dose times are 0,12,…,324 h, with observation ending at 564 h. Each probe is given once at either 0 h or the start of treatment day 14, 312 h. AUC comparisons use identical control/inhibited windows of 564 h and 252 h respectively. The day-14 probe is followed 12 h later by the final perpetrator dose; residual inhibitor PK and enzyme inactivation continue afterward.

## 5. Executed results

Every value below is conditional on the constructed inputs. The static calculation freezes the grid-average hepatic/intestinal concentration over the pre-final-dose day-14 cycle and equilibrates enzyme abundance at that constant concentration.

| Scenario | Day-14 midazolam dynamic AUCR | Generalized static AUCR | Metoprolol dynamic AUCR |
|---|---:|---:|---:|
| Ketoconazole | 11.6746 | 11.7784 | 1.0555 |
| Clarithromycin | 6.4009 | 6.5416 | 1.0094 |
| Ritonavir | 16.9622 | 17.1904 | 1.1569 |
| Fluconazole | 2.1647 | 1.5792 | 1.0224 |
| Quinidine | 1.4253 | 1.1207 | 6.2996 |
| Amoxicillin near-null control | 1.0000 | 1.0000 | 1.0000 |

Dynamic and mean-concentration static results need not agree. Intestinal peaks, probe dosing phase, enzyme history, and inhibitor decline after withdrawal all affect the comparison. The fluconazole-scenario difference of 2.1647 versus 1.5792 is a timing/model result, not a correction to or validation against that drug's clinical DDI data.

Recovery is the last upward crossing of 90% after the last dose, with every later 0.1 h grid sample remaining ≥90%; the final crossing bracket is root-refined. In the clarithromycin scenario, hepatic CYP3A4 recovers in 4.9976 days and intestinal CYP3A4 in 3.3666 days. Both ritonavir CYP3A4 pools remain below 90% at ten days and are right-censored, not assigned a ten-day recovery. Reversible-only scenarios have zero enzyme-abundance recovery time, which does **not** mean absent residual drug or complete functional-clearance recovery.

`ln(2)/kdeg` is the enzyme-turnover half-life without continuing inactivation, not the time to 90% recovery. Halving, retaining, and doubling clarithromycin-scenario turnover gives midazolam AUCR 8.9494, 6.4009, and 4.4326, with hepatic recovery respectively not reached by ten days, 4.9976 days, and 2.5349 days. Turnover uncertainty can therefore change a threshold band. These are not patient washout recommendations.

![Dynamic probe exposure](results/figures_task5/fig2_dynamic_ddi_midazolam_pk.png)

![Recovery with residual inhibitor PK](results/figures_task5/fig4_cyp3a4_resynthesis_timeline.png)

## 6. Lead optimization and falsifiable experiments

Furan bioactivation can generate reactive intermediates; direct work on bergamottin identifies CYP3A4 modification associated with such activation. [Primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3336797/) Candidate changes such as furan replacement with a saturated oxygen heterocycle or another isostere, aniline replacement with a less oxidation-prone heteroaryl/amide, and removal or replacement of a methylenedioxyphenyl group are **SAR hypotheses**. They can also change target potency, pKa, solubility, clearance, and alternate metabolism. Their benefit cannot be assumed, and a structural alert alone is not toxicity evidence.

A useful next experiment is a matched Ki/NADPH-preincubation matrix with measured unbound fraction, inhibitor depletion, and substrate dependence. Positive TDI requires drug-removal/recovery experiments, metabolite trapping, and protein/heme mass spectrometry. Hepatocyte induction, major circulating metabolites, real oral PK, and external victim-DDI data would then support calibration and validation. This module omits induction, perpetrator self-inhibition, inhibitory metabolites, transporters, and genetic phenotypes; particularly complex perpetrators cannot be quantitatively extrapolated from this demonstration.

## 7. Reproduction and numerical checks

See the [module README](README.md) for commands and artifact links. `--config` ingests the fixed six-name JSON panel and exports all resolved and inherited Task 4 settings. Unknown fields are rejected and assumption labels retained. `--out` requires a new directory, while `--self-test` does not mutate frozen outputs. Execution is offline and performs no package installation or Git operation.

Eleven tests cover concentration conversion, invalid kinetics, the regulatory factor 5, no/full inhibition limits, the constant-exposure enzyme solution, the unchanged Task 4 limit, mass conservation, analytical PBPK AUC, recovery/censoring, and category boundaries. A tighter coupled solve changes clarithromycin–midazolam AUC by 8.60×10⁻⁸ relatively. Maximum perpetrator mass error is 2.58×10⁻⁷ mg, maximum probe mass error 3.11×10⁻⁸ mg, and maximum terminal probe fraction 4.94×10⁻¹⁵. These establish numerical behavior, not clinical validity. Four 300-DPI PNGs were visually inspected; numerical artifacts have SHA256 hashes, and the unchanged Task 4 source hash is recorded in the summary.
