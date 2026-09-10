# Task 4 — PBPK and exposure-based dose screening

## Executive summary
Compound: **SYNTHETIC_LEAD_DEMO**. Synthetic inputs: **True**.
This is a reproducible, reduced seven-state research PBPK model, not a qualified clinical prediction platform. Input provenance, assumptions, numerical checks and uncertainty must be reviewed before interpreting exposure. No lead-specific measured clearance, tissue Kp values, toxicology dataset or human PK observations was supplied for the delivered demonstration. The configured toxicity threshold is a scenario assumption, not an established safe exposure limit. Setting `synthetic_inputs=false` changes the label only; it does not qualify the model.

The selected **research scenario** is: 15 mg BID (30 mg/day), coverage 90.443%, Cmax 0.2802 mg/L. A clinical recommended dose is deliberately left null in the machine-readable results. QD is every 24 h and BID every 12 h. Doses refer to active compound mass; salt correction is not inferred.

## Translational rationale and model scope
PBPK joins compound properties, organ capacities, circulation and clearance to connect an administered dose with time-varying exposure. Comparing free plasma exposure with potency can rank formulations or regimens and identify measurements that most affect a development decision. Here, IC90 is treated as a free-concentration threshold. If an assay reports a nominal concentration, binding and assay-to-tissue translation must be resolved separately. Time above IC90 is an exposure surrogate; it is not a measurement of receptor occupancy or clinical efficacy.

The requested seven compartments are blood/plasma, liver, GI depot, kidney, brain, lung and remaining tissues. Blood is an amount in whole blood with a mixed-venous plasma observation obtained using Rb. GI is a lumen depot; gut wall, portal blood and arterial blood do not have separate states. Lung is in series with the systemic organs. Hepatic flow lumps hepatic arterial and portal supply. The renal sink is placed in the central state as specified by the brief; the kidney tissue state is distribution-only. This is a reduced circulation approximation. Remaining tissues use an effective composition and act as a target proxy, not a separately validated target organ.

## Inputs and physiology
All full inputs and units are in `parameters_used.json`; an editable template is `compound_example.json`. MW=450 g/mol, cLogP=3, ionization=base, acidic/basic pKa=None/7.8, fu,p=0.1, Rb=1, CLint,mic=8 uL/min/mg, fu,mic=1, ka=1 h^-1, Fa=0.9, Fg=1. MW is used for unit conversion; it does not independently determine Kp in this model. Blood volume=5 L. Tissue volumes and flows are fixed adult scenario values, not individualized allometric predictions; 1,800 g is liver mass, distinct from the 1.8 L liver distribution volume.

| Tissue | Volume L | Flow L/h | Kp tissue:plasma | Fractions: neutral lipid, phospholipid, water |
|---|---:|---:|---:|---|
| liver | 1.8 | 90 | 2.74764 | 0.0138,0.0303,0.705 |
| kidney | 0.31 | 66 | 2.38782 | 0.0121,0.024,0.783 |
| brain | 1.4 | 45 | 6.25609 | 0.0391,0.0533,0.77 |
| lung | 0.5 | 360 | 0.925614 | 0.003,0.009,0.811 |
| rest | 60.99 | 159 | 6.13669 | 0.05,0.015,0.65 |

Tissue partitioning uses the Poulin–Theil composition form [1], with a clearly labeled Henderson–Hasselbalch logD7.4 substitution for ionization. Protein partitioning uses the configured interstitial/plasma ratio. This adaptation assumes equal pH, passive equilibration and binding surrogates. It is **not** a full Rodgers–Rowland or Berezhkovskiy implementation. It does not model acidic phospholipid binding, lysosomal trapping, active transport or BBB permeability. Tissue fractions are explicit illustrative assumptions. Brain and pooled-rest predictions particularly require empirical verification; measured Kp overrides are supported. Ampholytes use a simplified ionization approximation.

## IVIVE and equations
Microsomal clearance is interpreted as apparent clearance based on total incubation concentration and corrected by fu,mic. Set fu,mic=1 if the supplied clearance is already unbound. The explicit factor 60/1e6 converts uL/min to L/h. At the current inputs, CLint,u=38.88 L/h, CLH,b=3.72699 L/h and low-concentration FH=0.958589. FH is hepatic escape, not overall oral availability. Use it once through hepatic metabolism; multiplying the liver input by FH again would double-count first pass.

```text

Let M_B be mixed-venous blood drug mass; Cp=M_B/(V_B Rb).
For tissue i, C_i=M_i/V_i and C_vi=Rb*C_i/Kp_i (venous blood).
C_a=C_v,lung; Qc=sum_i Q_i, i=liver,kidney,brain,rest.

    dM_g/dt = -ka*M_g
    dM_lung/dt = Qc*(Rb*Cp - C_a)
    dM_i/dt = Qi*(C_a-C_vi)                      [kidney, brain, rest]
    dM_L/dt = QH*(C_a-C_vL) + ka*Fa*Fg*M_g - H
    dM_B/dt = sum_i(Qi*C_vi) - Qc*Rb*Cp - CLrenal,p*Cp
    H_linear = CLint,u * fu,p * C_L/Kp,L
    H_saturable = Vmax*Cu,L/(Km+Cu,L), Cu,L=fu,p*C_L/Kp,L
    Vmax = CLint,u*Km; CLrenal,p = GFR*fu,p

Bookkeeping: dE_H/dt=H; dE_R/dt=CLrenal,p*Cp;
dE_pre/dt=(1-Fa*Fg)*ka*M_g; dAUC/dt=Cp.
M_B+M_g+sum(M_tissues)+E_H+E_R+E_pre = cumulative administered dose.

    CLint,u [L/h] = (CLint,mic / fu,mic) * MPPGL * liver_mass_g * 60/1e6
    fu,b = fu,p/Rb
    CLH,b = QH*fu,b*CLint,u / (QH+fu,b*CLint,u)
    EH = CLH,b/QH; FH = 1-EH; F_linear = Fa*Fg*FH
    Cfree [nM] = Cp [mg/L] * fu,p * 1e6 / MW [g/mol]

Kp approximation (common pH 7.4):
    D = 10^cLogP / (1 + acid_term + base_term)
    acid_term = 10^(7.4-pKa_acid); base_term = 10^(pKa_base-7.4)
    [include only terms appropriate to the selected ionization class]
    fu,t = 1/(1 + protein_ratio*(1/fu,p-1))
    A(f) = D*(f_nl+0.3*f_ph)+f_water+0.7*f_ph
    Kp,t = A(tissue)/A(plasma) * fu,p/fu,t

Linear model dM/dt=A*M:
    AUC_0_inf = e_B^T*(-A)^(-1)*M0/(V_B*Rb)
    M_ss,post = (I-exp(A*tau))^(-1)*dose_vector

```

The default ODEs are linear. Optional `km_unbound_mg_l` enables saturable metabolism with Vmax=CLint,u*Km. Km must be supplied or explicitly assumed; it cannot be inferred from microsomal CLint alone. Well-stirred clearance then describes only the low-concentration limit. With nonlinear kinetics the equal-dose oral/IV AUC ratio remains an apparent exposure ratio and is not generally the absolute fraction reaching circulation; it can exceed 100%. No dose-linear scaling is used in nonlinear dose screening.

## Numerical methods and results
Radau or BDF integrates every dosing interval separately, with exact state jumps at administration. The 7-day simulation gives 14 oral doses at 0,12,...,156 h; 168 h is pre-dose. Seven amount states and four integral counters are distinguished. Matrix solutions independently check linear integration, integrate the entire AUC tail and solve the periodic steady-state boundary. Nonlinear AUC integrates until residual mass is below 1e-9 of the dose, then adds a low-concentration tail approximation. Finite limits trigger errors rather than silently accepting unconverged results. Reported half-lives are asymptotic low-concentration eigenmode rates; the oral rate allows flip-flop absorption and is not a fitted half-life from the 48 h graph. Cmax is searched on 0-48 h for the single doses.

| Endpoint | Value |
|---|---|
| IV Cmax / Tmax | 20 mg/L / 0 h |
| Oral Cmax / Tmax | 0.393085 mg/L / 0.389014 h |
| IV AUC 0-inf | 22.3364 mg h/L |
| Oral AUC 0-inf | 19.2703 mg h/L |
| Oral apparent terminal half-life | 62.5233 h |
| AUC-ratio oral F | 86.273% |
| Oral AUC after 48 h | 57.946342612108964% |
| Day 7 pre-dose trough | 1.25124 mg/L |
| True SS peak / pre-dose trough | 1.86809 / 1.48126 mg/L |
| True SS within-interval minimum | 1.48126 mg/L |
| True SS 24h coverage | 100% |
| Cmax accumulation ratio | 4.752378630997845 |
| Day 7 vs SS state relative error | 0.155287 |
| Day 7 reached SS (1% state/peak/trough criterion) | False |

The day-7 result is tested against a separate periodic steady state; seven days is not assumed sufficient. The pre-dose trough and actual within-interval minimum are both reported because continued absorption can move the minimum away from the dosing boundary. Accumulation uses the peak over the first 12 h interval, not an arbitrarily long single-dose window. Complete metrics, mass errors and periodic residuals are in `results_summary.json`.

## Dose rationale, uncertainty and reproducibility
For each QD/BID regimen, the script finds the minimum on the configured per-administration lattice (10–800 mg, step 1 mg), evaluates free-plasma coverage >= 90% and enforces **total-plasma** Cmax < 5 mg/L. The daily-dose objective can therefore span 10–800 mg/day for QD and 20–1600 mg/day for BID with default bounds. Root-refined threshold crossings determine continuous covered time, including strict equality handling. Regimen minima are compared by total daily dose, then peak. The plotted shaded region is illustrative on the displayed grid; the tabulated selected dose is separately verified on the search lattice.

`sensitivity_analysis.csv` changes fu,p, microsomal clearance, absorption, and all Kp values one at a time. These scenario changes are not confidence intervals, a population simulation or global sensitivity analysis. Key missing determinants include formulation/dissolution, intestinal metabolism evidence, transporter effects, organ disease, interindividual variability, toxicology and observed human PK. Independent clinical data would be needed to evaluate predictive performance and clinical dose selection. FDA guidance [3] distinguishes intended use and the evidence supporting a PBPK report; passing code tests alone does not establish clinical validity.

Run `python run_task4_pbpk_pharmacokinetics_dose_prediction.py --self-test --out results`. To replace assumptions, export a template with `--write-example compound.json`, edit it, then use `--config compound.json`. `verification.json` records numerical checks; `manifest.json` contains file hashes and versions. Four figures are saved at 300 DPI. The self-contained script regenerates this report and all artifacts. Task 1 descriptors can be transferred through the JSON template; they do not replace measured binding, clearance or toxicology inputs.

## References
1. Poulin & Theil (2002), tissue-composition partitioning: https://doi.org/10.1002/jps.10005
2. Poulin & Theil (2002), generic PBPK models: https://doi.org/10.1002/jps.10128
3. FDA (2018), PBPK report format and content: https://www.fda.gov/regulatory-information/search-fda-guidance-documents/physiologically-based-pharmacokinetic-analyses-format-and-content-guidance-industry
