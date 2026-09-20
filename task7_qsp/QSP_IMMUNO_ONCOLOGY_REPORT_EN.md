# Task 7: Tumor–immune combination QSP technical report

**Executed 2026-09-20. Evidence: synthetic model calculations; no calibration, external validation, patient data, or experimental measurements.**

This module couples oral small-molecule and finite-infusion antibody PK to a Simeoni-type tumor growth/damage system and CTL recruitment/exhaustion. The executed study contains 50 paired virtual people under four regimens for 60 days. It is a reproducible research prototype, not an industrially qualified or clinically validated product. The requested 50 mg QD and 10 mg/kg Q2W regimens are hypothetical simulation inputs, not recommended doses of an identified medicine.

## 1. Sources and necessary scientific corrections

The Simeoni framework supplies a preclinical precedent for exponential-to-linear growth and delayed drug-induced tumor damage. Here the user-specified proliferating-mass growth equation is extended by direct immune lysis. This does not reproduce every structural choice, parameter, or dataset in the original work. The large-mass parameter `lambda1` is a linear growth rate in mg/day, not a carrying capacity. Each of three damage compartments has residence time `tau`; their total mean delay is `3*tau`, here 4.5 days. [Simeoni et al., 2004](https://pubmed.ncbi.nlm.nih.gov/14871843/)

Positive excess over Bliss indicates departure from a specified independent-action reference, not proof of true clinical synergy. Loewe requires equivalent single-agent doses at the combination effect, and may be unsupported when effects lie beyond either monotherapy curve. [Bliss statistical methodology](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0224137), [Loewe consistency study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5808155/)

The requested mass increase above 120% of baseline is retained as a **model progression endpoint**. RECIST 1.1 instead uses lesion diameter sums, nadir-based relative increase of at least 20% with at least 5 mm absolute increase for target-lesion progression, and additional rules including new lesions. Our mass endpoint is not RECIST PFS. [RECIST Working Group](https://recist.eortc.org/recist-1-1/)

## 2. Exposure and system equations

Time is in days, tumor mass and administered amounts in mg, concentrations in mg/L, clearance in L/day, and CTL in arbitrary normalized density units. Linear one-compartment oral PK is calculated by superposition:

\[
C_s(t)=\sum_{t_d\leq t}\frac{F D k_a}{V_s(k_a-k_e)}
\left(e^{-k_e(t-t_d)}-e^{-k_a(t-t_d)}\right),\qquad k_e=CL_s/V_s.
\]

A continuous limit handles equal absorption and elimination rates. Antibody infusions last 0.5 hours, or 1/48 day. During an infusion, its concentration contribution is `R/CL_m*(1-exp(-k_e*u))`, where `R=D_m*70/T_inf` mg/day; after infusion it decays exponentially. There are 60 oral doses at days 0–59 and five antibody infusions at days 0, 14, 28, 42, and 56. A dose at the final study instant is excluded. Integration is segmented at every dose and infusion end. Infusion is not approximated as an instantaneous bolus.

\[
G(w_0)=\frac{\lambda_0w_0}{[1+(\lambda_0w_0/\lambda_1)^\psi]^{1/\psi}},\qquad
\dot w_0=G(w_0)-k_2C_sw_0-k_{kill}Ew_0.
\]

\[
\dot w_1=k_2C_sw_0-w_1/\tau,\quad
\dot w_2=(w_1-w_2)/\tau,\quad
\dot w_3=(w_2-w_3)/\tau,\quad W=\sum_{j=0}^3w_j.
\]

Mass balance is `dW/dt=G-w3/tau-k_kill*E*w0`. Entry into the damaged pool does not immediately remove tumor mass. Immune killing directly removes viable mass. The antigen proxy is `J_ag=w3/tau`; direct immune lysis is deliberately not included in antigen release, an acknowledged structural simplification rather than an experimentally established exclusion.

\[
RO=\frac{C_m}{K_D+C_m},\qquad
\dot E=s_{basal}+\alpha\frac{J_{ag}}{K_{ag}+J_{ag}}
-[\mu_E+k_{exh}PDL1(1-RO)]E.
\]

Basal recruitment is `(mu_E+k_exh*PDL1)*E_baseline` and `E(0)=E_baseline`. Untreated CTL is therefore nonzero and at equilibrium, allowing antibody monotherapy to act on pre-existing immunity. Setting both recruitment and baseline CTL to zero would make antibody monotherapy degenerate; that would be a modeling artifact. PD-L1 is fixed at relative value 1. There is no dynamic PD-L1 regulation, antigen-presentation compartment, tumor antibody distribution, or target-mediated drug disposition. Immunogenic death and antigen-driven CTL recruitment remain hypotheses requiring measurement.

## 3. Explicit assumptions and cohort design

All values in [parameters_used.json](parameters_used.json) are synthetic, not literature-calibrated drug parameters. Initial tumor mass is 100 mg; intrinsic doubling-time mean is 12 days; asymptotic growth is 8 mg/day; transition exponent is 20; and per-compartment damage residence is 1.5 days. Oral F=0.7, ka=12/day, CL=15 L/day, V=40 L, and drug damage coefficient=0.009 L/(mg·day). Antibody CL=0.2 L/day, V=3 L, and occupancy KD=5 mg/L. Baseline CTL=1; antigen stimulation maximum=0.9 CTL/day; half-maximal death flux=4 mg/day; CTL loss=0.08/day; exhaustion=0.35/day; immune killing=0.015/(CTL·day).

Small-molecule clearance, intrinsic doubling time, and baseline CTL use independent lognormal variability with CVs 30%, 25%, and 40%. The distribution parameters are `sigma=sqrt(log(1+CV²))` and `mu=log(arithmetic_mean)-sigma²/2`, so the requested mean and CV describe the generating population. A finite 50-person sample is not forced to have exactly those moments. Seed=20260920. Antibody clearance and body weight have no additional variability, and cross-parameter correlations are absent.

Each of the same 50 people receives all four counterfactual regimens: 200 trajectories do not mean 200 independent patients. SEM is sample SD/sqrt(50), reflecting finite synthetic-cohort sampling only. Vehicle includes basal immune killing: it means untreated, not an immune-disabled growth curve. Neither parameter-estimation uncertainty nor structural uncertainty is represented in these bands.

## 4. Executed outputs

Sampling every 0.25 days produces 241 time points and 48,200 state rows. All numbers below come from [results_summary.json](results_summary.json). Parameters were not fitted to efficacy measurements.

| Arm | Day-60 mass, mean ± SEM (mg) | Model progression events /50 | Administrative censoring | Median model progression (day) |
|---|---:|---:|---:|---:|
| Vehicle | 347.151 ± 8.051 | 50 | 0 | 4.292 |
| Targeted | 140.959 ± 6.154 | 43 | 7 | 4.628 |
| Anti-PD-1 | 84.198 ± 10.115 | 28 | 22 | 7.698 |
| Combination | 7.744 ± 1.179 | 14 | 36 | Not reached |

The first threshold crossing remains an event even if tumor mass later regresses. Thus a small late tumor burden can coexist with early model progression. Delayed damage also permits early total-mass growth despite inhibition of viable cells. This behavior depends strongly on the artificial endpoint; it is not a prediction of clinical resistance or treatment failure.

For each person, `f=1-W_arm/W_vehicle`; the Bliss score is computed within-person before averaging. Day-60 mean excess over Bliss is **0.06825 ±0.00686 SEM**. A separate nominal-person 6×6 dose landscape spans oral doses 0/10/25/50/100/200 mg QD and antibody doses 0/0.1/0.3/1/3/10 mg/kg Q2W. Its largest EOB is **0.16002**, at 50 mg QD plus 1 mg/kg Q2W. This is the largest value on a finite, hypothetical lattice without toxicity constraints, not an optimal clinical dose.

Loewe calculation uses 78 nominal single-agent support simulations and genuine piecewise-linear inverse doses: `CI=dA/DA(f_combo)+dB/DB(f_combo)`. Monotonicity, effect range, and inverse uniqueness are checked. Only **9 of 25** interior cells have supported indices; **16** remain `NA` because at least one monotherapy cannot attain the combination effect within the simulated support. Single-agent axes have CI=1; the origin is undefined. No extrapolation is substituted for missing support. Finite interpolation error remains, and mechanistically different therapies need not satisfy Loewe dose-equivalence assumptions.

## 5. Survival calculations respecting pairing

The first mass crossing above 120% baseline is linearly interpolated from the 0.25-day observation grid; otherwise the observation is administratively censored at day 60. Kaplan–Meier handles events before censoring at tied times. There is no mortality or dropout model, so the plot represents freedom from **model mass progression**, not clinical PFS.

The log-rank score uses pooled risk sets and is tested by swapping labels within each virtual-person pair, with 4095 permutations. Monte Carlo p-values use `(1+extreme_count)/(4095+1)`. Both exploratory comparisons give **0.0002441**, the simulation's p-value resolution floor, not zero or clinical evidence. No multiple-comparison adjustment was applied.

| Combination vs | Restricted mean delay through day 60 | Paired-bootstrap 95% interval | Descriptive Cox HR | Person-cluster-robust 95% interval |
|---|---:|---:|---:|---:|
| Targeted | 31.493 days | [24.289, 38.862] | 0.1949 | [0.1155, 0.3290] |
| Anti-PD-1 | 15.073 days | [8.722, 22.371] | 0.4323 | [0.2863, 0.6528] |

With common administrative censoring, the restricted mean contrast equals `mean(min(T_combo,60)-min(T_mono,60))`. Its interval uses 4000 bootstrap resamples of complete person pairs. Cox estimation uses Breslow ties and a sandwich variance formed from per-person sums of score residuals, with a G/(G−1) cluster correction. Proportional hazards is not established; delayed immune activity may violate it. HR is descriptive, while restricted mean delay is a more direct finite-horizon summary. No-event arms or separation produce an explicitly non-estimable HR, never an artificial zero/infinity. Every interval conditions on the assumed model and distribution, excluding uncertainty in their validity.

## 6. Verification, figures, and reproducibility

SciPy RK45 uses rtol=2e−7, atol=1e−9, and maximum step 0.2 days, with dosing discontinuities segmented explicitly. Nominal-person refinement of all four arms reduces tolerances tenfold and halves the solver step and observation grid. Maximum day-60 relative mass difference is **1.07e−6**, and maximum event-time difference is **0.000590 days**. These establish numerical stability on the checked scenarios, not biological validity. Tests also cover mass balance, the three-stage Erlang delay, dose/AUC identities, finite infusion, lognormal moments, paired null behavior, unsupported Loewe inversion, HR boundary behavior, and stored hashes.

The four requested 300-DPI PNGs have been visually inspected for axes, labels, legends, and interpretation boundaries:

1. [Tumor trajectories and SEM](figures_task7/fig1_simeoni_tgi_monotherapy_vs_combo.png)
2. [CTL and delayed-death flux](figures_task7/fig2_ctl_immune_infiltration_dynamics.png)
3. [Bliss and Loewe dose landscapes](figures_task7/fig3_bliss_synergy_matrix_heatmap.png)
4. [Model-endpoint Kaplan–Meier and paired statistics](figures_task7/fig4_virtual_cohort_kaplan_meier_pfs.png)

See [README](README.md) for offline CLI commands. A new output directory is required to preserve prior evidence. [manifest.json](manifest.json) records the script and generated output SHA256 hashes; [run_log.json](run_log.json) records versions and actual timing. The first complete study took approximately eight seconds using Python 3.12.14, NumPy 2.4.6, SciPy 1.18.0, and Matplotlib 3.11.1. Numerical outputs are deterministic within the same software environment; timing logs and version records are not expected to match across machines.

## 7. Path toward industry use

The prototype can support sampling-time design, diagnosis of delayed-response endpoint bias, monotherapy-support planning, and identification of unestimable interaction regions. Real development requires compound-specific total/unbound PK, tumor exposure and occupancy, a tumor size-to-mass mapping, longitudinal CTL/PD-L1 measurements, and both single-agent and combination dose experiments. Calibration, identifiability assessment, and held-out prediction must precede claims of translational performance.

Subsequent model comparisons should address dynamic PD-L1, antigen release from immune lysis, reversible exhaustion, antibody TMDD, dose sequencing, toxicity limits, and propagation of structural/parameter uncertainty. Combining a 100-mg tumor scale with a 70-kg dose conversion is an illustrative construction, not validated xenograft-to-human translation. Current results support no claim of true clinical synergy, PFS/OS benefit, or a recommended combination regimen.
