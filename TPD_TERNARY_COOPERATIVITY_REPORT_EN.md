# Task 2: TPD ternary cooperativity and hook-effect modeling

## Scope and reproducibility

This is an executable, synthetic biophysical model, not a validated predictor of a specific degrader or a clinical dosing tool. No measured affinity, protein structure, calibrated assay, or full PROTAC structure was supplied. The selected defaults are illustrative: E0=100 nM; T0=100 nM; Kd,E=10 nM; Kd,T=100 nM; Vmax=12 nM/h; Km=10 nM; endpoint=24 h; basal turnover=0 /h; seed=20260910. The requested protein levels are modeling assumptions, not universal intracellular abundances. All internal concentrations are nM (1 nM = 1e-9 M), time is hours, and conformer distances are angstroms. The dose sweep covers 1e-12–1e-3 M, including concentrations that should not be interpreted as feasible cellular exposures.

Run `python run_task2_tpd_ternary_cooperativity.py`; use `--help` for parameters. The single file regenerates all four 300 DPI figures, both reports, CSV/SDF records, validation, dependency versions, checksums and logging. Dependencies are NumPy, SciPy, Matplotlib and RDKit; no network calls or external templates are needed. Seeds and single-thread conformer embedding make repeated runs reproducible within the recorded environment. Different RDKit releases may produce different conformers. Software versions and file hashes are recorded in `manifest_task2.json`.

## Event-driven pharmacology

An occupancy-driven inhibitor generally needs sustained engagement to suppress target function. A degrader can instead recruit an E3 ligase and initiate a ubiquitination event, after which target removal can outlast an individual binding event. Structural studies show that induced contacts between the recruited proteins contribute to cooperative recognition [3]. Binding remains necessary, but equilibrium ternary abundance alone does not establish productive ubiquitination, proteasomal processing, permeability or selectivity. The present equations assume all ternary complex contributes equally to an effective degradation flux; catalytic PROTAC recycling and fixed E3 total are explicit assumptions.

## Mass-action derivation and the correct optimum

Lowercase p, e, t denote free concentrations. K_E and K_T are the binary dissociation constants. Alpha is a supplied thermodynamic parameter, not inferred from linker geometry. The two ternary association paths obey the same thermodynamic cycle.

$$EP=ep/K_E,\quad PT=tp/K_T,\quad C=EPT=\alpha ept/(K_EK_T).$$
$$P_0=p+EP+PT+C,\quad E_0=e+EP+C,\quad T_0=t+PT+C.$$
$$K_{D,ET}=K_T/\alpha,\qquad K_{D,TE}=K_E/\alpha.$$
$$q(p)=\frac{(K_E+p)(K_T+p)}{\alpha p},\quad Cq=(E_0-C)(T_0-C).$$
$$C(p)=\frac{2E_0T_0}{E_0+T_0+q+\sqrt{(E_0-T_0)^2+2q(E_0+T_0)+q^2}}.$$
$$P_0(p)=p+(E_0-C)\frac{p}{K_E+p}+(T_0-C)\frac{p}{K_T+p}+C.$$
$$p_* = \sqrt{K_EK_T},\quad
P_{0,*}=\sqrt{K_EK_T}+\frac{E_0\sqrt{K_T}+T_0\sqrt{K_E}}{\sqrt{K_E}+\sqrt{K_T}}.$$


This derivation is included explicitly so the script can be audited independently of published three-body frameworks [1,2]. At fixed free p, eliminating e and t gives a quadratic in C. The rationalized smaller root above avoids cancellation. The total-to-free mapping is monotonically increasing for this stable equilibrium binding network; therefore maximizing C in free p also locates its maximum in total P. Since q=(p+K_E+K_T+K_E K_T/p)/alpha, q is minimized at p=sqrt(K_E K_T). At that point the two binary occupancy fractions sum to one, so the coefficient of C in P0 cancels. The **optimal total P is independent of alpha**, although peak height and threshold-window width are not.

The requested expression sqrt((K_E+E0)(K_T+T0)) is retained in the exported table as `supplied_formula_p_nm`; it is not used as a general identity. Here it gives **148.324 nM**, versus the exact **131.623 nM**, with free optimum **31.6228 nM**. They coincide in the fully symmetric special case K_E=K_T and E0=T0. Asymmetric-protein test cases in the validation file establish why the distinction matters.

At p approaching zero, C ~ alpha E0 T0 p/(K_E K_T), and rises from zero. At sufficiently large p, C ~ alpha E0 T0/p, while EP and PT approach their respective protein totals. Thus every finite positive alpha in this **two-binary-arm model** eventually exhibits a hook. A high-alpha curve is a cooperative PROTAC scenario; it is not a mechanistically complete molecular-glue model. A glue with negligible binding to one free protein, or a pre-existing protein-protein complex, requires a different reaction network. Neither a universal hook for all molecular glues nor an obligatory hook inside every experimental dose range follows from these equations.

| alpha | Peak EPT (nM) | Optimal total P (nM) | 80% low (nM) | 80% high (nM) |
|---:|---:|---:|---:|---:|
| 0.01 | 0.570646 | 131.623 | 70.7272 | 236.327 |
| 1 | 29.0536 | 131.623 | 68.3461 | 275.137 |
| 10 | 66.1477 | 131.623 | 69.2893 | 439.603 |
| 100 | 87.6755 | 131.623 | 73.6266 | 1287.53 |

The 80% windows refer to each curve's own maximum, not a common efficacy threshold and not a therapeutic window. The heatmaps use the exact maximum over total ligand at each affinity pair (0.1–1000 nM), with four alpha panels and one shared logarithmic color scale. Numerically maximizing each representative curve independently verifies the analytic peak.

## Solver implementation and validation

The main sweep uses `scipy.optimize.root` with hybr, then Levenberg-Marquardt if necessary, in logarithms of the three free concentrations. Residuals are log reconstructed totals minus log specified totals, evaluated using log-sum-exp. This enforces positive species without mixing incompatible scales. Convergence flags are insufficient: mass balance, finiteness and stoichiometric bounds are checked. A separate scalar Brent solution of the eliminated quadratic is the fallback and cross-check. The ODE uses vectorized log-bisection of the same exact mass balances at each time step, with changing target total.

Validation: **PASS**; 120 randomized positive systems, zero-dose and zero-protein cases, four invalid-input checks, both equilibrium routes, peak checks with asymmetric proteins, and a tighter-tolerance ODE comparison. Maximum relative mass-balance error: 3.232e-11. Root-versus-Brent species error: 2.135e-11. ODE endpoint fraction difference: 2.181e-13. These checks support the tested parameter regime; they do not establish the absence of every possible numerical or biological error.

## Linker ensemble and entropy interpretation

The two explicit, capped **linker proxies** are `[CH3:1]COCCOCCOC[CH3:2]` and `[CH3:1]C#CC#CC#C[CH3:2]`. Mapped atoms 1 and 2 define attachment-point positions. These fragments do not contain thalidomide, a VHL binder or a POI warhead; consequently the output is an attachment-point distance, not a measured protein-to-protein distance. Real warheads and protein structures are required to define exit-vector orientations and steric compatibility. The implementation chooses a rigid alkynyl example; piperazine chemistry is not modeled.

Each design receives 100 ETKDGv3 samples with explicit hydrogens and MMFF94s minimization. RMS pruning is disabled to preserve the requested sample count; repeated or symmetry-related minima are retained. These are algorithmic samples, not 100 guaranteed distinct conformers or a Boltzmann population. Nonconverged force-field calculations are retried and then cause a hard failure. Every conformer and its energy are exported as SDF and CSV.

After iterative least-squares heavy-atom superposition, per-atom RMSF is sqrt(mean_k ||x_ki-mean_k x_ki||²). The summary is the root mean square of those RMSFs over non-anchor linker heavy atoms. Rg is mass-weighted over all explicit atoms. Distance SD is a separate statistic and is not substituted for atomic RMSF. The fixed 8–12 A window is illustrative and was not tuned to the observed distribution.

| Linker | n | Mean r (A) | SD r (A) | Linker RMSF (A) | Mean Rg (A) | Compatible fraction |
|---|---:|---:|---:|---:|---:|---:|
| Flexible PEG | 100 | 9.6715 | 0.9577 | 0.7760 | 3.3518 | 96.00% |
| Rigid alkynyl | 100 | 9.6050 | 0.0000 | 0.0000 | 3.2967 | 100.00% |

Rigidification can reduce the number of accessible unbound conformations, but can also lock an unfavorable geometry. A compatible-state probability p_comp would imply a geometric selection cost -RT ln(p_comp) only if sampled populations were thermodynamic and the compatibility criterion adequately represented binding. Those conditions are not met here, so the script does **not** turn histogram widths or sample fractions into binding entropy, alpha, or potency. For an independently measured alpha, the coupling free energy is ΔG_coop=-RT ln(alpha). Positive cooperativity reflects the net balance of conformational restriction, strain, solvation and intermolecular contacts; one distance distribution does not determine it [3,4].

## Synthetic Western blot and HiBiT digital twin

The requested rate law is implemented as dT_total/dt = -Vmax C/(Km+C). The symbol k_deg in a concentration-rate law must have concentration/time units; here it is explicitly named Vmax in nM/h, not a first-order constant. Optional basal turnover adds k_base(T_initial-T_total), representing constant synthesis plus first-order turnover. Its default is zero to reproduce the requested no-synthesis equation. Integration uses log remaining target to preserve positivity and re-equilibrates C as target is consumed. It assumes rapid binding relative to degradation, constant total intracellular PROTAC, immediate recycling, fixed E3 and no drug clearance. At sufficiently long times the no-synthesis model can exhaust target even at weak effective doses; all reported DC50/Dmax values are time-specific.

Dmax is the numerically optimized degradation maximum at 24 h over the requested concentration range. DC50 is the **ascending** concentration giving half of that maximum, not necessarily 50% absolute target loss. A second crossing on the hook side is reported separately; missing crossings are null/outside range. The optimal dynamic endpoint dose can differ from the initial equilibrium peak because T_total evolves.

| alpha | Dmax (%) | DC50 rising (nM) | Half-max hook-side (nM) | Optimal endpoint dose (nM) |
|---:|---:|---:|---:|---:|
| 0.01 | 14.4574 | 35.5407 | 423.919 | 129.829 |
| 1 | 99.8679 | 5.6381 | 3328.34 | 111.455 |
| 10 | 100 | 2.45349 | 32695.6 | 110.723 |
| 100 | 100 | 2.13609 | 326418 | 113.187 |

HiBiT signal is 500 + 100000 × target fraction RLU. Three synthetic replicates use Gaussian SD = 3% of expected RLU + 100 RLU. Background-corrected normalization uses a known synthetic vehicle signal; noisy values are not clipped. A Hill function is fitted only to the ascending limb for alpha=10, with Dmax fixed to the mechanistic maximum. It is an empirical, conditional fit, not a mechanistic global fit to a bell curve: fitted DC50=2.2442 nM, Hill slope=1.9615, and rising-limb RMSE=3.5418 percentage points. No confidence interval is claimed from three generated replicates. The full-dose curve is always the dynamic mass-action prediction.

The Western blot is explicitly labeled **SYNTHETIC**, with target-band darkness proportional to remaining target and a constant loading control; vehicle and high-dose recovery are shown. It contains no fabricated claim of an actual experiment. Geometry sampling and assay kinetics are separate sensitivity modules; no unsupported linker-to-alpha calibration connects them.

## Figures

![Figure 1](figures_task2/fig1_ternary_hook_effect_curves.png)

![Figure 2](figures_task2/fig2_cooperativity_alpha_heatmap.png)

![Figure 3](figures_task2/fig3_linker_conformational_histogram.png)

![Figure 4](figures_task2/fig4_synthetic_western_blot_hibit.png)

## References

1. [Douglass et al. (2013), A Comprehensive Mathematical Model for Three-Body Binding Equilibria](https://doi.org/10.1021/ja311795d).
2. [A suite of mathematical solutions to describe ternary complex formation (2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7650257/).
3. [Gadd et al. (2017), Structural basis of PROTAC cooperative recognition](https://doi.org/10.1038/nchembio.2329).
4. [Affinity and cooperativity modulate ternary complex formation (2023)](https://www.nature.com/articles/s41467-023-39904-5).
5. [RDKit ETKDG and embedding API](https://www.rdkit.org/docs/source/rdkit.Chem.rdDistGeom.html).
6. [RDKit MMFF conformer optimization API](https://www.rdkit.org/docs/source/rdkit.Chem.rdForceFieldHelpers.html).

