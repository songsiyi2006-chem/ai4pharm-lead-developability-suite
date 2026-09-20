# Task 3 — Covalent inhibitor kinetics, residence and GSH reactivity

## Scope and evidence status

This is a reproducible computational demonstration, not a calibrated drug discovery
prediction. Eight hypothetical warheads share an arylaminopyrimidine core; no target
protein, binding pose or measured kinetic dataset was supplied. They are not
sotorasib, osimertinib or ibrutinib. Orbital descriptors are calculated; default
binding and chemical rates are illustrative. The barrier model is an uncalibrated
descriptor surrogate, not a located transition state. No clinical efficacy,
therapeutic window, safety or DILI probability is established by these outputs.

## 3A — Electronic descriptors and thiolate barrier scenario

RDKit validates structures, tracks the reacting atom, samples ETKDGv3 conformers,
and chooses the lowest-energy converged MMFF conformer. EHT/YAeHMOP is a genuine
semiempirical orbital calculation, but its absolute energies are not DFT energies
and must not be mixed with xTB values in a calibrated correlation. The optional
xTB mode optimizes with GFN2-xTB/ALPB water and computes neutral/anion single points
at the same neutral geometry. Logs, geometries and JSON are retained.

Using the Koopmans-style approximation I≈−EHOMO and A≈−ELUMO:

$$\mu=(E_H+E_L)/2,\quad\eta=(E_L-E_H)/2,\quad
\omega=\mu^2/(2\eta),\quad S=1/(2\eta).$$

The hardness/softness convention is explicitly fixed above. EHT local f+ is the
atom-resolved LUMO population (a frozen-orbital proxy); xTB f+ is q(N)−q(N+1), a
vertical finite difference. Local softness is S f+. Neither is a direct validation
of a reaction pathway. The beta carbon is compared with alpha, and its rank among
all heavy atoms is exported. Negative population values are retained, not clipped
to fabricate localization. Chloroacetamide undergoes SN2 at the carbon bearing Cl;
calling that atom a Michael beta carbon would be incorrect. Neutral microstates,
one selected conformer and unverified anion binding limit these descriptors;
protonation of the dimethylamino group requires separate pH-dependent modeling.

The optional barrier scenario is fully specified:

$$\Delta G^\ddagger_{scenario}=B_{warhead}+
\mathrm{clip}[2\,(E_L-\overline{E_L}),-4,4]\quad\mathrm{kJ/mol}.$$

The coefficient is 2 kJ mol⁻¹ eV⁻¹. The eight baselines are declared in the parameter
JSON, not learned from data. For a 1 M standard concentration and transmission
coefficient one:

$$k_{thiolate}=(k_BT/h)/C^\circ\exp[-\Delta G^\ddagger/(RT)],\quad
f_{thiolate}=1/(1+10^{pK_a-pH}),\quad k_{GSH}=f_{thiolate}k_{thiolate}.$$

An explicit kGSH override bypasses this Eyring scenario. The retained barrier column
then remains a separate scenario and is not inferred from that measured rate.
Real barrier validation requires a solvent/microstate-consistent thiolate addition
or SN2 pathway, transition-state optimization, one appropriate imaginary frequency,
IRC confirmation, free-energy corrections and experimental calibration.

## 3B — Two-step inactivation and regression

For free target F, noncovalent complex B and covalent complex C:

$$E+I\rightleftharpoons B\rightleftharpoons C,$$
$$\dot F=-k_{on}IF+k_{off}B,$$
$$\dot B=k_{on}IF-(k_{off}+k_{inact})B+k_{rev}C,$$
$$\dot C=k_{inact}B-k_{rev}C.$$

For the irreversible limit krev=0, rapid pre-equilibrium gives
B/(F+B)=I/(KD+I), where KD=koff/kon. Thus for unmodified target U=F+B,
Udot=−kinact I U/(KD+I), yielding the requested saturating kobs equation.
This reduction needs constant free inhibitor, negligible depletion, and binding
equilibration faster than chemical inactivation. KD is not generally the fitted KI.
The quasi-steady-state kinetic constant is KI,QSSA=(koff+kinact)/kon. The exact slow
decay eigenvalue at clamped I is:

$$a=k_{on}I,\quad S_r=a+k_{off}+k_{inact},\quad
\lambda_{slow}=\frac{2ak_{inact}}{S_r+\sqrt{S_r^2-4ak_{inact}}}.$$

Its low-I slope is kinact/KI,QSSA. At high I it tends to kinact. An exact transient
is biexponential; a hyperbola is not generally its exact slow eigenvalue.

The code numerically solves full mass action at 1 nM–100 µM. After 12 fast-mode time
constants, nonlinear least squares fits A exp(−kobs t) to U. It then fits a positive
hyperbola in log-rate residuals, giving equal relative weighting across doses.
The CSV distinguishes KD, KI,QSSA, fitted KI, fitted kinact, efficiency definitions
and approximation error. These are noiseless synthetic trajectories, so statistical
confidence intervals would not represent experimental uncertainty and are omitted.
Acquisition times adapt to slow rates and can exceed realistic assay durations.
Experiments require fixed observation windows, replicates and identifiability checks.

For comparable forward reactivity, figure 1 sets krev=0 for all warheads; cyano
compounds' reversibility is restored in dynamic exposure and washout. The forward
efficiency does not quantify their equilibrium occupancy or residence. Instantaneous
inhibition is B+C, whereas U=F+B is the covalently unmodified fraction.

A separate finite-bolus simulation includes actual inhibitor depletion, clearance,
GSH trapping and reverse GSH release. With target normalized by Etotal,
I+ B+C+SG+cleared is conserved. This assay has no protein turnover; the washout
experiment below handles turnover separately. No substrate competition is included.
The requested efficiency categories (<10³, 10³–10⁴, 10⁴–10⁶, >10⁶ M⁻¹s⁻¹) are
screening labels, not universal thresholds for clinical optimality.

## 3C — Residence, target turnover and rapid washout

Affinity contributes to recognition and exposure-dependent occupancy; residence
alone does not dictate clinical efficacy. Pharmacodynamics also depend on unbound
exposure, target synthesis, target vulnerability, tissue distribution and selectivity.
Long-lived covalent occupancy can persist after free drug clearance, as demonstrated
experimentally for reversible covalent kinase inhibitors [Bradshaw et al., 2015].

All target states degrade at kdeg=ln2/t1/2,protein; synthesis supplies free target
at ksyn=kdeg Etotal,baseline. Accordingly:

$$\dot F=k_{syn}-k_{deg}F-k_{on}IF+k_{off}B,$$
$$\dot B=k_{on}IF-(k_{off}+k_{inact}+k_{deg})B+k_{rev}C,$$
$$\dot C=k_{inact}B-(k_{rev}+k_{deg})C.$$

The cells receive a 2-hour constant free-drug pulse, followed by an instantaneous
perfect sink for free inhibitor. Bound drug remains; released drug is immediately
removed. The experiment has no residual intracellular depot or external rebinding.
Reversible C→B can still re-form C within the bound complex.

The requested reduced irreversible model is also integrated:
Fdot=ksyn−(kdeg+kobs)F. After washout,
F(t)=Etotal−[Etotal−F(0)]exp(−kdeg t).
For a single reversible noncovalent state intrinsic residence is 1/koff, while the
bound cellular lifetime is 1/(koff+kdeg). An irreversible bond has infinite chemical
lifetime in this model; 1/kdeg is the mean protein/adduct lifetime, not a chemical
off-rate, and ln2/kdeg is the recovery half-life. For reversible covalent states,
1/krev is just the bond-opening lifetime and differs from complete target release.

The full bound-state matrix Q is used to compute the occupancy-weighted mean
first-passage lifetime after washout: tau=1ᵀ(−Q)⁻¹p_bound(0). This accounts for
bond reversal, re-covalentization, dissociation and turnover. The numerical recovery
half-time is measured relative to inhibition at washout; it is not the same as tau.
An unreached half-time is exported as null rather than guessed or extrapolated.

## 3D — GSH twin and screening trade-off

At buffered 5 mM GSH (configurable), the irreversible pseudo-first-order model is
I(t)=I0 exp(−kGSH[GSH]t), so t1/2=ln2/(kGSH[GSH]). The code numerically integrates
this twin and checks against its analytic solution. For reversible GSH adducts,
the adduct fraction is p[1−exp(−(a+b)t)], with a=kGSH[GSH], b=kreverse,GSH,
p=a/(a+b). It exports the forward half-life, relaxation half-life, equilibrium
fraction and actual time to 50% adduct. The latter may not exist if p≤0.5.

The requested safety ratio (kinact/KD)/kGSH is dimensionless; a QSSA ratio is also
provided. It compares two rate efficiencies, not concentrations, selectivity across
proteins, exposure or clinical safety. Figure 3 maximizes both on-target efficiency
and forward GSH half-life, marks the Pareto frontier, and shades the requested
10⁴–10⁶ M⁻¹s⁻¹ / ≥30 min screening region. Its historical filename contains 'radar',
but it is the requested Pareto scatter plot.

A forward GSH half-life below 30 min triggers a hyperreactivity screening flag.
This is a user-specified heuristic, not a validated idiosyncratic DILI classifier.
GSH conjugation can detoxify electrophiles; GSH depletion and protein modification
depend on exposure, regeneration, metabolism and compartment. A buffered-GSH
assay does not itself simulate cellular GSH depletion. GSH reactivity measurements
provide intrinsic-reactivity information [Flanagan et al., 2014], not clinical outcomes.

## Reproducibility and limits

No random noise is added to kinetic data. The seed controls conformer generation;
library versions and SHA-256 hashes are recorded. Numerical checks cover mass
balance, positivity, ODE/analytic agreement, limiting kinetics, turnover recovery
and PNG DPI. Determinism is expected within the same software build, not necessarily
bitwise across platforms. EHT and xTB remain approximate; no reaction barrier,
binding potency or rate is claimed as experimentally validated. Figure 4 is an
exploratory OLS relationship across mixed mechanisms, with independently specified
kinetic inputs; neither a strong nor weak correlation establishes causality.

From the repository root, use `python projects/task03_covalent_kinetics/run_task3_covalent_kinetics_residence_time.py --help` for all options.
Run with `--output-dir work/task03_reproduction` to keep fresh outputs separate from this archived project.
Edit a copy of [data_task3/parameters_template.json](data_task3/parameters_template.json) and pass `--parameters` to supply rates
with provenance. Partial overrides retain the other illustrative values; provenance
must describe this. GSH rate overrides are apparent rates at the configured assay
conditions and are not automatically pH-corrected. `--self-test-only` checks core
mathematics. `--git-sync` requires an existing main checkout and origin remote;
generated paths alone are staged and pushed, without force or branch switching.
The script can be rerun: its README section is replaced idempotently and artifacts
are overwritten only at their dedicated paths. Avoid concurrent runs to one folder.

## Results of this run

- Electronic method: EHT/YAeHMOP; frozen-orbital LUMO population proxy.
- Temperature: 310.15 K; pH 7.4; GSH pKa assumption 8.7; GSH 5 mM.
- Pareto candidates: W02, W06, W08.
- Hyperreactivity flags: W03, W04, W05.
- Beta localization unsupported or weak: W05, W06.
- LFER R²: 0.2992 (illustrative).
- Maximum ODE kobs / eigenvalue relative error: 7.81e-08.

| ID | Warhead | LUMO (eV) | kinact/KD (M^-1 s^-1) | GSH forward t1/2 (min) | Ratio | Flag <30 min |
|---|---|---:|---:|---:|---:|---|
| W01 | Acrylamide | -9.479 | 3e+04 | 68.73 | 8.92e+05 | no |
| W02 | Dimethylaminomethyl acrylamide | -9.364 | 2e+04 | 240.32 | 2.08e+06 | no |
| W03 | Alpha-cyanoacrylamide | -9.780 | 6e+04 | 3.60 | 9.36e+04 | yes |
| W04 | Chloroacetamide | -9.329 | 4.5e+04 | 3.47 | 6.76e+04 | yes |
| W05 | Vinyl sulfone | -9.330 | 1.07e+04 | 11.09 | 5.12e+04 | yes |
| W06 | Methacrylamide | -9.329 | 400 | 790.61 | 1.37e+05 | no |
| W07 | Beta-methyl cyanoacrylamide | -9.583 | 9e+04 | 43.02 | 1.68e+06 | no |
| W08 | Beta-isopropyl cyanoacrylamide | -9.509 | 1.2e+05 | 145.82 | 7.57e+06 | no |

### Input provenance by compound

- W01: kinetics = Illustrative scenario; no experimental calibration; GSH = Eyring/LFER scenario; uncalibrated
- W02: kinetics = Illustrative scenario; no experimental calibration; GSH = Eyring/LFER scenario; uncalibrated
- W03: kinetics = Illustrative scenario; no experimental calibration; GSH = Eyring/LFER scenario; uncalibrated
- W04: kinetics = Illustrative scenario; no experimental calibration; GSH = Eyring/LFER scenario; uncalibrated
- W05: kinetics = Illustrative scenario; no experimental calibration; GSH = Eyring/LFER scenario; uncalibrated
- W06: kinetics = Illustrative scenario; no experimental calibration; GSH = Eyring/LFER scenario; uncalibrated
- W07: kinetics = Illustrative scenario; no experimental calibration; GSH = Eyring/LFER scenario; uncalibrated
- W08: kinetics = Illustrative scenario; no experimental calibration; GSH = Eyring/LFER scenario; uncalibrated

## References

- [Flanagan et al. (2014), experimental GSH reactivity and computational methods](https://doi.org/10.1021/jm501412a)
- [Bradshaw et al. (2015), reversible covalent inhibitors and tunable residence](https://doi.org/10.1038/nchembio.1817)
- [RDKit: rdEHTTools / YAeHMOP API](https://www.rdkit.org/docs/source/rdkit.Chem.rdEHTTools.html)
- [xTB: orbital properties and machine-readable JSON output](https://xtb-docs.readthedocs.io/en/latest/properties.html)
- [Assay Guidance Manual: mechanism-of-action assays](https://www.ncbi.nlm.nih.gov/books/NBK92001/)
