# ADC translational engineering: a conservative multiscale computation

## Evidence and objective

This project executes the four requested numerical modules using publicly checked mechanistic literature and hypothetical parameters. It produces inspectable trajectories, mass ledgers, four 300-DPI figures and numerical checks. It does not claim a production-qualified assay, clinically optimal DAR, actual DXd/MMAF kinetic constants, or experimentally demonstrated tumor killing. Sources, access limitations and parameter origins are registered in [SOURCES.md](SOURCES.md) and [parameter_provenance.csv](parameter_provenance.csv).

Industrial relevance lies in connecting a manufacturing attribute—the distribution of attached payloads—to a separations model, cellular payload supply and spatial exposure. This connection can help decide which measurements to collect. It cannot replace release specifications, validated chromatography, clinical exposure-response work or product-specific safety assessment.

## 8A: partial reduction and finite-reagent conjugation

An IgG1 is represented by four interchain disulfide bonds. Reduction is a completed upstream stage with independent bond probability p = 1 − exp(−k_red t_red); the default probability is 0.950213. The number of reduced bonds is binomial, and a molecule with b reduced bonds has capacity m = 2b. The simulated state is the joint fraction f_(m,n), where m is 0, 2, 4, 6 or 8 and 0≤n≤m is the actual DAR. There are 25 joint states and one finite free-linker pool L, expressed as equivalents per initial antibody. Reduction reagent depletion, differential disulfide accessibility, intrachain reduction and subsequent oxidation are not modeled.

The n→n+1 flux is k₀(m−n) exp(−αn) L f_(m,n). Each event consumes one linker and moves one antibody between states. The penalty α=0.18 represents an assumed steric/accessibility effect. Summing over capacities yields f_DAR0…f_DAR8. Both total antibody fraction and L + Σn f_DARn are conserved. Odd DAR states remain because individual thiols react sequentially; grouping them into even peaks would destroy this model's probability and payload balance. The result is not claimed to reproduce the strong even-species patterns possible under different conjugation conditions.

At 8 h, input linker equivalents 2, 4, 6 and 10 give average DAR 2.000000, 3.999320, 5.912516 and 7.523958, respectively. Their upper capacity is set by partial reduction as well as finite reagent. For the 4-equivalent condition, 24 independent Gillespie simulations of 1,000 antibodies yield mean DAR 3.999500, Monte Carlo standard error 0.0001474 and an absolute difference from the master-equation limit of 0.0001803. The largest species-fraction discrepancy is 0.005599. Finite population and rounded integer reagent counts differ from the continuum limit; exact identity is not expected. Seed 82026 and individual replicate outputs are archived.

Site-specific conjugation is a distinct comparator, not the reaction simulated here. Junutula and colleagues demonstrated engineered-cysteine conjugation with controlled attachment in their tested constructs. Hamblett and colleagues found a loading-dependent tradeoff for cAC10–MMAE; the favored DAR in that preclinical context cannot establish a universal clinical optimum. Ogitani's studied DS-8201a construct had approximately DAR7–8, further illustrating why a universal DAR4 rule is inappropriate. See the primary papers in [SOURCES.md](SOURCES.md).

## 8B: empirical HIC twin

The default descending salt gradient is 1.5 M to zero over 35 min. Assumed hydrophobic retention follows t_R(n)=2+2n^1.15 minutes. This maps each DAR to its elution time and the corresponding salt concentration; it is an empirical retention rule, not a mechanistic adsorption isotherm. A Voigt kernel convolves Gaussian σ=0.35 min with Lorentzian γ=0.12 min. Each component is normalized over the recorded time window, with equal detector response per antibody. Species-specific extinction coefficients, drift, aggregates, non-DAR impurities, peak asymmetry and isomeric separation are absent.

All nine DAR components contribute to the total trace. Even components are emphasized visually without discarding odd components. FWHM is measured numerically; the reported Gaussian-equivalent resolution is Rs=1.18 Δt/(FWHM₁+FWHM₂). Because the kernels contain Lorentzian tails, this is an operational resolution index, not a universal baseline-separation criterion. DAR3/4 and DAR4/5 have Rs 1.71395 and 1.78010. A hypothetical DAR4 collection window from 10.4620 to 13.2898 min gives 96.9147% DAR4 antibody-area purity and 94.6290% DAR4 recovery. Purity counts all tail contamination present in the nine-species model, but says nothing about absent process impurities.

Nonnegative least squares recovers the known noiseless kernel weights to 1.11×10⁻¹⁶. This is an internal arithmetic check. It does not demonstrate identifiability, robustness to noise, suitability for release testing, or accurate analysis of measured chromatograms.

## 8C: finite-bath trafficking and payload stoichiometry

The downstream material is the entire 4-equivalent reaction product, with mean DAR 3.999320; the HIC DAR4 collection fraction is **not** silently substituted as the administered material. Equal binding, trafficking and effective cleavage are assumed across DAR species, allowing mean DAR to set payload stoichiometry. DAR-specific PK or uptake would require separate species states.

Ten states track external ADC X, free surface receptor R, bound ADC B, endosomal ADC E, lysosomal ADC L_y, recycling receptor R_i, degraded antibody Z, cytoplasmic payload P, cumulative exported payload Q and metabolized payload M. Molecular counts are per target-positive cell. The finite bath starts with 10 nM ADC in 1 nL per cell, or 6.02214 million ADC molecules. The task-specified receptor abundance is 10⁶ per cell and internalization rate is 0.05 h⁻¹; neither is asserted to be a measured property of a selected HER2/TROP2 cell line.

Binding is k_on C_ext R − k_off B. Internalization moves bound ADC into E while the receptor enters a recycling pool. ADC sorting transfers E to L_y; receptor recycling is represented independently. Effective lysosomal processing has v=Vmax L_y/(Km+L_y), with hypothetical Vmax=50,000 ADC/cell/h and Km=50,000 ADC/cell. Antibody degradation releases meanDAR×v payload molecules per hour. Cytoplasmic payload then undergoes first-order export and metabolism. The model lumps linker cleavage, self-immolation and antibody processing; it is not a detailed enzyme or bond-by-bond reaction mechanism.

The invariants are X+B+E+L_y+Z=X₀; R+B+R_i=R₀; and meanDAR(X+B+E+L_y)+P+Q+M=meanDAR X₀. Across the two permeability runs, maximum relative antibody/payload errors are below 7.1×10⁻¹⁴ and receptor errors below 4.3×10⁻¹³. Numerical positivity and zero-binding/zero-cleavage limits are checked independently.

The effective protease is deliberately not uniquely cathepsin B. Primary knockout evidence shows other proteases can support processing of tested Val-Cit systems. The generic Val-Cit-like module also does not assert that DS-8201a uses the same linker chemistry. Plasma premature release is not included: the systemic stability-versus-intracellular-release tradeoff remains a measurement and modeling gap.

## 8D: spherical tissue transport and conditional response

A 200-µm-radius sphere has an antigen-positive source core of radius 20 µm. The extracellular fraction is 0.3; the remaining core volume divided by the assumed 2,000 µm³/cell defines an effective 11.73 source cells. Fractions represent volume-averaged tissue, not integer reconstructed cells. The Ag-negative intracellular volume is the non-core cellular fraction. Exact shell/core intersections preserve geometry during mesh refinement.

Extracellular and neighbor-intracellular amounts are stored separately. In shell i, C_e=N_e/V_e and C_i=N_i/V_i. Neighbor uptake is k_perm V_i C_e and outward exchange is k_perm N_i, so intracellular gain is exactly matched by extracellular loss. Shared-face diffusion uses conductance D_eff ε4πr²/Δr. Symmetry imposes zero flux at the origin without evaluating a singular 1/r expression. At 200 µm, an absorbing boundary represents removal into unmodeled surrounding tissue; that outward flux is explicitly counted as a sink. Extracellular clearance and intracellular metabolism are also counted. The supplied amount is the integrated cellular export Q multiplied by effective core cell count.

An implicit Euler solve of the conservative finite-volume operator is positive and stable for this linear transport system. The default grid has 80 shells, Δr=2.5 µm, Δt=0.1 h and output every 0.5 h. It integrates 72 h. The Ag-positive model supplies payload in one direction: tissue payload does not diffuse back into that source cell model. Neighbor exchange is bidirectional. Source death does not alter payload release, receptor abundance or geometry; the exposure-response calculation is an output observation model rather than feedback into transport.

The hypothetical surviving fraction is exp[−∫k_max C_i/(EC50+C_i)dt], with EC50=50 nM and k_max=0.05 h⁻¹. The displayed extracellular threshold contour uses the same assumed concentration, but extracellular threshold crossing is not equivalent to intracellular cytotoxicity. The core response uses its own cytoplasmic payload concentration. The high/low scenarios use permeability rates 0.2/0.002 h⁻¹; they are mechanistic examples inspired by permeability-sensitive bystander findings, not calibrated DXd/MMAF products.

| Default 72-h result | High permeability | Low permeability |
|---|---:|---:|
| Mean Ag-negative surviving fraction, volume weighted | 0.958800 | 0.999848 |
| Ag-positive source surviving fraction, conditional model | 0.049318 | 0.038334 |
| Furthest sampled Ag-negative shell center above extracellular 50 nM | 23.75 µm | None |
| Maximum tissue mass relative error | 1.36×10⁻¹⁴ | 4.20×10⁻¹⁵ |

The high-permeability case produces localized bystander exposure yet only about 4.12% volume-averaged Ag-negative loss across the entire sphere. The computation therefore does not demonstrate broad tumor eradication. The high source response and weak average bystander response illustrate geometry and transport dilution. Retained intracellular payload explains why the low-permeability source has a slightly lower modeled survival while its neighbors remain largely unaffected.

## Verification, sensitivity and what remains

The run includes four conjugation conditions, 24 stochastic replicates, two nominal cellular/transport conditions, two transport refinements and four one-at-a-time diffusion/clearance perturbations. The 23,200-row spatial dataset and mass ledger are in [data](data/). Numerical self-checks include the independent constant-reagent binomial limit, reduction/conjugation blocks, cell stoichiometry, finite-volume column balance, uniform concentration equilibrium, positivity and refinement. Independent tests additionally check closed uniform exponential decay using a matrix-exponential reference.

Halving the time step changes the mean Ag-negative survival by 7.05×10⁻⁶; doubling shell count at the fine time step changes it by 2.12×10⁻⁶. Peak-normalized final extracellular profile error is recorded in [summary.json](summary.json), rather than judging spatial convergence solely from a scalar average. These are discretization checks for this model, not biological validation. [sensitivity.json](sensitivity.json) changes D_eff and extracellular clearance individually by factors 0.5 and 2 while holding the source fixed. It provides conditional sensitivity, not a joint probabilistic uncertainty interval.

To support actual manufacturing or development decisions, the next evidence must include reduction-site distributions and DAR-resolved chemistry, measured HIC traces with response corrections and impurity controls, finite-exposure binding/internalization data, released-species LC-MS kinetics, plasma stability across matrices, permeability and efflux measurements, and spatially resolved coculture or spheroid viability. Parameter estimation and external holdouts would follow. The delivered code and figures make those assumptions and balances reviewable now; they do not claim those missing studies have been performed.
