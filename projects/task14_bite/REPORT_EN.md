# Task 14 — BiTE molecular crosslinking and target-cell loss

For the nominal asymmetric affinities (CD3 KD=100 nM, TAA KD=1 nM), the sampled molecular-trimer maximum occurs at **26.1016 nM total BiTE**, with **572.752 complexes per target cell**. At 100 µM BiTE, crosslinking falls to 0.144815% of that maximum. At a maintained 10 nM exposure, the model leaves **3.82331%** of initial target cells after 48 h. These are conditional model outputs, not clinical response or safety thresholds. Raw data are in [binding grid](outputs/binding_grid.csv), [lysis grid](outputs/lysis_grid.csv) and [summary](outputs/summary.json).

## Three-body equilibrium and units

Let $e,t,b$ be free CD3, TAA and BiTE concentrations in nM. Binary complexes are $EB=eb/K_E$ and $TB=tb/K_T$; trimer is $X=\alpha etb/(K_EK_T)$, with nominal $\alpha=1$. Three exact finite balances are enforced:

$$E_0=e+EB+X,\quad T_0=t+TB+X,\quad B_0=b+EB+TB+X.$$

Eliminate $t=T_0/[1+b/K_T+\alpha eb/(K_EK_T)]$. For each trial free $b$, a bounded root determines free $e$ from the first balance; a second bounded root determines $b$ from total BiTE. No approximation identifies free BiTE with its total dose. Reversing the two terminal receptor labels gives the same trimer, independently tested. With one receptor absent the remaining binary equilibrium agrees with the stable quadratic solution.

The effective contact volume is an assumed $v=10^{-12}$ L per initial target cell. Accessible receptor copies are 5,000 CD3 per effector and 10,000 TAA per target; these are declared coarse-grained copy numbers, not measured receptor abundance. For surviving target fraction $n=N/N_0$ and fixed effector/initial-target ratio $r$, concentrations are

$$E_0=\frac{r\,n_{CD3}}{N_Av}10^9\ \mathrm{nM},\qquad
T_0(n)=\frac{n\,n_{TAA}}{N_Av}10^9\ \mathrm{nM},$$
$$S(n)=X(n)10^{-9}N_Av/n\quad\text{molecular trimers per surviving target cell}.$$

$N_A$ is Avogadro's exact constant. The calculation preserves finite receptor totals even while targets disappear. It reports molecular crosslinks, not a count of whole cell–cell immunological synapses. A conversion to actual cell pairs would require a spatial contact and synapse-population model.

## Lysis dynamics and controls

$$\dot n=-k_{max}\frac{S(n)^h}{K_S^h+S(n)^h}n,\qquad
\dot d=+k_{max}\frac{S(n)^h}{K_S^h+S(n)^h}n.$$

The illustrative parameters are $k_{max}=0.12$ h⁻¹, $K_S=500$ complexes/cell, $h=2$, $n(0)=1$, $d(0)=0$. Thus $n+d=1$ exactly. Effector count and externally maintained total BiTE remain fixed. Re-equilibration releases ligand from lost receptor sites into the maintained pool; antibody PK or irreversible antibody consumption is not modeled. The source compartment maintaining exposure lies outside the finite-receptor model.

The run computes **981 binding equilibria** (109 doses × three CD3 affinities × three TAA affinities) and **270 lysis ODEs** (30 doses × three E:T ratios × three TAA expression levels). Binding-grid mass error is below $8.62\times10^{-11}$ nM, cell loss-ledger error below $4.45\times10^{-16}$, and solver refinement changes cell fractions by less than $5.68\times10^{-7}$. Zero drug, zero effector or zero crosslink cooperativity gives exactly zero lysis. Full assumptions and checks are in [config](outputs/config.json) and [verification](outputs/verification.json).

The high-dose hook follows competition between binary saturation and ternary bridging in this particular reversible equilibrium. Its position moves with affinities and receptor abundance; a universal clinical self-inhibition window above 1 µM is not demonstrated. Conditional TCE architectures, nonequilibrium binding, drug PK, cell motility, serial killing, activation thresholds, exhaustion, antigen shedding and cytokine toxicity can alter the observed response. Industry-facing validation needs quantitative receptor copy counts, binding kinetics, cell-pair imaging and matched dose/time killing assays before this model can guide a specific molecule. Primary mechanistic studies are in [sources](sources.md).
