# Task 12 — Competitive chromatin-target degradation

At the assumed 100 nM degrader dose, 48 h total target remains **51.2557 nM** (51.2557% of initial), comprising a 22.7889 nM free-context pool and a 28.4668 nM chromatin-context pool. These pools include their drug-bound and ternary states; they are not the unliganded species alone. In the matched untreated model, the initial 30:70 split relaxes to 37.7776:62.2224 nM at 48 h. This baseline redistribution must not be mislabeled as drug-induced degradation. Exact results are in [summary](outputs/summary.json).

## Binding thermodynamics and explicit kinetics

Let unliganded target pools be $F,C$, degrader $P$, and ligase $E$. Binary species are $PF,PC,EP$; ternary species are $EFP,ECP$. Equilibrium binding alone obeys

$$PF=p f/K_T,\quad PC=p c/K_T,\quad EP=e p/K_E,$$
$$EFP=\alpha_F e p f/(K_EK_T),\quad ECP=\alpha_C e p c/(K_EK_T),$$
$$\alpha_C=\alpha_F\exp(-\Delta G_{occlusion}/RT).$$

Here $K_E=50$ nM, $K_T=100$ nM, $\alpha_F=10$, and $RT=0.61633$ kcal mol⁻¹. The nominal hypothetical penalty is 2 kcal mol⁻¹, giving $\alpha_C=0.38968$. For fixed context totals, eliminate $f=F_{tot}/[1+p/K_T+\alpha_Fe p/(K_EK_T)]$ and analogously $c$. Nested bounded roots solve the remaining finite $E_{tot}=20$ nM and $P_{tot}$ balances. This separate binding-only calculation does not assume the degrading system remains at equilibrium.

The kinetic simulation explicitly integrates both association paths: $EP+F\rightleftharpoons EFP$ and $PF+E\rightleftharpoons EFP$, with corresponding chromatin reactions. The shared association constant is 0.01 nM⁻¹ h⁻¹; binary off-rates are $k_{on}K_D$, and ternary off-rates divide the relevant binary off-rate by $\alpha_F$ or $\alpha_C$. Both paths therefore have the same equilibrium product, satisfying the binding thermodynamic cycle. For each reaction, the derivative is assembled as stoichiometry times net forward-minus-reverse flux, in nM/h.

Only unliganded targets exchange contexts, with $J_{chrom}=k_{cap}F-k_{rel}C$, where $k_{cap}=0.14$ h⁻¹ and $k_{rel}=0.06$ h⁻¹. Every target-containing species undergoes constitutive turnover at $k_0=0.025$ h⁻¹, releasing surviving degrader/ligase. New target enters $F$ at $s=k_0T_0=2.5$ nM/h. Ternary degradation occurs at $k_{cat}=0.5$ h⁻¹ and also releases $E$ and $P$. Adding all target-containing derivatives cancels binding/exchange:

$$\dot T=s-k_0T-k_{cat}(EFP+ECP).$$

Integrated ledgers satisfy $T+D-S=T_0=100$ nM; the sums of all degrader-containing and all ligase-containing species remain exactly their initial totals. Without induced degradation and starting at $T_0$, total target remains 100 nM. Without drug, $C(t)=C_\infty+[C(0)-C_\infty]e^{-\lambda t}$, with $\lambda=k_{cap}+k_{rel}+k_0$ and $C_\infty=k_{cap}T_0/\lambda$. This analytic baseline is independently tested.

## Computation, use and limits

The run evaluates **152 dynamic scenarios** (38 doses including zero × four penalties 0, 1, 2, 4 kcal mol⁻¹) and 152 corresponding binding equilibria, plus no-drug, no-induced-degradation and tighter-tolerance controls. Largest grid mass-balance error is $1.29\times10^{-9}$ nM; refinement changes species by less than $7.82\times10^{-7}$ nM. Data and units are in [trajectories](outputs/trajectories.csv), [dose grid](outputs/dose_penalty_grid.csv), [binding grid](outputs/binding_equilibria.csv) and [config](outputs/config.json).

This model can compare experimental designs for time-resolved soluble/chromatin fractionation, ligase limitation, washout, and chromatin accessibility perturbations. It does not ingest an actual nucleosome structure, infer steric occlusion from coordinates, model nucleosome-site saturation, distinguish individual complex subunits, or calculate a nucleosome uncoupling energy. Its chromatin penalty is a declared assumption rather than a fitted or quantum-computed energy. Literature motivation is recorded separately in [sources](sources.md); no clinical selectivity or efficacy is inferred from the simulated depletion.
