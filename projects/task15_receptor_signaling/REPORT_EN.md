# Task 15 — Finite-receptor proofreading and cytokine signaling

This executed model links engineered receptor dwell-time discrimination to a bounded immune-cell activation network. It predicts normalized cytokine trajectories and IL6R signal inhibition. Signal reduction does not guarantee lower ligand concentration or restoration of a defined homeostatic set point. See the [source register](sources.md), [numerical summary](outputs/summary.json), [verification](outputs/verification.json) and [figure](outputs/figures/fig15_receptor_proofreading_cytokine_balance.png).

## Model and derivation

For antigen class a, free receptor R enters C(a,0) with pseudo-first-order association k(a)R. Five phosphorylation transitions occur at rate kp. Every bound state dissociates at koff(a), returning the receptor to R. With fixed antigen reservoirs:

$$\dot C_{a,0}=k_aR-(k_p+k_{off,a})C_{a,0},$$
$$\dot C_{a,j}=k_pC_{a,j-1}-(k_p+k_{off,a})C_{a,j},\quad j=1,\ldots,4,$$
$$\dot C_{a,5}=k_pC_{a,4}-k_{off,a}C_{a,5}.$$

The free-receptor derivative is the sum of all dissociation fluxes minus all binding fluxes. Adding equations cancels every internal transition, giving R+sum(C)=1. Ligand itself is not conserved: its concentration is a specified external reservoir, with the target reservoir declining after 24 h.

At steady state, R=1/[1+sum(k(a)/koff(a))]. Summing each antigen's bound states gives k(a)R/koff(a), while the terminal signaling fraction among bound receptors is [kp/(kp+koff(a))]^5. This analytical relation supplies an independent check of the simulated reaction topology. The default target/self terminal-state ratio is **59.499**, conditional on assumed rates; it is not experimentally measured antigen specificity.

Let T and M be active fractions, and I6, N, F and I1 be normalized IL6, TNF, IFN and IL1 concentrations. All four cytokines start at zero, as do T and M. Define trigger q=(Ctarget,5+Cself,5)/(0.15+Ctarget,5+Cself,5), antagonist occupancy B=X/(1+X), and IL6R signal s=I6/(1+I6)(1-B). The cytokine network is:

$$\dot T=0.22q(1-T)-0.06T,$$
$$\dot M=[0.35F/(1+F)+0.20s](1-M)-0.10M,$$
$$\dot I_6=0.18T+1.10M-[0.12+0.08(1-B)]I_6,$$
$$\dot N=0.55T+0.60M-0.30N,\quad\dot F=0.90T+0.10M-0.25F,$$
$$\dot I_1=0.08T+0.80M-0.20I_1.$$

Time is hours. T and M remain between zero and one because their boundary derivatives point inward. The concentration unit and all rate constants are hypothetical. Antagonist X is zero before 12 h and then X0 exp[-ln(2)(t-12)/48]. The numerical solver integrates the intervention boundary in separate segments.

## Executed results and checks

Four 96-hour scenarios contain 3,844 time rows; a 65-point self-antigen off-rate scan evaluates discrimination. IL6R signal AUC is **75.201, 51.158, 17.106 and 6.781 h** for initial X0 of 0, 1, 10 and 100. The corresponding IL6 concentration peaks are **5.108, 5.944, 7.186 and 7.587 normalized units**. A lower receptor signal can coexist with higher ligand concentration because receptor-linked removal is blocked. The result prevents the incorrect assumption that anti-IL6R action must eliminate IL6 production.

Checks cover receptor conservation, finite nonnegative concentrations, analytical steady-state probabilities, tighter-tolerance integration, equal-antigen symmetry and distinction between ligand concentration and receptor signal. Solver verification is evidence of numerical consistency, not biological accuracy.

## Translational use and evidence gap

The model is suitable for designing receptor dwell-time and cytokine time-course measurements. Calibrating it would require construct-specific binding/phosphorylation data, antigen density, macrophage activation assays, cytokine units and an independently evaluated antagonist exposure/occupancy model. Vascular leak, endothelial response, organ injury, infections and neurotoxicity are absent. Consequently `homeostasis_restoration_established` and `clinical_recommendations` remain null. Antagonist C/KD is not a tocilizumab dose.
