# Task 13 — TfR shuttle transport with finite antibody and receptor pools

The nominal scenario gives a **0–72 h extravascular brain/plasma concentration AUC ratio of 0.00227242** and brain AUC of **6.69954 nM·h**. The affinity sweep peaks at **13.3352 nM apical KD**, outside the requested 50–500 nM range. This is retained as an actual model outcome. It is not a measured affinity optimum or human dose prediction. See [summary](outputs/summary.json), [grid](outputs/affinity_sensitivity.csv) and [trajectories](outputs/trajectories.csv).

## Amount balances and pH-dependent binding

State variables are amounts (nmol): plasma antibody $A_p$, unoccupied surface receptor $R_s$, surface complex $C_s$, endosomal receptor $R_e$, endosomal complex $C_e$, free endosomal antibody $A_e$, and free brain antibody $A_b$. Three cumulative losses record systemic clearance, lysosomal degradation and brain clearance. Concentrations are amount/volume: nM = nmol/L. Volumes are assumed $V_p=3$ L, $V_e=0.001$ L and $V_b=0.3$ L; initial plasma concentration is 100 nM, giving 300 nmol antibody. Total receptor is 0.03 nmol and starts at its unoccupied recycling/internalization steady state.

Surface and endosomal net association fluxes are

$$J_s=k_{on}[(A_p/V_p)R_s-K_{D,s}C_s],\qquad
J_e=k_{on}[(A_e/V_e)R_e-K_{D,e}(t)C_e].$$

$k_{on}=0.02$ nM⁻¹h⁻¹; both fluxes have units nmol/h. The common endosomal pH closure is $pH_e(t)=5.8+1.6e^{-t/(1h)}$, and $K_{D,e}=K_{D,s}10^{q(7.4-pH_e)}$. The assumed slope $q=0.7$ log10(KD)/pH makes acidification weaken binding. It is varied independently; it is not a measured histidine-dependent off-rate profile.

Let $I_C=k_iC_s$, $I_R=k_{ir}R_s$, $Q=k_rR_e$, $D_C=k_{dc}C_e$, $D_F=k_{df}A_e$, $T_b=k_bA_e$, $T_p=k_pA_e$, $L_p=k_{cl,p}A_p$, and $L_b=k_{cl,b}A_b$. The seven inventory equations are

$$\dot A_p=-J_s+T_p-L_p,\quad \dot R_s=-J_s-I_R+Q,\quad \dot C_s=J_s-I_C,$$
$$\dot R_e=I_R-Q-J_e+D_C,\quad \dot C_e=I_C+J_e-D_C,$$
$$\dot A_e=-J_e-D_F-T_b-T_p,\quad \dot A_b=T_b-L_b.$$

Each loss receives its named clearance flux; the lysosome ledger receives $D_C+D_F$. Adding inventories yields exact conservation of initial antibody plus losses. Adding receptor species yields $R_s+C_s+R_e+C_e=R_{tot}$. Cargo degradation releases receptor back into the endosomal pool. This deliberately omits receptor downregulation, which must be added before comparing to systems in which high-affinity antibodies degrade TfR itself. No antibody is injected from an infinite reservoir.

## Exposure definition and checks

$$K_{p,0-72h}=\frac{\int_0^{72h}(A_b/V_b)dt}{\int_0^{72h}(A_p/V_p)dt}.$$

There is no brain target-binding state, so all brain antibody is free within this model. The result is not a general target-binding-antibody $K_{p,uu}$ estimate. A deliberately contaminated measurement adds vascular concentration $C_pV_{vascular}/V_b$, with $V_{vascular}=0.015$ L; the apparent ratio becomes **0.0522724**, exactly 0.05 above the extravascular ratio. Surface and endothelial complexes are excluded from the parenchymal measurement. Their contribution would require an additional assay-specific correction.

**328 ODE scenarios** vary KD over 0.1–10,000 nM, pH slopes 0/0.7/1.4, initial doses 1/10/100/1000 nM, and receptor amounts 0.01/0.03/0.1 nmol using one-factor sensitivity slices, not a full Cartesian design. Zero receptor or zero complex internalization produces no brain delivery. The independent zero-receptor test matches $A_p(t)=A_p(0)e^{-0.03t}$. Grid antibody-balance error is below $3.83\times10^{-10}$ nmol; tighter-solver relative state difference is $1.34\times10^{-8}$. See [verification](outputs/verification.json) and all [parameters](outputs/config.json).

The industry-relevant output is the sensitivity of release/retention tradeoffs to assay definitions and pH-dependent kinetics. Actual selection needs species-specific TfR abundance, valency/avidity, endogenous transferrin competition, receptor turnover, brain target engagement, and capillary-depleted exposure measurements. Primary observations supporting the mechanism, without parameter transfer, are in [sources](sources.md).
