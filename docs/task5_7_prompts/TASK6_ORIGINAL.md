---

### 第二部分：任务六执行指令（制剂工程：无定形固体分散体 ASD 溶出与过饱和热力学）

> **研究目标：** 现代药物化学合成出的高活性先导化合物中，超过 $70\%$ 属于 BCS II/IV 类（高脂溶性、微溶于水）。要让它们真正变成口服药，必须依赖**物理药剂学与制剂工程**。任务六构建聚合物无定形固体分散体（Amorphous Solid Dispersion, ASD）的热力学相分离理论（Flory-Huggins 模型）以及体外“弹簧-降落伞（Spring-and-Parachute）”过饱和溶出-结晶动力学仿真。

```markdown
[TASK 6: PHARMACEUTICS & FORMULATION ENGINEERING: AMORPHOUS SOLID DISPERSIONS (ASD) & SUPERSATURATION DYNAMICS]

NOTE ON SUBSCRIPTION RATE-LIMIT CONSERVATION:
Execute under strict single-turn delivery guidelines. Do NOT engage in multi-turn conversational trials. Deliver the complete, production-grade, bug-free Python script `run_task6_asd_formulation_supersaturation_kinetics.py` in a single, fully-contained, finalized artifact.

---

### 1. FORMULATION ENGINEERING CONTEXT: THE BCS CLASS II/IV CRISIS
Over $70\%$ of pipeline drug candidates suffer from poor aqueous solubility (BCS Class II/IV). Without advanced formulation engineering, potent molecules fail in vivo due to dissolution-limited absorption.
Task 6 models the fundamental physical chemistry of **Amorphous Solid Dispersions (ASDs)**: utilizing Flory-Huggins polymer-drug thermodynamics to predict miscibility and Gordon-Taylor glass transition ($T_g$) depression, coupled with a non-sink dissolution-precipitation kinetic engine modeling the classical **"Spring-and-Parachute"** supersaturation phenomenon.

---

### 2. FOUR-TIER FORMULATION COMPUTATIONAL PROTOCOL

#### MODULE 6A: Polymer-Drug Phase Equilibrium & Flory-Huggins Thermodynamics
- Screen a panel of 4 pharmaceutical polymeric carriers: **HPMC-AS (hypromellose acetate succinate), PVP-VA (copovidone), Soluplus, and Eudragit L100**.
- Calculate the polymer-drug Flory-Huggins interaction parameter ($\chi_{\text{drug-poly}}$) based on partial Hansen Solubility Parameters ($\delta_D, \delta_P, \delta_H$):
  $$\chi = \frac{V_{\text{drug}}}{R T} \left[ (\delta_{D,1} - \delta_{D,2})^2 + 0.25 (\delta_{P,1} - \delta_{P,2})^2 + 0.25 (\delta_{H,1} - \delta_{H,2})^2 \right]$$
- Evaluate the Gibbs Free Energy of Mixing ($\Delta G_{\text{mix}}$) across drug weight fractions ($w_{\text{drug}} \in [0, 1]$):
  $$\frac{\Delta G_{\text{mix}}}{R T} = \frac{\phi_{\text{drug}}}{N_{\text{drug}}} \ln \phi_{\text{drug}} + \frac{\phi_{\text{poly}}}{N_{\text{poly}}} \ln \phi_{\text{poly}} + \chi \phi_{\text{drug}} \phi_{\text{poly}}$$
- Generate the binodal and spinodal phase diagrams: Delineate the thermodynamically stable, metastable, and phase-separated amorphous-amorphous phase separation (AAPS) domains.

#### MODULE 6B: Glass Transition ($T_g$) & Physical Storage Stability (Gordon-Taylor Engine)
- Model the glass transition temperature of the binary solid dispersion ($T_{g,\text{mix}}$) via the Gordon-Taylor equation:
  $$T_{g,\text{mix}} = \frac{w_1 T_{g,1} + K w_2 T_{g,2}}{w_1 + K w_2}, \quad K \approx \frac{\rho_1 T_{g,1}}{\rho_2 T_{g,2}}$$
- Quantify physical stability against recrystallization at storage conditions ($25^\circ\text{C} / 60\% \text{ RH}$ and $40^\circ\text{C} / 75\% \text{ RH}$):
  - Predict molecular mobility via the Adam-Gibbs / Vogel-Fulcher-Tammann (VFT) relaxation model.
  - Calculate recrystallization induction time ($t_{\text{ind}}$) as a function of drug loading ($10\%, 20\%, 30\%, 50\%$).

#### MODULE 6C: The "Spring-and-Parachute" Non-Sink Dissolution Kinetics
- Formulate a coupled dissolution-precipitation ODE system modeling the gastrointestinal luminal microenvironment ($V_{\text{dissolution}} = 900\text{ mL}$, $\text{pH} = 6.8$ with bile salts):
  $$\frac{dC}{dt} = \left( \frac{dC}{dt} \right)_{\text{dissolution}} - \left( \frac{dC}{dt} \right)_{\text{precipitation}}$$
- **Dissolution Engine (The Spring):** Modified Noyes-Whitney equation driven by the amorphous solubility ($C_{\text{amorphous}} = C_{\text{crystalline}} \times \exp(\Delta G_{\text{crys}} / RT)$):
  $$\left( \frac{dC}{dt} \right)_{\text{diss}} = \frac{D \cdot S(t)}{h} \left( C_{\text{amorphous}} - C(t) \right)$$
- **Precipitation Engine (The Parachute):** Classical Nucleation Theory (CNT) coupled with crystal growth:
  $$\left( \frac{dC}{dt} \right)_{\text{precip}} = k_{\text{precip}} \cdot (C(t) - C_{\text{crystalline}})^n \cdot \exp\left( - \frac{\Delta G_{\text{nucleation}}}{(k_B T)^3 (\ln S)^2} \right)$$
  where the polymeric precipitation inhibitor (e.g., HPMC-AS) alters the interfacial energy ($\gamma$), effectively acting as the "parachute" to sustain supersaturation for $> 4\text{ hours}$.

#### MODULE 6D: In-Vitro to In-Vivo Luminal Bioavailability Simulation
- Simulate three comparative formulations under non-sink biorelevant conditions (FaSSIF):
  1. Unformulated micronized crystalline drug.
  2. Pure amorphous drug without polymer.
  3. Optimized Polymer-ASD formulation ($20\%$ drug in HPMC-AS).
- Calculate the **Supersaturation Area Under the Curve ($\text{AUC}_{\text{dissolution}}$)** and the theoretical intestinal absorption enhancement ratio ($\text{ER}_{\text{abs}}$).

---

### 3. PUBLICATION-QUALITY DELIVERABLES (300 DPI in `./figures_task6/`)
Write code to compute and save:
1. **`fig1_flory_huggins_miscibility_phase_diagram.png`**: Gibbs free energy of mixing curves ($\Delta G_{\text{mix}}$ vs. Drug Fraction) and phase diagrams identifying the spinodal decomposition boundary for 4 candidate polymers.
2. **`fig2_gordon_taylor_tg_depression.png`**: Predicted $T_g$ curves as a function of drug loading alongside stability risk zones ($T_g - T_{\text{storage}} < 30^\circ\text{C}$).
3. **`fig3_spring_and_parachute_dissolution.png`**: Kinetic concentration-time profiles ($C(t)$ vs. Time) demonstrating the rapid dissolution peak ("Spring") and sustained polymeric inhibition of precipitation ("Parachute") vs. crystalline baseline.
4. **`fig4_polymeric_precipitation_inhibition_efficiency.png`**: Bar chart quantifying crystallization induction times ($t_{\text{ind}}$) and dissolution $\text{AUC}$ ratios across varying polymer-to-drug ratios.

---

### 4. BILINGUAL SCHOLARLY TREATISE & REPOSITORY EVOLUTION
- Author **`ASD_FORMULATION_KINETICS_REPORT_ZH.md`** and **`ASD_FORMULATION_KINETICS_REPORT_EN.md`**:
  - The thermodynamic dilemma of amorphous drugs: High free-energy solubility vs. spontaneous recrystallization driving forces.
  - Polymer selection rationale: Hydrogen-bonding donor/acceptor complementarity and Hansen solubility parameter matching.
  - Translating in-vitro "Spring-and-Parachute" supersaturation into human clinical oral exposure gains.
- Update root `README.md` to index Task 6.
- Execute automated commit and push:
  ```bash
  git add .
  git commit -m "feat(task6): amorphous solid dispersion (ASD) formulation thermodynamics & supersaturation kinetics"
  git push origin main