[TASK 5: CYTOCHROME P450 DRUG-DRUG INTERACTIONS (DDI) & TIME-DEPENDENT INHIBITION (TDI/MBI) DYNAMICS]

NOTE ON SUBSCRIPTION RATE-LIMIT CONSERVATION:
Execute under strict single-turn delivery guidelines. Do NOT engage in multi-turn conversational trials. Deliver the complete, production-grade, bug-free Python script `run_task5_cyp_ddi_mechanism_based_inhibition.py` in a single, fully-contained, finalized artifact.

---

### 1. CLINICAL TRANSLATIONAL CONTEXT: AVOIDING BLACK-BOX WARNINGS
Drug-Drug Interactions (DDIs) caused by the inhibition of Cytochrome P450 (CYP) enzymes represent a leading cause of late-stage clinical attrition and post-market withdrawals. 
Task 5 constructs a multi-enzyme mechanistic DDI simulation framework quantifying both reversible competitive inhibition ($K_i$) and irreversible Time-Dependent Inhibition / Mechanism-Based Inactivation (TDI/MBI: $k_{\text{inact}}, K_I$).

---

### 2. FOUR-TIER DDI PHARMACOKINETIC PROTOCOL

#### MODULE 5A: CYP Isoform Panel & In-Vitro Kinetic Parameterization
- Model a clinical panel of the three major hepatic metabolizing isoforms: **CYP3A4, CYP2D6, and CYP2C9**.
- Ingest compound parameters across 6 representative clinical test molecules (Strong inhibitor: Ketoconazole; MBI inactivator: Clarithromycin/Ritonavir; Moderate/weak inhibitors; Non-inhibitor):
  - Reversible competitive inhibition constant: $K_i$ ($\mu\text{M}$).
  - Time-dependent inactivation parameters: Maximum inactivation rate $k_{\text{inact}}$ ($\text{h}^{-1}$) and apparent inactivation constant $K_I$ ($\mu\text{M}$).
  - Natural degradation/turnover rate of hepatic CYP enzymes: $k_{\text{deg}}$ (e.g., $k_{\text{deg, 3A4}} = 0.019\text{ h}^{-1}$, $t_{1/2} = 36\text{ h}$).

#### MODULE 5B: Mechanistic Static DDI Prediction (FDA Regulatory Criteria)
- Implement official FDA/EMA guidance static risk models evaluating the area under the curve ratio ($\text{AUCR}$):
  - **Reversible Inhibition Ratio ($R_1$):**
    $$R_1 = 1 + \frac{[I]_{\max, u}}{K_i}$$
  - **Mechanism-Based Inactivation Ratio ($R_2$):**
    $$R_2 = \frac{k_{\text{obs}} + k_{\text{deg}}}{k_{\text{deg}}} = \frac{\frac{k_{\text{inact}} \times [I]_{\text{gut/in vivo}}}{K_I + [I]_{\text{gut/in vivo}}} + k_{\text{deg}}}{k_{\text{deg}}}$$
  - **Total Systemic & Intestinal Mechanistic Net Ratio:**
    $$\text{AUCR}_{\text{net}} = \left( \frac{f_m}{R_{\text{hepatic}} \times (1 - E_H) + E_H} + (1 - f_m) \right) \times \frac{1}{R_{\text{gut}} \times (1 - F_g) + F_g}$$
    where $f_m$ is the fraction of victim probe drug cleared by the specific CYP isoform (e.g., $f_m = 0.9$ for Midazolam via CYP3A4).

#### MODULE 5C: Dynamic PBPK-Coupled DDI Co-Administration Kinetics
- Couple the Task 4 PBPK multi-compartment solver to simulate the simultaneous co-administration of the perpetrator lead compound with clinical victim index drugs:
  - **Probe Drug 1:** Midazolam ($2\text{ mg}$ oral, sensitive CYP3A4 substrate).
  - **Probe Drug 2:** Metoprolol ($50\text{ mg}$ oral, CYP2D6 substrate).
- Solve the dynamic depletion and de novo resynthesis of active functional enzyme abundance:
  $$\frac{d[\text{CYP}]_{\text{active}}}{dt} = k_{\text{syn}} - k_{\text{deg}} [\text{CYP}]_{\text{active}} - \frac{k_{\text{inact}} [I(t)]_{\text{liver}, u}}{K_I + [I(t)]_{\text{liver}, u}} [\text{CYP}]_{\text{active}}$$
- Simulate repeated BID dosing of perpetrator for 14 days, observing the progressive suppression of active enzyme capacity and the resulting plasma concentration spike of the victim drug.

#### MODULE 5D: Clinical DDI Classification & Washout Window
- Automatically categorize the DDI severity index according to regulatory thresholds:
  - **Strong DDI:** $\text{AUCR} \ge 5.0$ (High regulatory barrier; contraindicated).
  - **Moderate DDI:** $2.0 \le \text{AUCR} < 5.0$.
  - **Weak DDI:** $1.25 \le \text{AUCR} < 2.0$.
- Quantify the **Enzyme Recovery Washout Half-Life**: Calculate the number of days required post-cessation of perpetrator dosing for CYP active pools to regenerate back to $90\%$ baseline.

---

### 3. PUBLICATION-QUALITY DELIVERABLES (300 DPI in `./figures_task5/`)
Write code to compute and save:
1. **`fig1_cyp_inactivation_curves.png`**: Time- and concentration-dependent pre-incubation enzyme activity loss curves identifying $k_{\text{inact}}$ and $K_I$ for MBI compounds.
2. **`fig2_dynamic_ddi_midazolam_pk.png`**: Simulated clinical plasma concentration-time curves ($C_p - t$) of oral Midazolam in the absence (Control) vs. presence of perpetrator lead drug at steady state.
3. **`fig3_fda_ddi_risk_heatmap.png`**: 2D regulatory risk matrix mapping $[I]_{\max, u} / K_i$ vs. $k_{\text{inact}} / K_I$ with color-coded safety boundaries (Green: Safe, Yellow: Clinical Study Warranted, Red: Severe Liability).
4. **`fig4_cyp3a4_resynthesis_timeline.png`**: Active enzyme recovery trajectory post-perpetrator withdrawal tracking the regeneration of clearance capacity over 10 days.

---

### 4. BILINGUAL SCHOLARLY TREATISE & REPOSITORY EVOLUTION
- Author **`CYP_DDI_KINETICS_REPORT_ZH.md`** and **`CYP_DDI_KINETICS_REPORT_EN.md`**:
  - The biochemical mechanism of Cytochrome P450 heme coordinate destruction vs. apoprotein covalent adduct formation in MBI.
  - Mathematical comparison between static $R$-value estimation and dynamic PBPK co-administration models.
  - Mitigating DDI liabilities during lead optimization: Structure-activity relationships (SAR) to reduce metabolic activation (e.g., furan, aniline, or methylenedioxyphenyl ring replacements).
- Update root `README.md` to index Task 5.
- Execute automated commit and push:
  ```bash
  git add .
  git commit -m "feat(task5): CYP450 DDI kinetics, mechanism-based inactivation (MBI) & clinical AUCR prediction"
  git push origin main