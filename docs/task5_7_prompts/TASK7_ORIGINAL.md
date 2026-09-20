---

### 第三部分：任务七执行指令（定量系统药理学 QSP：肿瘤免疫联合用药动态协同）

> **研究目标：** 从单纯的“体内浓度（PK）”跨入“真实疾病系统药效（PD）”。任务七构建工业级**定量系统药理学（Quantitative Systems Pharmacology, QSP）**模型，联立 Simeoni 肿瘤生长抑制动力学（TGI）与免疫检查点（PD-1/PD-L1）细胞动力学微分方程，模拟小分子激酶抑制剂与大分子单抗联合用药，利用 **Bliss 独立模型与 Loewe 加和模型**在虚构病人队列中求解真正的抗癌协同增效窗口。

```markdown
[TASK 7: QUANTITATIVE SYSTEMS PHARMACOLOGY (QSP) & IMMUNO-ONCOLOGY COMBINATION SYNERGY]

NOTE ON SUBSCRIPTION RATE-LIMIT CONSERVATION:
Execute under strict single-turn delivery guidelines. Do NOT engage in multi-turn conversational trials. Deliver the complete, production-grade, bug-free Python script `run_task7_qsp_tumor_immune_pkpd_synergy.py` in a single, fully-contained, finalized artifact.

---

### 1. TRANSLATIONAL SYSTEMS PHARMACOLOGY CORE
Modern oncology rarely relies on monotherapy. Clinical cures require dynamic synergies between targeted small molecules (e.g., KRAS/EGFR inhibitors) and macromolecular immunotherapy (e.g., Anti-PD-1 antibodies).
Task 7 constructs a multi-scale **Quantitative Systems Pharmacology (QSP)** model translating drug exposure directly into dynamic tumor regression, integrating the classical Simeoni Tumor Growth Inhibition (TGI) differential framework with cytotoxic T-lymphocyte (CTL) immune infiltration dynamics.

---

### 2. FOUR-TIER QSP PHARMACODYNAMIC PROTOCOL

#### MODULE 7A: Simeoni Tumor Growth Inhibition (TGI) Differential Framework
- Model unperturbed tumor growth incorporating exponential growth transitioning into linear growth governed by carrying capacity:
  $$\frac{dw_0}{dt} = \frac{\lambda_0 w_0}{\left[ 1 + (\frac{\lambda_0}{\lambda_1} w_0)^\psi \right]^{1/\psi}}$$
  where $w_0$ is proliferating tumor mass ($\text{mg}$), $\lambda_0$ is initial exponential rate, and $\lambda_1$ is linear growth rate.
- Implement the drug-induced cell transit damage model: Small molecule targeted drug exposure induces cell death passing through 3 successive non-proliferating apoptotic transit compartments ($w_1, w_2, w_3$):
  $$\frac{dw_0}{dt} = \text{Growth}(w_0) - k_2 C_{\text{drug}}(t) w_0$$
  $$\frac{dw_1}{dt} = k_2 C_{\text{drug}}(t) w_0 - \frac{1}{\tau} w_1$$
  $$\frac{dw_2}{dt} = \frac{1}{\tau} (w_1 - w_2), \quad \frac{dw_3}{dt} = \frac{1}{\tau} (w_2 - w_3)$$
  $$\text{Total Tumor Weight: } W_{\text{total}}(t) = w_0 + w_1 + w_2 + w_3$$
  where $\tau$ is the apoptotic transit mean lifetime.

#### MODULE 7B: Immune Checkpoint & Cytotoxic T-Lymphocyte (CTL) Modulation
- Couple small-molecule tumor cell killing with host anti-tumor immune response:
  - Model tumor antigen release proportional to cell death flux ($\text{Flux}_{\text{antigen}} = w_3 / \tau$).
  - Infiltration of active cytotoxic T-lymphocytes ($E_{\text{CTL}}$):
    $$\frac{dE_{\text{CTL}}}{dt} = \alpha_{\text{stim}} \frac{\text{Flux}_{\text{antigen}}}{K_{\text{ag}} + \text{Flux}_{\text{antigen}}} - \mu_E E_{\text{CTL}} - k_{\text{exh}} \cdot \text{PDL1}_{\text{tumor}} \cdot E_{\text{CTL}}$$
- **Anti-PD-1 Antibody Mechanism:**
  - Model monoclonal antibody pharmacokinetics ($C_{\text{mAb}}(t)$) and target receptor occupancy:
    $$\text{RO}(t) = \frac{C_{\text{mAb}}(t)}{K_D + C_{\text{mAb}}(t)}$$
  - Anti-PD-1 binding competitively suppresses T-cell exhaustion ($k_{\text{exh}} \to k_{\text{exh}} \times (1 - \text{RO}(t))$), restoring CTL-mediated tumor lysis:
    $$\left( \frac{dw_0}{dt} \right)_{\text{immune kill}} = - k_{\text{kill}} \cdot E_{\text{CTL}} \cdot w_0$$

#### MODULE 7C: Synergy Quantification (Bliss Independence vs. Loewe Additivity)
- Simulate 4 comparative clinical arms across a 60-day therapeutic regimen:
  1. **Arm 1 (Vehicle Control):** Unperturbed aggressive tumor growth.
  2. **Arm 2 (Targeted Monotherapy):** Small molecule oral daily dosing ($50\text{ mg}$ QD).
  3. **Arm 3 (Immuno Monotherapy):** Anti-PD-1 IV infusion every 2 weeks ($10\text{ mg/kg}$ Q2W).
  4. **Arm 4 (Combination Therapy):** Simultaneous small molecule + Anti-PD-1.
- Calculate the **Excess Over Bliss (EOB)** synergy score across dynamic time points:
  $$\Delta \text{Bliss}(t) = f_{\text{combo}}(t) - \left[ f_{\text{targeted}}(t) + f_{\text{immuno}}(t) - f_{\text{targeted}}(t) \cdot f_{\text{immuno}}(t) \right]$$
  where $\Delta \text{Bliss} > 0$ strictly demonstrates true synergistic tumor suppression exceeding additive expectations.

#### MODULE 7D: Virtual Patient Cohort Clinical Trial Simulation ($n = 50$)
- Ingest biological inter-individual variability (IIV) using log-normal distributions:
  - Drug systemic clearance ($\text{CV} = 30\%$).
  - Tumor intrinsic doubling time ($\text{CV} = 25\%$).
  - Baseline tumor-infiltrating lymphocyte density ($\text{CV} = 40\%$).
- Simulate a virtual trial of 50 heterogeneous oncology patients, generating clinical survival outcomes:
  - Progression-Free Survival (PFS) Kaplan-Meier survival curves (defining progression as $> 20\%$ increase in tumor mass over baseline, RECIST 1.1 criteria).
  - Hazard Ratio ($\text{HR}$) and log-rank statistical test comparing Monotherapy vs. Combination.

---

### 3. PUBLICATION-QUALITY DELIVERABLES (300 DPI in `./figures_task7/`)
Write code to compute and save:
1. **`fig1_simeoni_tgi_monotherapy_vs_combo.png`**: Dynamic tumor burden trajectories ($W_{\text{total}}$ vs. Days) comparing Vehicle, Monotherapies, and Combination therapy with standard error bands.
2. **`fig2_ctl_immune_infiltration_dynamics.png`**: Active immune effector cell ($E_{\text{CTL}}$) population expansion and tumor cell death flux over the 60-day window.
3. **`fig3_bliss_synergy_matrix_heatmap.png`**: Dose-response synergy landscape ($\Delta \text{Bliss}$ surface) over a matrix of targeted inhibitor doses vs. anti-PD-1 exposures.
4. **`fig4_virtual_cohort_kaplan_meier_pfs.png`**: Kaplan-Meier progression-free survival (PFS) curves for the 50 virtual patients, featuring median PFS and hazard ratios.

---

### 4. BILINGUAL SCHOLARLY TREATISE & REPOSITORY SEALING
- Author **`QSP_IMMUNO_ONCOLOGY_REPORT_ZH.md`** and **`QSP_IMMUNO_ONCOLOGY_REPORT_EN.md`**:
  - Foundations of Quantitative Systems Pharmacology (QSP) bridging molecular targets to clinical endpoints.
  - Mathematical formulation of transit-compartment tumor cell death kinetics.
  - Mechanistic rationale for kinase-inhibitor and checkpoint-blockade combination: Immunogenic cell death triggering antigen presentation.
- Update root `README.md` to establish the complete 7-Task Core Architecture of the Developability & Translational Suite.
- Automatically stage, commit, and push:
  ```bash
  git add .
  git commit -m "feat(task7): quantitative systems pharmacology (QSP), Simeoni TGI model & combination synergy"
  git push origin main