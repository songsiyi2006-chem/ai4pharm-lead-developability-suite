[TASK 8: ANTIBODY-DRUG CONJUGATE (ADC) MULTI-SCALE ENGINEERING: DAR DYNAMICS, LYSOSOMAL CLEAVAGE & BYSTANDER DISPERSION]

NOTE ON PRO TIER UTILIZATION:
Deploy full chain-of-thought reasoning without computational shortcuts or toy simplifications. Deliver the complete, production-grade, bug-free Python script `run_task8_adc_dar_cleavage_bystander_dynamics.py` in a single, fully-contained monolithic artifact alongside all figures and reports.

---

### 1. THE ADC TRANSLATIONAL COMPLEXITY
Antibody-Drug Conjugates (ADCs) combine the exquisite targeting specificity of monoclonal antibodies (mAb) with the potent cytotoxic killing of small-molecule payloads via chemical linkers. 
Task 8 constructs an end-to-end multi-scale biophysical and pharmacological framework capturing:
1. Stochastic chemical conjugation kinetics determining the Drug-to-Antibody Ratio (DAR) distribution.
2. In-silico Hydrophobic Interaction Chromatography (HIC) analytical deconvolution.
3. Receptor-mediated internal trafficking, lysosomal Cathepsin-B cleavage, and 3D tumor bystander payload diffusion.

---

### 2. FOUR-TIER INTEGRATED ADC PHARMACOLOGICAL ENGINE

#### MODULE 8A: Stochastic Chemical Conjugation Kinetics & DAR Distribution
- Model the interchain disulfide reduction of an IgG1 antibody followed by cysteine-maleimide conjugation:
  - An IgG1 possesses 4 interchain disulfide bonds yielding up to 8 reactive thiol groups upon TCEP/DTT reduction.
  - Formulate a stochastic chemical kinetic scheme (Monte Carlo Gillespie / master differential equations) tracking the sequential addition of payload-linker intermediates:
    $$\text{mAb} \xrightarrow{k_1} \text{DAR}_1 \xrightarrow{k_2} \text{DAR}_2 \dots \xrightarrow{k_8} \text{DAR}_8$$
    incorporating steric hindrance and accessibility penalties for high-order additions ($k_{n+1} < k_n$).
- Evaluate the emergent binomial/Poisson-like DAR species fractions ($f_{\text{DAR}_0}$ through $f_{\text{DAR}_8}$) and calculate the Average Weighted DAR:
  $$\langle\text{DAR}\rangle = \sum_{i=0}^8 i \times f_{\text{DAR}_i}$$

#### MODULE 8B: Analytical Hydrophobic Interaction Chromatography (HIC) Twin
- Simulate the clinical-grade HIC analytical separation used in quality control:
  - Map each DAR species to an effective hydrophobic retention surface area ($A_{\text{hydrophobic}} \propto \text{DAR}^{1.15}$).
  - Model the descending linear salt-gradient elution (e.g., $1.5\text{ M}$ to $0\text{ M}$ Ammonium Sulfate) coupled with Gaussian-Lorentzian chromatographic band broadening.
  - Quantify the analytical resolution ($R_s$) and purity of the clinically optimal $\text{DAR}_4$ fraction.

#### MODULE 8C: Intracellular Receptor Trafficking & Lysosomal Cleavage Kinetics
- Construct the system of coupled non-linear ODEs tracking the cellular life-cycle of the ADC in a target tumor cell expressing surface antigen (e.g., HER2 or TROP2, density $B_{\max} = 10^6\text{ receptors/cell}$):
  - Surface receptor binding ($k_{\text{on}}, k_{\text{off}}$) and dynamic receptor internalization ($k_{\text{int}} = 0.05\text{ h}^{-1}$).
  - Endosomal sorting and lysosomal degradation: Cleavable dipeptide linker (Val-Cit) undergoes enzymatic cleavage driven by lysosomal Cathepsin-B:
    $$-\frac{d[\text{ADC}]_{\text{lyso}}}{dt} = \frac{V_{\max, \text{cleave}} [\text{ADC}]_{\text{lyso}}}{K_{m, \text{cleave}} + [\text{ADC}]_{\text{lyso}}}$$
  - Release rate of active, uncharged, membrane-permeable free payload ($[D]_{\text{free}}$) into the intracellular cytoplasm.

#### MODULE 8D: 3D Radial Tumor Tissue Bystander Diffusion Dynamics
- Model the **Bystander Killing Effect**: Free payload molecules diffuse out of the dying target-positive cell ($Ag^+$) through interstitial fluid to eradicate neighboring antigen-negative tumor cells ($Ag^-$):
  $$\frac{\partial C_{\text{payload}}(r, t)}{\partial t} = D_{\text{eff}} \frac{1}{r^2} \frac{\partial}{\partial r} \left( r^2 \frac{\partial C_{\text{payload}}}{\partial r} \right) - k_{\text{clear}} C_{\text{payload}} - \text{CellularUptake}(C_{\text{payload}})$$
- Simulate dynamic radial cytotoxic concentration contours across a micro-tumor cord ($r \in [0, 200\text{ }\mu\text{m}]$):
  - Track killing efficacy in $Ag^+$ cells vs. $Ag^-$ bystander cells over 72 hours.
  - Demonstrate how membrane permeability of the payload (e.g., Dxd vs. polar MMAF) governs the bystander therapeutic radius.

---

### 3. PUBLICATION-QUALITY DELIVERABLES (300 DPI in `./figures_task8/`)
Write code to compute and save:
1. **`fig1_adc_conjugation_dar_distribution.png`**: Multi-panel plot: Time-evolution of discrete DAR species ($0$ to $8$) and the final weighted DAR distribution histogram for varying payload equivalents.
2. **`fig2_analytical_hic_chromatogram_twin.png`**: Simulated High-Performance Hydrophobic Interaction Chromatography (HIC) chromatogram resolving all individual DAR peaks ($0, 2, 4, 6, 8$) with peak integration baselines.
3. **`fig3_intracellular_lysosomal_release_ode.png`**: Multi-compartment dynamic concentration curves: Cell-surface bound ADC, internalized endosomal ADC, lysosomal cleavage flux, and free cytoplasmic payload accumulation over 48 hours.
4. **`fig4_bystander_killing_spatiotemporal_contour.png`**: 2D/3D spatial heatmap showing the radial concentration gradient of bystander payload diffusing from antigen-positive cores into surrounding antigen-negative tissue, with labeled cytotoxic threshold boundaries ($C > \text{IC}_{50}$).

---

### 4. BILINGUAL SCHOLARLY TREATISE & REPOSITORY UPGRADE
- Author **`ADC_TRANSLATIONAL_ENGINEERING_REPORT_ZH.md`** and **`ADC_TRANSLATIONAL_ENGINEERING_REPORT_EN.md`**:
  - Biophysical principles of site-specific vs. stochastic conjugation in ADC developability.
  - The mathematics of mass-action and diffusion in tumor microenvironments governing the bystander effect.
  - Linker stability trade-offs: Preventing systemic premature payload leakage while ensuring rapid lysosomal enzymatic cleavage.
- Update root `README.md` to index Task 8.
- Automatically stage, commit, and push:
  ```bash
  git add .
  git commit -m "feat(task8): ADC stochastic DAR dynamics, HIC analytical twin & bystander spatial diffusion"
  git push origin main