---

### 第三部分：任务十执行指令（靶向 RNA 的小分子药物设计与前体剪接调控热力学）

> **前沿背景：** 全球只有不到 $2\%$ 的基因组编码蛋白质，而超过 $70\%$ 的基因组转录为非编码 RNA。靶向 RNA 的小分子药物（如治疗脊髓性肌萎缩症 SMA 的重磅炸弹药物利司扑兰 Risdiplam）是当今药学界最革命性的全新疆域。任务十构建动态 RNA 发卡/突起基序三级结构系综、小分子配体诱导契合与可变剪接调控的完整生物物理闭环。

```markdown
[TASK 10: RNA-TARGETED SMALL-MOLECULE DRUG DISCOVERY & PRE-MRNA SPLICING MODULATION THERMODYNAMICS]

NOTE ON PRO TIER UTILIZATION:
Deploy full chain-of-thought reasoning without computational shortcuts. Deliver the complete, production-grade, bug-free Python script `run_task10_rna_targeted_small_molecule_dynamics.py` in a single, fully-contained monolithic artifact alongside all figures and reports.

---

### 1. THE REVOLUTION: DRUGGING THE NON-CODING TRANSCRIPTOME
Over $85\%$ of the human genome is transcribed into RNA, but historically less than $0.05\%$ of drug discovery targets RNA directly. Breakthroughs like Risdiplam prove that small molecules can selectively bind dynamic RNA secondary/tertiary structural motifs to modulate alternative pre-mRNA splicing.
Task 10 constructs a state-of-the-art computational biophysics framework modeling:
1. Dynamic RNA secondary structure ensembles and base-pair opening fluctuations.
2. 3D major/minor groove pocket geometry and non-canonical base-pairing.
3. Thermodynamic mass-action shift in alternative splicing equilibrium (Exon Inclusion %).

---

### 2. FOUR-TIER RNA-CADD COMPUTATIONAL PROTOCOL

#### MODULE 10A: Secondary Structure Ensemble & Base-Pair Opening Fluctuations
- Target RNA Motif: A model pre-mRNA exon-intron 5'-splice site junction featuring an asymmetric internal bulge loop (inspired by SMN2 / Risdiplam binding site: $5'\text{-GGCAGU...-3'}$):
  - Evaluate the base-pairing probability matrix $\mathbf{P} \in \mathbb{R}^{L \times L}$ via the McCaskill dynamic programming partition function algorithm:
    $$Q = \sum_{\text{structures } s} \exp\left( -\frac{\Delta G^\circ(s)}{k_B T} \right)$$
  - Calculate the positional Shannon entropy ($H_i = - \sum_j P_{ij} \log_2 P_{ij}$) to identify flexible, non-canonical bulge regions ($H_i > 1.2\text{ bits}$).
- Model the spontaneous base-flipping equilibrium (e.g., adenine or guanine flipping out of the helical stack into the solvent) extracting the unliganded free energy penalty: $\Delta G_{\text{flip}}^\circ \in [4.0, 7.5]\text{ kcal/mol}$.

#### MODULE 10B: 3D RNA Tertiary Groove Modeling & Pocket Featurization
- Construct the 3D atomistic representation of the RNA double-helix containing the internal bulge loop:
  - Enforce canonical A-form RNA helical parameters (deep/narrow major groove, shallow/wide minor groove) alongside non-canonical bulge distortions.
  - Quantify groove dimensions along the helical axis: Major groove width ($d_{\text{major}} \in [3.0, 10.0\text{ \AA}]$) and depth.
  - Featurize the RNA binding cleft: Hydrogen-bonding donor/acceptor arrays from nucleobase edges (Watson-Crick, Hoogsteen, and Sugar faces), electrostatic negative surface potential from phosphate backbones, and inter-strand purine-purine $\pi$-stacking platforms.

#### MODULE 10C: Small-Molecule Induced-Fit Thermodynamic Binding Engine
- Screen a panel of 4 small-molecule chemical scaffolds targeting the bulge junction:
  - Molecule 1: Risdiplam-like planar heteroaromatic derivative (intercalates between non-canonical base pairs).
  - Molecule 2: Aminoglycoside-like highly charged basic ligand.
  - Molecule 3: Inactive mismatch control compound.
  - Molecule 4: Flexible linker-connected bis-intercalator.
- Model the **Induced-Fit Coupled Thermodynamic Cycle**:
  $$\text{RNA}_{\text{closed}} \stackrel{\Delta G_{\text{conf}}^\circ}{\rightleftharpoons} \text{RNA}_{\text{competent}}^* + L \stackrel{\Delta G_{\text{bind}}^\circ}{\rightleftharpoons} [\text{RNA}^* \cdot L]$$
  Calculate the net effective dissociation constant:
  $$K_{D, \text{eff}} = \frac{K_D}{\text{Fraction}_{\text{competent}}} = K_D \left( 1 + \exp\left(\frac{\Delta G_{\text{conf}}^\circ}{k_B T}\right) \right)$$
- Evaluate the specific binding energetics: Stacking free energy ($\Delta G_{\text{stack}}$), hydrogen-bonding network with exocyclic functional groups, and electrostatic phosphate neutralization.

#### MODULE 10D: Alternative Splicing Equilibrium & Clinical Dose-Response Twin
- Couple RNA target engagement to the multi-component spliceosome assembly equilibrium (U1 snRNP binding to the 5' splice site):
  - Unliganded baseline: 5' splice site possesses weak intrinsic U1 snRNP affinity, resulting in exon skipping ($f_{\text{inclusion}} < 10\%$).
  - Ligand-stabilized state: Small molecule $[\text{RNA}^* \cdot L]$ forms a composite structural interface with U1 snRNP protein components (e.g., U1A / U1-C), boosting spliceosome binding affinity by $\Delta\Delta G_{\text{splice}} = -3.5\text{ kcal/mol}$.
- Formulate the non-linear mass-action dose-response equation predicting clinical **Exon Inclusion %**:
  $$\text{Exon Inclusion}([\text{Drug}]) = \text{Baseline} + \frac{\text{Max} - \text{Baseline}}{1 + \left( \frac{\text{EC}_{50}}{[\text{Drug}]} \right)^{\text{Hill}}}$$
- Extract pharmacological endpoints: $\text{EC}_{50}$ (potency), maximum therapeutic exon inclusion percentage ($\text{Efficacy} \ge 85\%$), and therapeutic dosage window avoiding global off-target pre-mRNA missplicing.

---

### 3. PUBLICATION-QUALITY DELIVERABLES (300 DPI in `./figures_task10/`)
Write code to compute and save:
1. **`fig1_rna_secondary_basepair_entropy_map.png`**: Heatmap of the RNA base-pairing probability matrix $\mathbf{P}_{ij}$ alongside positional Shannon entropy $H_i$, identifying the localized flexible bulge loop.
2. **`fig2_rna_3d_groove_electrostatic_cleft.png`**: 3D geometric trajectory of A-form RNA major/minor groove width fluctuations and negative electrostatic phosphate potential profiles.
3. **`fig3_coupled_thermodynamic_binding_cycle.png`**: Free energy decomposition bar chart comparing the 4 ligand scaffolds: Conformational opening penalty ($\Delta G_{\text{conf}}^\circ$) vs. Intrinsic binding energy ($\Delta G_{\text{bind}}^\circ$) vs. Net $K_{D, \text{eff}}$.
4. **`fig4_alternative_splicing_exon_inclusion.png`**: Pharmacodynamic sigmoidal dose-response curve of Exon Inclusion % as a function of drug concentration, comparing wild-type vs. mutated splice sites.

---

### 4. BILINGUAL SCHOLARLY TREATISE & REPOSITORY COMPLETION
- Author **`RNA_TARGETED_CADD_REPORT_ZH.md`** and **`RNA_TARGETED_CADD_REPORT_EN.md`**:
  - The paradigm shift: From targeting the proteome to drugging the transcriptome.
  - Biophysical challenges of RNA-CADD: Extreme structural charge density, non-canonical base-pair plasticity, and water-mediated binding networks.
  - Molecular mechanism of splicing modulation: How small molecules create composite interfaces with macromolecular complexes (spliceosomes).
- Update root `README.md` to index all 10 Tasks, establishing the definitive **10-Task Decade Architecture of the AI4Pharm Suite**.
- Automatically stage, commit, and push:
  ```bash
  git add .
  git commit -m "feat(task10): RNA-targeted small molecule CADD, secondary ensemble & pre-mRNA splicing thermodynamics"
  git push origin main