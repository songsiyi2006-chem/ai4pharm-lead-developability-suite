# Phase 1 · Lead developability and multi-parameter optimization

**AI4Pharm Lead Developability Suite · September 2026**

[中文版本](DEVELOPABILITY_MPO_REPORT_ZH.md) · [Model card](MODEL_CARD.md) · [Results](results_task1/developability_results.csv)

## Executive interpretation

This reproducible benchmark evaluates 30 parent compounds in three equally sized clinical archetypes. It connects molecular structure to continuous desirability functions, solubility estimation, exploratory ADMET scenarios, and conformational hydrogen-bond geometry. Its purpose is to expose trade-offs and guide experimental prioritization. It does not establish clinical safety, oral bioavailability, or externally validated prediction performance.

The implemented CNS-MPO median is **4.990** for the oral reference group, **3.647** for toxicity references, and **1.617** for bRo5 modalities. These are calculations on a deliberately selected panel with approximate ionization inputs. The distributions overlap; they do not justify a diagnostic threshold or a general claim about all drugs. In particular, CNS property preferences systematically disadvantage large, polar, peripherally acting molecules.

## 1. From potency to developability

A strong biochemical interaction is only one step toward a useful medicine. Candidate selection also depends on whether a practical formulation produces adequate exposure, whether unbound drug reaches its target, whether clearance and interactions are manageable, and whether the effective exposure is separated from harmful exposure. Optimizing potency alone can move a series toward increased size or lipophilicity while leaving these questions unanswered.

The workflow therefore treats potency, exposure, formulation and safety as connected design constraints. Continuous scoring can preserve intermediate designs that improve one property while temporarily compromising another. A score is a decision aid; the underlying property profile and missing evidence remain visible. The original CNS-MPO publication provides a precedent for this combined-desirability approach. [Wager et al., 2010](https://doi.org/10.1021/cn100008c).

## 2. Panel curation and boundaries

Ten oral references satisfy all four operational Ro5 conditions: MW ≤500 Da, RDKit cLogP ≤5, HBD ≤5 and HBA ≤10. Imatinib and atorvastatin were not included in that strict stratum because their parent molecular weights exceed 500 Da. Sildenafil, gefitinib, aspirin and diazepam preserve requested examples that meet the chosen definition.

The toxicity stratum mixes historical withdrawals, restrictions and investigational failures. Their mechanisms differ. Cerivastatin's rhabdomyolysis, rofecoxib's thrombotic risk and fialuridine's hepatic/mitochondrial toxicity cannot be collapsed into a single hERG label. [FDA cerivastatin review](https://www.accessdata.fda.gov/drugsatfda_docs/nda/2003/21-366_Crestor_Admindocs_P1.pdf), [APPROVe trial](https://pubmed.ncbi.nlm.nih.gov/15713943/), [fialuridine trial](https://pubmed.ncbi.nlm.nih.gov/7565947/).

The bRo5 stratum includes macrocycles, nonmacrocyclic large molecules and a PROTAC. Paclitaxel is a challenging conventional IV comparator; it is not counted as evidence of oral success. ARV-110 is represented by an investigational bavdegalutamide record with unspecified stereochemistry, an explicit limitation for geometric interpretation. [PubChem depositor record](https://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?sid=441604884).

All structures are frozen with PubChem provenance and checked by molecular formula and InChIKey. Salt/mixture inputs and invalid identities are rejected. Archetype labels are excluded from every scoring equation.

## 3. Mathematical formulation

For a property preferred at low values, define

$$L(x;a,b)=\min\left(1,\max\left(0,\frac{b-x}{b-a}\right)\right).$$

Let `U=1−L` and `W(x;a,b,c,d)=min(U(x;a,b),L(x;c,d))`. The six-term CNS-MPO sum uses the following continuous ramps:

$$M=L(P;3,5)+L(D;2,4)+L(MW;360,500)+W(TPSA;20,40,90,120)+L(HBD;0.5,3.5)+L(pK_a;8,10).$$

The sum lies in [0,6]. The original cLogP preference is ≤3; the requested optimum of 2–4 is delivered as a **separately named custom variant**. These transformations are documented in the [independent implementation comparison](https://pmc.ncbi.nlm.nih.gov/articles/PMC8260158/).

RDKit cLogP and chemical-class pKa assumptions replace unavailable original input calculations. A dominant acid/base neutral-fraction approximation provides logD₇.₄. Values are labeled approximate, with optional sourced overrides and a ±1 pKa scenario envelope. A null basic pKa indicates no recognized basic site and contributes one under the documented convention; missed ionization chemistry can invalidate that assumption.

The oral score is `100 × mean(d_solubility, d_permeability, 1−hERG_proxy, d_TPSA, d_Fsp3, d_aromatic_rings)`. It is a custom preference score, not an oral-bioavailability model. Fixed anchors and every component are available in [MODEL_CARD.md](MODEL_CARD.md).

## 4. ADMET and conformational interpretation

The solubility equation is

$$\log_{10}S=0.16-0.63P-0.0062MW+0.066RotB-0.74AP.$$

Here S is mol/L and AP is aromatic heavy-atom fraction. The conversion to µg/mL multiplies `10^logS` by `MW×1000`. This implements published ESOL coefficients with RDKit descriptors; it does not reproduce the original descriptor software or validate extrapolation to macrocycles. [Delaney, 2004](https://doi.org/10.1021/ci034243x).

hERG uses a transparent combination of lipophilicity, assumed cation fraction, aromaticity, basic-nitrogen/ring graph distance, and a neutral-form Gasteiger electrostatic proxy. Caco-2 uses TPSA and charge-weighted volume. HIA uses a sigmoid of a desolvation-energy-like score. **These three formulas are unfitted heuristics**, not validated QSAR estimators; their numerical units or bounded ranges do not establish predictive validity. A low hERG score is not overall safety.

ETKDGv3/MMFF ensembles count internal donor–hydrogen vectors using distance, angle and topological constraints. The range of internal vector counts is reported as a geometric flexibility proxy. Unpaired sites are not proven solvent-exposed; gas-phase sampling cannot establish environment-dependent chameleonicity. The requested absolute index is therefore explicitly null. Convergence counts, retry methods, timeouts, and missing geometry are recorded rather than silently replaced by successful-looking values.

## 5. Terfenadine: exposure turns a liability into clinical harm

Terfenadine illustrates why target activity and adequate absorption do not guarantee an acceptable safety margin. Its parent compound blocks hERG, reducing the repolarizing potassium current and providing a mechanism for QT-related arrhythmia. The original channel study found a much smaller effect for the carboxylate metabolite and implicated parent-drug accumulation. [Roy et al., 1996](https://pubmed.ncbi.nlm.nih.gov/8772706/).

Impaired CYP3A-mediated clearance can increase parent exposure, making the interaction between metabolism and channel activity central to the case. The relevant experimental program needs channel concentration–response measurements, unbound exposure and interaction studies. The repository heuristic only highlights a structural concern; it does not model exposure or reconstruct the clinical withdrawal decision. The broader connection between terfenadine and metabolic interactions was demonstrated in a [clinical interaction study](https://pubmed.ncbi.nlm.nih.gov/8445813/).

## 6. Venetoclax: oral delivery beyond Ro5

The computed parent MW of venetoclax is approximately **868.5 Da**, with three Ro5 violations under this implementation. Its low CNS-MPO result is unsurprising and does not refute its oral therapeutic use: CNS-MPO addresses a different design objective.

Venetoclax's tablet uses an amorphous solid dispersion to address its delivery challenge, and food materially influences exposure. Oral success therefore reflects the compound together with formulation and administration conditions; it cannot be inferred from Ro5 alone. [Human formulation/bioavailability study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9338003/).

It would be unsupported to attribute this success solely to chameleonic hydrogen bonding from the present calculations. The practical lesson is to couple bRo5 design with solubility, precipitation, permeability and formulation experiments. The model's strong lipophilicity penalty and high hERG heuristic for venetoclax illustrate why an uncalibrated score must not become a clinical safety conclusion.

## 7. Figures, verification and next experiments

![MPO comparison](figures_task1/fig1_mpo_distribution_violin.png)

Figure 1 combines violins, boxes and individual observations. Figure 2 uses fixed desirability anchors for six axes; its safety axis covers hERG only. Figure 3 plots cLogP versus MW, with RotB-scaled bubble areas and the custom oral score as color. Its shaded rectangle omits HBD/HBA and is not an oral-success boundary.

The software checks numerical breakpoints, ESOL coefficients and unit conversion, structure identity, duplicate rejection, ionization behavior, absence of label leakage, deterministic small-molecule geometry, timeout handling, CLI validation, image resolution and manifest hashes. [VALIDATION.md](VALIDATION.md) records the delivered run. These checks validate software behavior, not ADMET biology.

An experimental follow-up should measure pKa/logD, kinetic and thermodynamic solubility under specified conditions, bidirectional Caco-2 transport with recovery, hERG concentration–response with exposure margins, metabolic stability and formulation performance. Prospective model development additionally needs independent assay labels, scaffold-separated evaluation, uncertainty calibration and explicit applicability domains. The historical panel is a useful demonstration, not that validation dataset.
