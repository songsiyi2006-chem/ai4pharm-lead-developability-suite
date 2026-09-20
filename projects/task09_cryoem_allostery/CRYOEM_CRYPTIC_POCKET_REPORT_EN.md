# Task 9 — Public-structure anchored synthetic heterogeneity and mechanical allostery

This executable benchmark combines two public X-ray endpoints with synthetic coordinate interpolation, Gaussian scattering proxies, a five-state Markov generator, geometric cavity measurements and an ANM/PRS calculation. It does **not** establish a new druggable pocket, a cancer-target mechanism, a physical kinetic rate, or a validated allosteric route. The model protein is bacteriophage T4 lysozyme L99A, not a human kinase or GPCR.

## Executed results

| Quantity / 指标 | Executed result / 已执行结果 |
|---|---:|
| Conformers / 合成构象 | 500 |
| Common heavy atoms / 匹配重原子 | 1278 |
| Cα residues / 残基 | 164 |
| Core-aligned Cα endpoint RMSD / 端点偏差 | 0.644236 Å |
| State counts 0–4 / 状态计数 | [225, 129, 79, 61, 6] |
| Estimated open-minus-closed population ΔG | 2.057351 kcal/mol |
| Parametric-bootstrap 95% interval | [1.101320, 3.334307] kcal/mol |
| Assumed-clock closed-to-open MFPT | 2630.063587 ns |
| Conditional opening rate 1/MFPT | 0.000380219 ns⁻¹ |
| MFPT bootstrap 95% interval | [1078.4541766969403, 10896.389935623289] ns |
| Bootstrap draws missing a state / 缺态重抽样 | 22/200 |
| Cavity volume range / 空腔体积 | [76.78125, 209.671875] Å³ |
| Uncalibrated bounded score range / 未校准评分 | [0.45128280616522265, 0.5115793317265658] |
| Joint volume >500 Å³ and score >0.7 / 同时达阈值 | 0 frames |
| ROI boundary-contact frames / 触边帧 | 135 |
| Indeterminate joint-gate frames / 联合阈值未知帧 | 135 |
| Adjacent Cα distance range / 相邻 Cα 距离 | [3.419749282972558, 3.867489756283791] Å |
| Minimum nonadjacent Cα distance / 非相邻最近距 | 4.319542 Å |
| ANM null modes / 零模 | 6 |
| Relative H H⁺ H − H norm / 伪逆残差 | 1.02e-13 |
| PRS force direction convergence / 相对误差 | [(64, 0.069655), (256, 0.036533), (1024, 0.013366)] |
| Contact path residue IDs / 耦合路径 | 106 → 11 |
| Path contour length / 路径长度 | 12.624658 Å |
| Endpoint separation / 端点直线距离 | 12.624658 Å |


No sampled, fully contained component reaches the requested joint cavity gate; 135 frames have unknown full volume and indeterminate gate classification because they touch the ROI boundary. This does not exclude a larger or druggable pocket. The contact path does not span 30 Å end to end in the archived default run. Only 6 open frames were sampled; the generator's declared opening penalty is 1.6 kcal/mol. Sampling noise and explicit pseudocount regularization explain the difference from the estimate. No seed, state energy or score calibration was tuned to force the requested outcome. The estimated ΔG is a **state population difference**, not an activation barrier or a transition-state free energy.

## 9A. Structures, synthetic coordinates and density

Chain A heavy atoms common to 4W51 and 4W59 are matched by residue number and atom name. Alternate A is selected coherently; alternate B, hydrogens, waters, ions and all ligands are excluded from the protein calculation. The bound ligand 3GZ supplies only a cavity-center coordinate and is then removed. Core Cα atoms outside residues 100–120 define a proper Kabsch rotation. The 4W59 alternate-A endpoint is ligand-conditioned; it does not establish a spontaneously populated apo-open structure.

For each generated state s, q=(s+u)/5 with u uniform on (0.001,0.999). Coordinates are (1−q)x_closed+q x_bound plus a tapered random loop displacement, largest around residue 110. The loop amplitude is 0.25 Å times sin(πq). This produces continuous geometric variability and an explicit ordering, but it is neither an atomistic dynamical integrator nor an energy-minimized transition pathway. Adjacent Cα compression to the reported minimum reveals an interpolation artifact; the coarse contact check does not test covalent stereochemistry or all heavy-atom clashes.

Kernel PCA uses the centered RBF kernel of aligned Cα Cartesian distances. Its two latent axes are geometry-derived nonlinear coordinates, not learned cryoDRGN embeddings. Atomic-number weighted coordinate deposits are Gaussian-filtered on a 1.0 Å grid, with **2.5 Å FWHM** point spread. Pixel spacing and PSF width are different quantities. There are no electron micrographs, orientation inference, CTF, detector noise, FSC resolution estimates or Coulomb-potential calibration. Maps average the actual generated snapshots, with state-count weighting for the global mean. Density smoothing/attenuation is an illustrative consequence of coordinate averaging; no experimentally observed loop disappearance is claimed.

## 9B. Reversible synthetic MSM and uncertainty

The declared generator uses nearest-neighbor Metropolis transitions with proposed step probability 0.35/2 per direction and Boltzmann weights from the configured state energies at 298.15 K. Reflecting endpoints are implemented by rejected proposals. The initial state is drawn from the generator's stationary distribution. One saved sample is assigned **10.0 ns by assumption**. These time labels cannot be recovered from unordered cryo-EM snapshots.

From 499 adjacent-frame transitions, the symmetric flux estimator is F=(C+Cᵀ)/2+αS, α=0.5, where S includes diagonal and nearest-neighbor entries. T_ij=F_ij/Σ_jF_ij and π_i=Σ_jF_ij/Σ_ijF_ij exactly satisfy detailed balance. This is a regularized moment estimator, not a reversible maximum-likelihood fit. The PMF is −RT ln(π_i/π_0). MFPT solves (I−T_nonabsorbing)m=τ1 with state 4 absorbing; k_open=1/m_0 is a conditional first-passage rate, not an elementary transition rate.

200 parametric replicate trajectories use the fitted T and are refitted with the same prior. Percentile intervals measure finite-sampling variation under this synthetic model; they exclude force-field, structural and clock uncertainty. Missing states in bootstrap draws remain visible and are regularized, never silently discarded. The 2τ Chapman–Kolmogorov Frobenius discrepancy is 0.131753; this is a finite-data diagnostic, not a successful physical Markovianity test.

## 9C. Geometry, area and the score boundary

The ligand-centered cubic ROI has half-width 9.0 Å. A voxel is accessible when an assumed 1.0 Å probe clears all vdW spheres (C/N/O/S radii 1.70/1.55/1.52/1.80 Å). At least five of six axis directions must encounter a protein obstacle within the ROI. The connected eligible component nearest the ligand centroid defines the seeded cavity. Volume is occupied-voxel count times spacing cubed. Thus this is a local sphere-filling/occlusion algorithm, not a global pocket detector, alpha-shape method or validated commercial score.

Lining atoms lie within their vdW radius +2.5 Å of cavity voxels. Fibonacci-sphere exposure with 96 points and a 1.4 Å probe supplies lining probe-accessible surface area. Enclosed surfaces are included; access from bulk solvent is not proven. Aromaticity is the fraction of lining atoms belonging to Phe/Tyr/Trp/His residues, **not** a measured aromatic surface fraction. Polarity is the N/O/S atom fraction. Enclosure is the mean occluded-direction fraction; hydrophobic depth is maximum component distance-transform radius multiplied by the lining-carbon fraction, an explicitly geometric proxy.

The supplied formula is dimensionally interpreted as raw=0.49 ln[V/(1 Å³)]+0.78 ln(aromaticity)−0.22 polarity−1.2, with a 10⁻⁶ aromaticity floor. It is unbounded. We retain raw and additionally export sigmoid(raw), a bounded **uncalibrated heuristic** with no binding-probability interpretation. Empty cavities receive zero score and null raw. No exact transition threshold can be inferred from these samples. Boundary-contact volumes are ROI-limited descriptions of the detected component and cannot establish a closed pocket volume; full volume is null and classification indeterminate. The 0.5/0.75/1.0 Å endpoint audit measures voxel sensitivity. A separate 9/12/15 Å half-width audit evaluates closed/mid/open representative structures. Since changing ROI can also change directional enclosure and connectivity, a numerically larger component is not automatically the same physical cavity. The per-frame censoring remains, and fine-grid agreement is not global-cavity convergence.

## 9D. ANM, PRS and a contact-constrained hypothesis

The closed-endpoint Cα network uses unit springs and a 13.0 Å contact cutoff. Each pair adds uuᵀ diagonal blocks and −uuᵀ off-diagonal blocks to the 3N×3N Hessian. Six rigid translation/rotation modes are removed with a relative eigenvalue tolerance; all positive internal modes contribute to H⁺. The singular Hessian is never inverted directly. Absolute displacements have arbitrary spring/force units.

For a unit isotropic force at residue j, M_ij=||H⁺_ij||²_F/3. Monte Carlo force sets of [64, 256, 1024] independent directions per source are evaluated via their exact empirical covariance, algebraically identical to averaging squared individual displacement responses. Source self-response normalizes each column. Symmetric coupling uses M_ij/sqrt(M_ii M_jj). Contact edges cost −ln(coupling)+0.01 distance/Å, an explicitly chosen routing metric. Dijkstra selects the cheapest route from any detected lining residue to reference active-site residue 11; this endpoint selection can favor a nearby lining residue. The resulting path is not energy transport, causality, mutational validation or proof of long-range allostery.

## Industry relevance and next acceptance gates

This reusable calculation can expose pipeline failure modes before spending resources on structure-enabled hit discovery: conformational support, cavity truncation, score calibration and pathway sensitivity are separate questions. A real program needs target-specific cryo-EM particles or validated MD, state population/kinetic measurements, side-chain and solvent refinement, unbiased pocket detection, ligandability benchmarks, multiple-structure ANM sensitivity and perturbation experiments. Only then should predicted pockets prioritize fragment screening or medicinal chemistry. The archived negative gates are retained as outcomes rather than converted into discovery claims.

## Primary sources / 原始来源

- [4W51 apo coordinates](https://www.rcsb.org/structure/4W51) and [4W59 n-hexylbenzene-bound coordinates](https://www.rcsb.org/structure/4W59); [Merski et al., 2015](https://doi.org/10.1073/pnas.1500806112). These are X-ray structures of T4 lysozyme L99A.
- [Zhong et al., cryoDRGN, 2021](https://doi.org/10.1038/s41592-020-01049-4): methodological context only; its neural model is not executed here.
- [Trendelkamp-Schroer et al., reversible MSM estimation, 2015](https://arxiv.org/abs/1507.05990): reference for reversibility and uncertainty; this code uses a simpler explicitly stated estimator.
- [Atilgan et al., ANM, 2001](https://doi.org/10.1016/S0006-3495(01)76033-X) and [Atilgan & Atilgan, PRS, 2009](https://doi.org/10.1371/journal.pcbi.1000544).
- [Kuroki et al., T4 lysozyme catalytic-site study](https://pmc.ncbi.nlm.nih.gov/articles/PMC17713/): supports the Glu11 active-site reference; does not validate the computed pathway.
- [wwPDB/RCSB data usage policy](https://www.rcsb.org/pages/usage-policy): deposited PDB coordinate data are CC0; retain structure attribution.

## Inspectable outputs / 可检查产物

- [Coordinate archive](results/synthetic_conformers.npz), [ordered latent trajectory](results/conformer_latent_trajectory.csv), [density voxel maps](results/synthetic_density_maps.npz).
- [MSM estimate](results/msm.json), [bootstrap replicates](results/msm_bootstrap.csv), [pocket geometry](results/pocket_geometry.csv), [grid convergence](results/grid_convergence.json), [ROI sensitivity](results/roi_sensitivity.json).
- [ANM/PRS matrices](results/anm_prs.npz), [path coordinates](results/allosteric_path.csv), [summary](results/summary.json), [numerical verification](verification.json), [byte manifest](manifest.json).
- [Figure 1](figures_task9/fig1_cryoem_latent_conformational_manifold.png), [Figure 2](figures_task9/fig2_msm_free_energy_pathway.png), [Figure 3](figures_task9/fig3_dynamic_pocket_volume_druggability.png), [Figure 4](figures_task9/fig4_prs_allosteric_network_matrix.png).
