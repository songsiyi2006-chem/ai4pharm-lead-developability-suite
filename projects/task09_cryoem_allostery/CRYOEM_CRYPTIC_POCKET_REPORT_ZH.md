# 任务九：公开结构锚定的合成构象异质性与机械变构计算

本项目已执行两端公开 X 射线结构配准、500 个合成坐标构象、非线性嵌入、三维高斯散射代理、五态可逆马尔可夫模型、网格空腔/SASA 及 ANM/PRS。它**没有证明新的可成药口袋、癌症靶点机制、真实动力学速率或已验证变构路径**。模型蛋白是噬菌体 T4 溶菌酶 L99A，不是人源激酶或 GPCR。

## 已执行结果

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


在完整包含于ROI的已采样连通区域中，没有观察到“体积大于500 Å³且评分大于0.7”的联合条件；另有135帧触边，其完整体积和联合阈值分类为未知，**不能排除更大或可成药口袋的存在**。默认路径端点未跨越30 Å。开态仅有6个样本，生成器预先设定的开态代价为1.6 kcal/mol。与估计的差别体现有限采样及显式伪计数正则化，不能宣称发现或验证目标能量区间。未为满足阈值调整随机种子、状态能量或评分参数。这里的ΔG是**状态占有率差异**，并非过渡态活化自由能或真实势垒。

## 9A：结构与连续几何构象

使用4W51和4W59的A链共有重原子，以残基编号和原子名逐一匹配。统一取alternate A构象；不混拼不同占有率原子。排除B构象、氢、水、离子和配体。4W59中的3GZ配体仅提供口袋中心，随后从几何计算删除。使用100–120残基以外的Cα做Kabsch配准，保证旋转矩阵行列式为+1。4W59端点属于配体条件下的结构，不能等同于已证实的无配体自发开态。

生成器先给出状态s，再取q=(s+u)/5，其中u服从(0.001,0.999)均匀分布。坐标为两端线性插值，叠加以110残基为中心逐渐衰减的随机环区位移；最大幅度为0.25 Å×sin(πq)。这使构象几何连续，但不能保证路径动力学可达。相邻Cα最短距离已显示插值压缩，且没有做力场能量最小化、共价立体化学或全重原子冲突审查；因此这些帧不是MD轨迹。

非线性嵌入使用配准Cα坐标距离构建RBF核，再做核PCA，输出两个潜变量。没有训练cryoDRGN，也没有真实单颗粒图像。各原子按原子序数加权投到1.0 Å体素，再用FWHM=2.5 Å的高斯点扩散函数平滑。**FWHM是设置的PSF宽度，不是体素间距，也不是FSC验证的实验分辨率**。未模拟CTF、取向推断、探测噪声或定量库仑势。各状态平均和总体平均均依据实际生成帧数加权，散射单位为原子序数加权代理/Å³。图示只能说明坐标平均造成的平滑和局部衰减，不能声称真实环区密度消失。

## 9B：五态可逆模型与有限采样

在298.15 K下用配置中的五个假设自由能构造Boltzmann权重，以每方向0.35/2的提议概率构造最近邻Metropolis转移，边界拒绝移出状态空间的步长。初态按平衡分布抽取。每帧**人为赋予10.0 ns**；未从结构或冷冻电镜推断时间。

从499对相邻帧计数得到C，定义对称流量F=(C+Cᵀ)/2+αS，α=0.5，S包含对角线和相邻态。归一化行得到T，行和归一化得到π，从代数上满足详细平衡。这是带正则的矩估计，不是假称的可逆最大似然估计。PMF为−RT ln(π_i/π_0)。将状态4吸收，解(I−T_nonabsorbing)m=τ1得到首次到达时间；1/m_0为假设时钟下的条件开启速率。

实际完成200次参数自助法：用拟合T重抽样等长轨迹，并使用同一正则重新拟合。95%百分位区间仅包含该合成模型的采样变化，未包括结构、力场、时间标尺等不确定性。缺态重抽样次数单独保留，不静默删除。2τ的Chapman–Kolmogorov Frobenius差为0.131753，只是有限样本诊断，不能据此宣布真实蛋白满足马尔可夫性。

## 9C：从几何计算空腔、表面积和评分

围绕配体中心建立半宽9.0 Å的局部立方网格。以C/N/O/S=1.70/1.55/1.52/1.80 Å为假设vdW半径，检测1.0 Å探针是否与蛋白相交。候选体素要求六个轴向至少五向在ROI内遇到蛋白；选择距离配体中心最近的候选连通区域。体积为体素数量×间距³。该算法是局部球填充/遮挡方法，不是全蛋白无偏口袋搜索、alpha-shape或商用可成药评分。

距空腔体素小于自身vdW半径+2.5 Å的原子定义衬里。使用96点Fibonacci球与1.4 Å探针计算衬里可接近面积，包含内部封闭表面，不证明与体相溶剂连通。芳香性定义为属于Phe/Tyr/Trp/His残基的衬里原子比例，非精确芳香表面积；极性为N/O/S比例。遮挡比例是六方向平均遮挡，疏水围合深度是连通区域最大距离变换半径×衬里碳比例，均明确为几何代理。

原提示词公式并不在[0,1]内。这里保留raw=0.49 ln[V/(1 Å³)]+0.78 ln(芳香比例)−0.22极性−1.2；对零芳香比例使用10⁻⁶下限，再额外计算sigmoid(raw)。后者仅为**未校准有界描述符**，不等于结合概率或成药概率。空腔为空时raw记录null、评分记0。有限帧数不足以给出精确连续阈值，未观察到达到的门槛不能擅自声明达到。

触及ROI边界的帧不能作为完整闭合口袋体积；完整体积置null、阈值分类标为indeterminate，仅保留ROI内计算值。已在两端完成0.5/0.75/1.0 Å网格比较，另在闭/中间/开三个代表结构上执行9/12/15 Å半宽ROI敏感性。ROI扩大可能同时改变方向遮挡与连通性，因此扩大后的数值并非自动代表同一个物理空腔。逐帧截断标志仍保留；体素结果接近不代表全空腔体积收敛。

## 9D：ANM伪逆与扰动响应

在闭端点Cα上，用13.0 Å接触截断和单位弹簧构造3N×3N Hessian。接触对在对角块加uuᵀ、交叉块减uuᵀ。按相对本征值阈值剔除六个平移/旋转零模，使用其余全部内部模构造H⁺，没有直接对奇异H求逆。弹簧常数和外力幅度未经物理标定，位移仅有任意响应单位。

对源残基j施加单位各向同性力，精确平均平方响应为M_ij=||H⁺_ij||²_F/3。每残基实际抽取[64, 256, 1024]个随机方向，通过经验力协方差计算，代数上等价于逐方向计算位移后求平方均值。闭式结果与随机方向结果比较，最后的相对误差见上表。每列按源自身响应归一化；图中列是施力残基，行是响应残基。

耦合定义为M_ij/sqrt(M_ii M_jj)。只保留接触边，代价为−ln(耦合)+0.01距离/Å，用Dijkstra从全部被检测衬里残基到参考活性位点11选择最低代价路径。该定义可能偏向近邻衬里端点；结果不能冒充跨越30 Å的信号传递，也不表示真实能量流、因果调控或突变验证。必须在报告保留实际端距和路径长度。

## 业界使用场景与验收边界

本项目适合作为结构驱动先导发现前的计算流程试验：分别审查状态支持度、ROI截断、评分校准和机械路径可靠性。要推进真实靶点，应取得原始冷冻电镜颗粒或经验证MD、人口与动力学测量、侧链和水环境优化、无偏口袋检测、片段结合基准，以及多结构/接触参数敏感性和突变实验。只有这些证据逐步补齐后，才能把结果用于片段筛选或药化优先级。本次未满足的门槛已作为负结果保留。

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
