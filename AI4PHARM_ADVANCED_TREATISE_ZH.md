# AI4Pharm 进阶模块：递送、生物情境与分子设计

**任务 11–20 · 综合技术报告 · 2026 年 9 月 21 日**  
[English version](AI4PHARM_ADVANCED_TREATISE_EN.md) · [任务 1–10](AI4PHARM_DECADE_TREATISE_ZH.md) · [仓库导航](README.md)

本报告汇总本轮八个交付模块的方程、本机实际计算、参数来源与验收要求；Task 11 和 Task 17 按调整后的范围跳过。证据类型并不相同：公开实验结构、来源支持的药物基因组分类，与假设动力学、有限总量模拟、短程原子级试运行和未经校准的分子生成分数同时存在。整合使这些差异可检查，不会自动形成经过临床校准的模型。

研究合规要求落实为可复现输入、正确单位、明确守恒、来源归属、数值检查，以及由实际输出支持的结论。请求中的最佳区间、机制或自由能目标，如果当前计算未证明，就保留为未完成的科学目标。本项目不假设超算资源；多个 Agent 可分工实现和独立审阅有界模块，本地 CPU 任务仍保留明确资源限制。

## 阅读导航

[11 LNP](#task11) · [12 染色质降解](#task12) · [13 BBB–TfR](#task13) · [14 BiTE](#task14) · [15 受体信号](#task15) · [16 转运体](#task16) · [17 肺部递送](#task17) · [18 药物基因组](#task18) · [19 原子级采样](#task19) · [20 分子生成](#task20) · [综合解释](#integration)

<a id="task11"></a>
## 11. LNP 模块——本轮跳过

按用户调整后的范围，本轮不交付 Task 11。本报告不包含该任务的模型、计算或图，仅保留编号便于定位。[任务状态](projects/task11_lnp/README.md)
<a id="task12"></a>
## 12. 染色质情境竞争与靶蛋白降解

### 结合热力学及动力学网络

未结合的游离情境和染色质情境靶点为 $F,C$，降解剂和连接酶为 $P,E$。二元态为 $PF,PC,EP$，三元态为 $EFP,ECP$。纯结合平衡满足

$$PF=pf/K_T,\quad PC=pc/K_T,\quad EP=ep/K_E,$$
$$EFP=\alpha_Fepf/(K_EK_T),\quad ECP=\alpha_Cepc/(K_EK_T),\quad\alpha_C=\alpha_Fe^{-\Delta G_{occ}/RT}.$$

各有限总量包括含该组分的所有物种。消去 $f=F_{tot}/[1+p/K_T+\alpha_Fep/(K_EK_T)]$ 及相应 $c$，再对游离 $p,e$ 有界求根。标称 $K_E=50$ nM、$K_T=100$ nM、$\alpha_F=10$，假设遮挡代价 2 kcal/mol，对应 $\alpha_C=0.38968$。

降解过程采用显式动力学，不假定始终处于平衡。每个可逆反应 $r$ 按 $\dot x=\sum_r\nu_r(v_{r,+}-v_{r,-})$ 组装；二元解离率为 $k_{on}K_D$，三元解离率除以相应 $\alpha_F$ 或 $\alpha_C$。两条路径 $EP+F\rightleftharpoons EFP$ 和 $PF+E\rightleftharpoons EFP$ 因而具有相同平衡乘积；染色质物种同样处理。

仅未结合靶点交换情境，$J_{chrom}=k_{cap}F-k_{rel}C$。所有含靶点物种以 $k_0$ 周转，释放仍存活的降解剂与连接酶；新靶点以 $s=k_0T_0$ 进入 $F$。三元特异降解通量为 $k_{cat}(EFP+ECP)$。求和消去结合与情境交换：

$$\dot T=s-k_0T-k_{cat}(EFP+ECP),\qquad T+D-S=T_0.$$

$D,S$ 为累计降解和合成账本，降解剂与连接酶总量有限且守恒。无药时 $T=T_0$，且

$$C(t)=C_\infty+[C(0)-C_\infty]e^{-\lambda t},\quad\lambda=k_{cap}+k_{rel}+k_0,\quad C_\infty=k_{cap}T_0/\lambda.$$

这个独立解析基线说明染色质比例改变未必意味着药物诱导降解。

### 实际结果与开发接口

100 nM 降解剂处理 48 h 后，总靶点为 **51.2557 nM**，包括 **22.7889 nM 游离情境池** 和 **28.4668 nM 染色质情境池**；每个池都包括药物结合物种。未处理时，30:70 初始分配松弛到 **37.7776:62.2224 nM**，总量并未下降。**152 次动态场景和 152 个结合平衡**覆盖 38 个剂量与4个能量代价，最大网格守恒误差约 **1.29×10⁻⁹ nM**。

全部生物数值参数是假设；没有加载核小体坐标、测量位阻能、表达有限染色质位点或多聚体解离。实验证实的选择性 BET 降解提供研究动机，不提供本模型的参数。[Zengerle 等](https://pmc.ncbi.nlm.nih.gov/articles/PMC4548256/) 下一步应开展匹配的可溶/染色质分级、合成周转、连接酶丰度、洗脱和可及性扰动。能量罚项是可被检验的输入，不是计算得到的核小体解耦自由能。

[参数](projects/task12_chromatin/outputs/config.json) · [结果](projects/task12_chromatin/outputs/summary.json) · [来源](projects/task12_chromatin/sources.md) · [完整报告](projects/task12_chromatin/REPORT_ZH.md)

![任务12：有限连接酶结合与情境相关靶点耗减](figures_advanced/fig12_epigenetic_chromatin_depletion.png)

<a id="task13"></a>
## 13. BBB 的 TfR 转运：亲和力、释放与测量定义

### 有限抗体和受体药量

药量单位 nmol、体积 L，因此药量/体积为 nM。状态包括血浆抗体 $A_p$、游离表面受体 $R_s$、表面复合物 $C_s$、内体受体 $R_e$、内体复合物 $C_e$、内体游离抗体 $A_e$、脑游离抗体 $A_b$。净结合通量为

$$J_s=k_{on}[(A_p/V_p)R_s-K_{D,s}C_s],\qquad J_e=k_{on}[(A_e/V_e)R_e-K_{D,e}(t)C_e].$$

共同内体 pH 闭合为 $pH_e=5.8+1.6e^{-t/(1h)}$，通过 $K_{D,e}=K_{D,s}10^{q(7.4-pH_e)}$ 改变亲和力，假设斜率 $q=0.7$。这是声明的酸敏感结合假说，不是实测组氨酸依赖释放。

定义复合物/受体内吞 $I_C=k_iC_s$、$I_R=k_{ir}R_s$，受体回收 $Q=k_rR_e$，结合/游离载荷降解 $D_C=k_{dc}C_e$、$D_F=k_{df}A_e$，脑/血浆释放 $T_b=k_bA_e$、$T_p=k_pA_e$，清除 $L_p=k_{cl,p}A_p$、$L_b=k_{cl,b}A_b$：

$$\dot A_p=-J_s+T_p-L_p,\quad\dot R_s=-J_s-I_R+Q,\quad\dot C_s=J_s-I_C,$$
$$\dot R_e=I_R-Q-J_e+D_C,\quad\dot C_e=I_C+J_e-D_C,$$
$$\dot A_e=-J_e-D_F-T_b-T_p,\quad\dot A_b=T_b-L_b.$$

三个累计损失分别接收 $L_p$、$D_C+D_F$、$L_b$，它们加抗体库存等于初始抗体量。受体满足 $R_s+C_s+R_e+C_e=R_{tot}$。载荷降解后受体回到内体池，因此当前明确不含受体下调；分析会降解 TfR 的体系时需要增加该机制。

### 暴露定义与实际结果

观测量限定为有限时间和血管外暴露：

$$K_{p,0-72h}=\frac{\int_0^{72h}A_b/V_b\,dt}{\int_0^{72h}A_p/V_p\,dt}.$$

标称比值 **0.00227242**，脑 AUC **6.69954 nM·h**。因为没有脑靶点结合池，这不是所有靶点结合抗体的通用 $K_{p,uu}$ 预测。假设残留血管污染项 $C_pV_{vascular}/V_b$ 会使表观比值变为 **0.0522724**，恰好多 **0.05**；这是测量定义偏差示例，不是脑实质递送改善。

**328 次 ODE** 是亲和力、pH 敏感性、剂量和受体量的单因素切片，而非全笛卡尔网格。标称亲和力峰值为 **13.3352 nM**，不在要求的 50–500 nM 区间。零受体或零复合物内吞产生零脑递送；零受体条件独立满足 $A_p=A_p(0)e^{-0.03t}$。

所有转运动力学数值是假设。特定小鼠体系中降低亲和力促进 BBB 摄取有实验先例，但不是普遍最优定律。[Yu 等](https://pubmed.ncbi.nlm.nih.gov/21613623/) [TfR 转运研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC3920563/) 转化需要物种特异受体丰度、内源转铁蛋白竞争、价态/亲合效应、受体周转、去毛细血管污染的暴露及靶点结合。本模型的价值是暴露释放—滞留及测量定义的权衡。

[参数](projects/task13_bbb_tfr/outputs/config.json) · [结果](projects/task13_bbb_tfr/outputs/summary.json) · [来源](projects/task13_bbb_tfr/sources.md) · [完整报告](projects/task13_bbb_tfr/REPORT_ZH.md)

![任务13：亲和力相关转运与血管外暴露](figures_advanced/fig13_bbb_tfr_transcytosis_profile.png)

<a id="task14"></a>
## 14. BiTE 分子交联与靶细胞损失

### 平衡计量与细胞单位

游离 CD3、肿瘤抗原和 BiTE 为 $e,t,b$，单位 nM。二元态 $EB=eb/K_E$、$TB=tb/K_T$，分子三元态 $X=\alpha etb/(K_EK_T)$。有限总量满足

$$E_0=e+EB+X,\qquad T_0=t+TB+X,\qquad B_0=b+EB+TB+X.$$

嵌套有界求根满足守恒，不把游离和总 BiTE 混同。交换两个受体标签应得到同样三元量，移除一个受体应回到二元二次方程极限，这些构成独立拓扑检查。

假设每初始靶细胞接触体积 $v=10^{-12}$ L，将可及受体拷贝数转换为浓度。令 $n=N/N_0$，效应细胞/初始靶细胞比例为 $r$，拷贝数为 $n_{CD3},n_{TAA}$：

$$E_0=\frac{r n_{CD3}}{N_Av}10^9,\qquad T_0(n)=\frac{n n_{TAA}}{N_Av}10^9\quad[\mathrm{nM}],$$
$$S(n)=X(n)10^{-9}N_Av/n.$$

$S$ 统计每存活靶细胞的分子三元交联，不是细胞间免疫突触数。真实突触群体还需空间接触、相互作用细胞数和装配动力学。声明 $k_{max}=0.12$ h⁻¹、$K_S=500$ complexes/cell、$h=2$ 时，

$$\dot n=-k_{max}\frac{S(n)^h}{K_S^h+S(n)^h}n,\qquad\dot d=-\dot n,\qquad n+d=1.$$

效应细胞数和外部维持的 BiTE 总浓度固定；丢失的靶受体将配体释放回这一维持池。抗体 PK 与不可逆药量损失未包含。

### 实际结果与开发接口

CD3 $K_D=100$ nM、TAA $K_D=1$ nM 时，采样峰值为 **26.1016 nM 总 BiTE**，对应 **572.752 个三元复合物/靶细胞**。100 µM 时仅剩峰值的 **0.144815%**。维持 10 nM 暴露，48 h 后剩 **3.82331%** 初始靶细胞。运行包含 **981 个结合平衡与 270 次裂解 ODE**；零药物、效应细胞或交联协同性时没有裂解。

高剂量钩状效应来自这一可逆有限受体模型，不是普遍的 >1 µM 临床效力下降。亲和力、拷贝数、接触体积和裂解参数均是假设，未包含细胞因子毒性、连续杀伤动力学、耗竭、抗原脱落或运动。原始机制模型支持研究问题，但没有提供本次速率。[Betts 等](https://pmc.ncbi.nlm.nih.gov/articles/PMC6531394/) [细胞级突触群体动力学](https://elifesciences.org/articles/83659) 下一步接口是受体定量、结合动力学、细胞配对成像与匹配浓度/时间的杀伤实验。

[参数](projects/task14_bite/outputs/config.json) · [结果](projects/task14_bite/outputs/summary.json) · [来源](projects/task14_bite/sources.md) · [完整报告](projects/task14_bite/REPORT_ZH.md)

![任务14：分子交联钩状曲线与条件性靶细胞损失](figures_advanced/fig14_bite_synapse_crosslinking_curve.png)

<a id="task15"></a>
## 15. 有限受体动力学校对与细胞因子信号

### 守恒受体池中的驻留时间辨别

抗原类别 $a$ 的外部储库给出伪一阶结合 $k_aR$。五步磷酸化以 $k_p$ 进行；各结合态均以 $k_{off,a}$ 解离并归还受体：

$$\dot C_{a,0}=k_aR-(k_p+k_{off,a})C_{a,0},$$
$$\dot C_{a,j}=k_pC_{a,j-1}-(k_p+k_{off,a})C_{a,j},\quad j=1,\ldots,4,$$
$$\dot C_{a,5}=k_pC_{a,4}-k_{off,a}C_{a,5},\qquad\dot R=\sum_a k_{off,a}\sum_{j=0}^5C_{a,j}-R\sum_a k_a.$$

相加得 $R+\sum_{a,j}C_{a,j}=1$。配体本身是外部储库，靶抗原在 24 h 后下降。稳态 $R=[1+\sum_ak_a/k_{off,a}]^{-1}$，抗原 $a$ 的总结合量为 $k_aR/k_{off,a}$；每步成功修饰的概率为 $k_p/(k_p+k_{off,a})$，故结合受体中终末信号比例为其五次方。这个解析极限独立检查动力学校对拓扑。[McKeithan](https://pmc.ncbi.nlm.nih.gov/articles/PMC41844/)

令 $T,M$ 为 T 细胞/巨噬细胞激活比例，$I_6,N,F,I_1$ 为归一化细胞因子。定义 $q=(C_{target,5}+C_{self,5})/(0.15+C_{target,5}+C_{self,5})$、拮抗剂占据 $B=X/(1+X)$、IL6R 信号 $s=I_6(1-B)/(1+I_6)$。时间单位小时，网络为

$$\dot T=0.22q(1-T)-0.06T,\quad\dot M=[0.35F/(1+F)+0.20s](1-M)-0.10M,$$
$$\dot I_6=0.18T+1.10M-[0.12+0.08(1-B)]I_6,$$
$$\dot N=0.55T+0.60M-0.30N,\quad\dot F=0.90T+0.10M-0.25F,\quad\dot I_1=0.08T+0.80M-0.20I_1.$$

激活比例在零和一的边界导数朝向区间内部。全部系数及浓度单位是假设；拮抗剂在 12 h 前为零，之后为 $X_0e^{-\ln2(t-12)/48}$，求解器对干预边界分段。

### 结果与有意义的终点

假设条件下靶/自身抗原终末态比例为 **59.499**，不是实测构建体选择性。四种 96 h 场景有 **3,844 行**，另有 **65 点**自身抗原解离率扫描。$X_0=0,1,10,100$ 时，IL6R 信号 AUC 为 **75.201、51.158、17.106、6.781 h**，IL6 峰浓度却为 **5.108、5.944、7.186、7.587 归一化单位**。阻断受体相关清除可同时降低信号、升高配体浓度，并不等于停止 IL6 产生。

本模型没有血压、内皮渗漏、休克、器官损伤或临床毒性概率；`homeostasis_restoration_established` 保留空值。无量纲 $C/K_D$ 不是托珠单抗剂量。特定实验体系中的细胞因子研究支持机制区分，不支持人体休克救治结论。[Norelli 等](https://pubmed.ncbi.nlm.nih.gov/29808007/) 下一步需要构建体特异的结合/磷酸化、抗原密度、真实细胞因子单位、巨噬细胞实验及独立检验的拮抗剂暴露模型。当前接口是实验设计与信号假说。

[参数](projects/task15_receptor_signaling/outputs/config.json) · [结果](projects/task15_receptor_signaling/outputs/summary.json) · [来源](projects/task15_receptor_signaling/sources.md) · [完整报告](projects/task15_receptor_signaling/REPORT_ZH.md)

![任务15：受体校对、细胞因子浓度与受体信号](figures_advanced/fig15_receptor_proofreading_cytokine_balance.png)

<a id="task16"></a>
## 16. 既有 PBPK 循环中的转运体限制细胞进入

### 实际复用与避免重复清除

模块导入 Task 4 的 `PBPK.derivative`，保留肺与系统串联循环及肾小球滤过。首先在肝导数与损失账本中同时撤销原肝代谢通量，使新细胞内通路替代原机制，而非叠加第二条重复清除。七个原始药量状态增加有限肝细胞、肾细胞和肠细胞池；继承的分配系数仍是近似。

游离供体浓度 $C_{u,ext}$ 下，肝载体摄取加被动交换为

$$J_H=\sum_{j\in\{1B1,1B3\}}\frac{V_{max,j}C_{u,ext}}{K_{m,j}(1+I/K_i)+C_{u,ext}}+PS_H(C_{u,ext}-C_{u,H}).$$

$C_{u,H}=f_{u,H}A_H/V_H$；饱和项为质量/时间，被动 $PS$ 为体积/时间。竞争抑制提高表观 $K_m$，不改变极限 $V_{max}$。被动项可以反向，不是单向汇。肝细胞方程为

$$\dot A_H=J_H-CL_{met}C_{u,H}-\frac{V_{max,bile}C_{u,H}}{K_{m,bile}+C_{u,H}}.$$

摄取从肝供体空间扣除，代谢和胆汁进入不同账本。肾细胞 $\dot A_R=J_R-J_{urine}$ 使用独立肾摄取/被动项和饱和 P-gp/BCRP 外排；原滤过仅计算一次，肾摄取不冒称肝 OATP。

移除旧的肠腔→肝直接吸收捷径。令肠腔量 $A_g$、肠细胞量 $A_E$、返回管腔的外排 $J_{eff}$、基底侧速率 $k_b$ 和转运离开速率 $k_f$：

$$\dot A_E=k_aA_g-J_{eff}-k_bA_E,\qquad\dot A_g=-k_aA_g+J_{eff}-k_fA_g.$$

基底侧输出进入肝，管腔离开进入粪便损失。每个交换在供体与受体上具有相反符号，故 $\sum A_i+L_{met}+L_{urine}+L_{feces}=Dose$。AUC 的单位为浓度×时间，不参与质量求和。完整实现见 [TransportPBPK.derivative](projects/task16_transporters/driver.py)。

### 实际暴露与限制

8 种详细 IV/口服场景在100 mg比较 $I/K_i=0,1,10,100$，另 **24 种**剂量/抑制组合检查饱和。IV AUC0–96 为 **7.266、10.250、16.399、18.687 mg·h/L**，AUCR 为 **1、1.411、2.257、2.572**；口服 AUCR 为 **1、1.495、2.533、2.919**。IV $I/K_i=100$ 时96 h仍有 **55.72 mg** 留在系统，因此这些是截断 AUC 比，不能写成完整暴露比。

转运体动力学和新增细胞体积全为假设，固定抑制剂比值没有致相互作用药 PK。增加有效细胞池并未建立经过解剖验证的体积拆分。数值检查覆盖计量、空供体、竞争及加密。产业使用还需转运 $V_{max}/K_m$、细胞内结合、丰度缩放和外部 PK。转运体相互作用属于开发评估范围，但指南本身不验证当前场景。[FDA/ICH M12](https://www.fda.gov/regulatory-information/search-fda-guidance-documents/m12-drug-interaction-studies)

[参数](projects/task16_transporters/outputs/config.json) · [结果](projects/task16_transporters/outputs/summary.json) · [来源](projects/task16_transporters/sources.md) · [完整报告](projects/task16_transporters/REPORT_ZH.md)

![任务16：饱和摄取、外排与有限时间PK](figures_advanced/fig16_transporter_oatp_pgp_kinetics.png)

<a id="task17"></a>
## 17. 肺部递送模块——本轮跳过

按用户调整后的范围，本轮不交付 Task 17。本报告不包含该任务的模型、计算或图，仅保留编号便于定位。[任务状态](projects/task17_pulmonary/README.md)
<a id="task18"></a>
## 18. 药物基因组分类与母药—代谢物暴露

### 来源支持的类别不等于普适清除规律

CYP2D6 活性分数分类为：AS=0 是 PM，$0<AS<1.25$ 是 IM，$1.25\le AS\le2.25$ 是 NM，AS>2.25 是 UM。已定相等位基因拷贝数先乘对应活性再求和；支持的子集中 *10、*41 均赋值0.25。未知等位基因返回不确定。这是有限分类器，不是完整测序/CNV 解释服务。[CPIC 2026 更新，表1](https://files.cpicpgx.org/data/guideline/publication/ondansetron/2026/41979467.pdf)

CYP2C19 按自身等位基因功能分类：*1/*17 是快速，*17/*17 是超快，*2/*17 是中间代谢；功能降低组合保留适当的“可能”类别。不能套用 CYP2D6 分数边界。来源中的特定药物建议也没有转为通用剂量指导。[CPIC CYP2C19 更新](https://files.cpicpgx.org/data/guideline/publication/clopidogrel/2022/35034351.pdf)

### 母药等价摩尔守恒与精确传播

药物库 $D$、母药 $P$、代谢物 $M$ 与已消除母药等价量 $L$，单位 nmol，满足

$$\dot D=-k_aD,\qquad\dot P=Fk_aD-(CL_{other}+CL_{CYP})P/V_P,$$
$$\dot M=yCL_{CYP}P/V_P-CL_MM/V_M,$$
$$\dot L=(1-F)k_aD+[CL_{other}+(1-y)CL_{CYP}]P/V_P+CL_MM/V_M.$$

四项导数相加为零，每次给药增加1000 nmol。账本是母药等价摩尔量，不要求母药和代谢物分子量相同。另外两态积分各自浓度。线性六状态系统在七次每24 h给药之间用矩阵指数传播，168 h终点在下一次给药之前。

假设途径为 $CL_{CYP}=10f_mg/(1+I/K_i)$ L/h，$CL_{other}=10(1-f_m)$，默认 $f_m=0.8$。PM/IM/NM/RM/UM 的因子 **0/0.35/1/1.4/1.8** 是假设底物参数，独立于活性分数算术。抑制改变功能，不重写遗传基因型。

### 已执行发现与下一步证据

9 个双倍型乘3个抑制条件得到 **27 种场景、18,171 行轨迹**，另有 **84 点**途径贡献与活性扫描。未受抑制 CYP2D6 PM/IM/NM/UM 的母药 AUC0–168 为 **2700.647/1217.852/594.118/362.775 nM·h**，活性代谢物为 **0/663.207/934.816/1032.604 nM·h**。母药蓄积和无法生物活化可能同时发生；仅凭代谢表型不能给出通用毒性或疗效变化方向。

暴露归一化倍数 $AUC_{reference}/AUC_{scenario}$ 来自模型线性。当代谢物生成为零，没有有限倍数能匹配非零参考代谢物暴露，所以保留空值。即使倍数有限，也不代表疗效或安全性匹配；临床剂量和死亡概率均为空。

产业转化需要真实药物及活性物种、可靠双倍型/CNV 判读、底物特异的基因型 PK、合并用药、疾病效应和相应药物指南。矩阵数值精度及正确表型标签都不能验证假设 PK 因子。

[参数](projects/task18_pgx/outputs/config.json) · [结果](projects/task18_pgx/outputs/summary.json) · [来源](projects/task18_pgx/sources.md) · [完整报告](projects/task18_pgx/REPORT_ZH.md)

![任务18：有来源的类别与假设母药/代谢物暴露](figures_advanced/fig18_pgx_cyp2d6_phenotype_kinetics.png)

## 19. 原子级三元界面采样

5T35 的 BRD4–MZ1–VHL/Elongin B/C 体系完成 150,838 原子、0.2 ps 无偏热化与 1.0 ps、20 个 hills 的 CPU 试运行。坐标和偏置有限且已审计，但短轨迹不收敛 PMF；协同性 α、ΔΔG 和生产性降解构象保留为空。[Task 19 报告](projects/task19_metadynamics/REPORT_ZH.md) · [图](figures_advanced/fig19_ternary_metadynamics_pmf_landscape.png)

## 20. 三维口袋条件分子生成

基于 6YB7 Mpro 口袋完成 100 个候选、20 代 NSGA-II、等预算随机对照、12 个 Vina 子集和两次 X77 重对接。几何、QED、MPO 和 SA 是模型排序量；重对接 RMSD 为 1.242 Å 和 1.055 Å，但没有把分数转为 nM，也没有声称新颖性、可合成性或临床效力。[Task 20 报告](projects/task20_denovo/REPORT_ZH.md) · [图](figures_advanced/fig20_denovo_pareto_lead_optimization.png)
