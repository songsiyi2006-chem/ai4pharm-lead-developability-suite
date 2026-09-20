# AI4Pharm：药物转化可开发性的十个计算维度

**综合学术报告 · 2026 年 9 月 21 日**

[English version](AI4PHARM_DECADE_TREATISE_EN.md) · [仓库与运行入口](README.md) · [十项项目档案](projects/)

文件名中的 “Decade” 对应十项模块，并不表示十年的研究积累。本报告统一整理各模块的数学基础、已归档计算、证据边界和实验开发接口，定位为**可在本地计算资源上复现的研究与实验优先级分析工具**。本次没有临床验证、新增药理实测、经过确认的生产工艺，也没有构建一个已经统一校准的跨尺度数字孪生。

## 阅读导航

[研究架构](#architecture) · [1 分子可开发性](#task1) · [2 靶向降解](#task2) · [3 共价动力学](#task3) · [4 PBPK](#task4) · [5 CYP 相互作用](#task5) · [6 ASD 制剂](#task6) · [7 肿瘤免疫 QSP](#task7) · [8 ADC 工程](#task8) · [9 结构与变构](#task9) · [10 RNA 药物设计](#task10) · [转化与验收](#translation)

<a id="architecture"></a>
## 研究架构：本次整合了什么

先导化合物的价值来自化学性质、递送、暴露、靶点作用与生物响应之间的配合。这些问题使用不同状态变量，跨越不同时间尺度。分子满意度、三元复合物浓度、肿瘤质量和剪接概率具有不同含义，不能直接相加成为一个具有物理解释的“总成药分数”。

统一入口 [run_ai4pharm_omnibus_suite.py](run_ai4pharm_omnibus_suite.py) 组织已有模块及可检查结果，完成的是软件与证据层面的整合。**Task 5 已实际调用 Task 4 的 PBPK 导数函数。** 其余跨模块关系仍是待开发接口，需要一致的化合物身份、单位、实验条件和校准依据。

| 任务 | 状态变量或核心对象 | 已归档计算 | 支持开发决策仍需的证据 |
|---|---|---|---|
| 1 | 分子结构、描述符、满意度 | 30 个整理后的母体结构 | 实测 pKa/logD、渗透、毒理与独立标签 |
| 2 | 游离、二元和三元复合物 | 九个数量级配体扫描、200 个连接子样本 | 亲和力、有效泛素化与靶蛋白周转 |
| 3 | 游离、结合和共价修饰靶点 | 8 个假设弹头的动力学 | 特定条件下的结合与化学反应速率 |
| 4 | 7 个药物质量状态 | 静注、口服、重复给药与周期解 | 清除、分配、结合与观测 PK |
| 5 | PBPK 与 6 个酶池 | 24 组配对动态 AUCR 比较 | 抑制、周转、代谢物与外部 DDI 数据 |
| 6 | 相组成与 4 个质量池 | 324 行相图、69 次积分 | API/聚合物身份、相行为与溶出测量 |
| 7 | PK、肿瘤损伤与 CTL | 50 个配对虚拟个体、4 种方案 | 肿瘤与免疫纵向观测、剂量反应 |
| 8 | DAR、细胞与组织药量 | 24 次随机重复、空间输运 | 分 DAR 分析、载荷释放与空间活性 |
| 9 | 坐标、构象状态、响应矩阵 | 500 个合成构象、5 状态 MSM | 目标特异的构象、动力学与可结合性 |
| 10 | RNA 配分函数与结合状态 | 14 次折叠、976 个平衡条件 | RNA 结合、剪接报告与转录组选择性 |

全文区分**文献事实、实际执行的计算、假设参数、模型推论、需要实验验证**五类信息。公开结构是已沉积坐标的实验依据；以它为起点生成的轨迹仍是合成轨迹。精确求解方程说明软件正确实现了模型，不能说明模型已经描述了真实患者。缺失、无效、删失与不可辨识的结果均不替换为零。

各章使用自己的单位约定。跨模块传递前必须转换：$1\ \mathrm{nM}=10^{-9}\ \mathrm{M}$，$1\ \mathrm{\mu M}=10^{-6}\ \mathrm{M}$，$C[\mathrm{nM}]=10^6C[\mathrm{mg/L}]/MW[\mathrm{g/mol}]$。摩尔自由能使用 $RT$，单分子能量使用 $k_BT$；对数与指数的自变量必须无量纲。总浓度不能未经说明替代游离浓度。

<a id="task1"></a>
## 1. 分子可开发性、MPO 与超五规则化学空间

### 从偏好函数构造可解释评分

本项目包含口服参考、毒性相关参考和 bRo5 三组，每组 10 个结构。结构身份经过核对，组别标签不进入计算公式。对于“越低越好”的性质，定义

$$L(x;a,b)=\min\{1,\max[0,(b-x)/(b-a)]\},\qquad U=1-L,$$
$$W(x;a,b,c,d)=\min\{U(x;a,b),L(x;c,d)\}.$$

实现的 CNS-MPO 为

$$M=L(P;3,5)+L(D;2,4)+L(MW;360,500)+W(TPSA;20,40,90,120)+L(HBD;0.5,3.5)+L(pK_a;8,10).$$

每项在 $[0,1]$ 内，因此总分在 $[0,6]$ 内。平台与斜坡表达性质偏好，不是概率。CNS-MPO 的原始目标是中枢药物发现；直接用其排序大型外周药物会改变原来的决策问题。另行输出的自定义 cLogP 版本也不能冒充原始评分。[Wager 等，CNS-MPO](https://doi.org/10.1021/cn100008c)。

单一主要碱性位点的中性比例为 $f_N=(1+10^{pK_a-pH})^{-1}$；单一主要酸性位点为 $f_N=(1+10^{pH-pK_a})^{-1}$。若离子形态进入辛醇相的比例可忽略，则 $D\approx Pf_N$，所以 $\log D\approx\log P+\log f_N$。这说明假设 pKa 如何传递到评分，也说明多微观态与离子分配会使该近似失效。

ESOL 方程为

$$\log_{10}(S/\mathrm{M})=0.16-0.63P-0.0062MW+0.066RotB-0.74AP,$$

其中 $AP$ 是芳香重原子比例。单位转换为 $S_{\mu g/mL}=10^{\log S}MW\times1000$。这里采用发表的系数与 RDKit 描述符；系数相同不表示完全复现原论文的软件预测，更不自动保证对大环的外推有效。[Delaney，ESOL](https://doi.org/10.1021/ci034243x)。

### 计算结果及其边界

归档结果中，口服、毒性、bRo5 三组 CNS-MPO 中位数分别是 **4.990、3.647、1.617**。这是特意选择的样本集上的重叠分布，不是经过外部验证的分类器。六轴口服雷达使用另一套自定义偏好：

$$S_{oral}=100\operatorname{mean}(d_{sol},d_{perm},1-hERG_{proxy},d_{TPSA},d_{Fsp^3},d_{aromatic}).$$

hERG、Caco-2 和 HIA 三项是**没有拟合实验标签的启发式公式**。数值范围、单位或图形分组不能赋予其心脏风险概率、真实渗透系数或人体吸收百分比的含义。$Fsp^3$ 和芳环数是结构描述符；ETKDG/MMFF 分子内氢键几何是算法构象代理，不能证明溶剂依赖的“变色龙性”。失败的几何结果和缺乏依据的绝对变色龙指数保留缺失。

**参数来源与下一步决策。** 结构具有冻结的公开来源；电离和部分 ADMET 变换依赖假设。产业接口是化合物登记、性质权衡和实验优先级表。下一步应测 pKa/logD、明确介质条件的溶解度、带回收率的双向渗透、代谢稳定性及 hERG 浓度反应。建立预测器还需要按骨架分离的外部评估与适用域。维奈克拉的低 CNS 分数不否定其外周口服用途。

[逐分子结果](projects/task01_lead_developability/results_task1/developability_results.csv) · [模型卡](projects/task01_lead_developability/MODEL_CARD.md) · [项目完整报告](projects/task01_lead_developability/DEVELOPABILITY_MPO_REPORT_ZH.md)

![任务1：MPO分布、自定义满意度雷达与化学空间](figures_omnibus/fig_task1_developability_mpo.png)

<a id="task2"></a>
## 2. 三元平衡、协同性与靶向蛋白降解

### 由守恒方程得到精确平衡与峰值

以 $p,e,t$ 表示游离 PROTAC、E3 与靶蛋白，$C=[EPT]$ 表示三元复合物，二元解离常数为 $K_E,K_T$，热力学协同性为 $\alpha>0$：

$$[EP]=ep/K_E,\quad[PT]=tp/K_T,\quad C=\alpha ept/(K_EK_T),$$
$$P_0=p+[EP]+[PT]+C,\quad E_0=e+[EP]+C,\quad T_0=t+[PT]+C.$$

两条装配路径的条件解离常数分别是 $K_T/\alpha$ 与 $K_E/\alpha$，以保持热力学循环闭合。消去游离蛋白后，

$$q(p)=\frac{(K_E+p)(K_T+p)}{\alpha p},\qquad Cq=(E_0-C)(T_0-C).$$

将二次方程较小根有理化，避免相近数相减的数值误差：

$$C=\frac{2E_0T_0}{E_0+T_0+q+\sqrt{(E_0-T_0)^2+2q(E_0+T_0)+q^2}}.$$

总配体与游离配体之间满足

$$P_0=p+(E_0-C)\frac{p}{K_E+p}+(T_0-C)\frac{p}{K_T+p}+C.$$

因为 $q=(p+K_E+K_T+K_EK_T/p)/\alpha$，其导数在 $p_* =\sqrt{K_EK_T}$ 为零，对应

$$P_{0,*}=\sqrt{K_EK_T}+\frac{E_0\sqrt{K_T}+T_0\sqrt{K_E}}{\sqrt{K_E}+\sqrt{K_T}}.$$

此处两个二元占据分数相加为一，使 $P_0$ 中 $C$ 的系数消去。因此在**这个指定反应网络内**，最优总浓度与 $\alpha$ 无关，但峰高与峰宽有关。任务原文给出的 $\sqrt{(K_E+E_0)(K_T+T_0)}$ 不是一般恒等式。

当 $p\to0$，$C\sim\alpha E_0T_0p/(K_EK_T)$；当 $p$ 足够大，$C\sim\alpha E_0T_0/p$。因此有限正协同性的双二元结合臂模型最终出现钩状效应。这不能推广为所有分子胶的普遍定律，也不能保证实验可达浓度内一定看到下降段。

### 已执行的平衡与降解场景

在 $E_0=T_0=100$ nM、$K_E=10$ nM、$K_T=100$ nM 下，精确峰值总浓度是 **131.623 nM**，原指定公式则给 **148.324 nM**。当 $\alpha=0.01,1,10,100$，三元峰值分别为 **0.570646、29.0536、66.1477、87.6755 nM**。各自峰值 80% 的区间是曲线形状指标，并非临床剂量窗口。

靶蛋白周转扩展为

$$\dot T_{total}=-V_{max}\frac{C}{K_m+C}+k_{base}(T_{initial}-T_{total}).$$

$V_{max}$ 的单位是浓度/时间，不是一阶速率常数。靶点逐渐减少时重新求平衡，隐含结合快于降解、E3 总量固定及 PROTAC 催化回收。有效泛素化和蛋白酶体处理被归入有效通量。默认 $k_{base}=0$ 也不是已经测得细胞没有再合成。

两种封端连接子代理各有 100 个 ETKDG/MMFF 样本。PEG 与刚性炔基连接子的平均连接点距离分别为 **9.6715 Å、9.6050 Å**。它们不是完整 PROTAC 的蛋白间距离。算法样本频数缺乏平衡布居权重时，不能用 $-RT\ln p$ 直接计算结合熵，更不能据此确定 $\alpha$。三元结构协同性有实验先例，但本项目中的数值仍是假设。[Gadd 等的三元复合物结构研究](https://doi.org/10.1038/nchembio.2329)。

**参数来源与下一步决策。** 默认亲和力、周转和蛋白浓度均为场景参数。应测二元/三元结合、细胞游离暴露、降解与恢复时间曲线、E3/靶点丰度及泛素化能力。真正的开发问题是亲和力、几何或周转哪个因素限制了已测降解剂，而不是把合成平衡峰值命名为人体最佳剂量。

[平衡指标](projects/task02_tpd/data_task2/equilibrium_metrics.csv) · [连接子样本](projects/task02_tpd/data_task2/linker_conformers.csv) · [项目完整报告](projects/task02_tpd/TPD_TERNARY_COOPERATIVITY_REPORT_ZH.md)

![任务2：三元钩状曲线、亲和力敏感性与连接子几何](figures_omnibus/fig_task2_tpd_hook_effect.png)

<a id="task3"></a>
## 3. 共价抑制、驻留时间与 GSH 竞争

### 完整动力学及其约化条件

游离靶点 $F$、非共价复合物 $B$ 与共价复合物 $C$ 的质量作用方程为

$$\dot F=-k_{on}IF+k_{off}B,$$
$$\dot B=k_{on}IF-(k_{off}+k_{inact})B+k_{rev}C,\qquad\dot C=k_{inact}B-k_{rev}C.$$

求和得到 $\dot F+\dot B+\dot C=0$。不可逆化学、游离抑制剂恒定且快速预平衡时，$B/(F+B)=I/(K_D+I)$，其中 $K_D=k_{off}/k_{on}$。令未共价修饰靶点 $U=F+B$，则 $\dot U=-k_{obs}U$，

$$k_{obs}\approx\frac{k_{inact}I}{K_D+I}.$$

快速平衡是条件而非恒真事实。准稳态约化给出 $K_{I,QSSA}=(k_{off}+k_{inact})/k_{on}$。不约化的不可逆瞬态中，令 $a=k_{on}I$、$s=a+k_{off}+k_{inact}$，精确慢衰减率为

$$\lambda_{slow}=\frac{2ak_{inact}}{s+\sqrt{s^2-4ak_{inact}}}.$$

其低浓度斜率为 $k_{inact}/K_{I,QSSA}$。中间浓度下拟合双曲线并不必然等于该特征值；输出明确区分 $K_D$、$K_{I,QSSA}$ 和拟合 $K_I$。

蛋白周转向 $\dot F$ 添加合成项 $k_{syn}=k_{deg}E_{baseline}$，每个靶点状态添加一阶损失。理想完全洗脱后，不可逆约化模型有

$$F(t)=E_{baseline}-[E_{baseline}-F(0)]e^{-k_{deg}t}.$$

非共价本征驻留 $1/k_{off}$、共价键打开寿命 $1/k_{rev}$ 与蛋白/加合物寿命 $1/k_{deg}$ 回答不同问题。多个结合态相互转化时，平均驻留是首达量 $\tau=\mathbf1^T(-Q)^{-1}p_{bound}(0)$，包含断键、复合物内部再次共价化、解离与蛋白周转。

### 电子结构情景与反应性权衡

八种假设弹头进行了 RDKit/EHT 计算。采用 $\mu=(E_H+E_L)/2$、$\eta=(E_L-E_H)/2$、$\omega=\mu^2/(2\eta)$ 和软度 $1/(2\eta)$ 的约定。LUMO 原子布居是冻结轨道局域描述符，不是已经实验验证的 Fukui 函数。场景能垒由声明的弹头基线和有界 LUMO 修正构成，没有寻找过渡态或完成 IRC。

在 1 M 标准浓度下，双分子 Eyring 场景为 $k_2=(k_BT/h)\exp(-\Delta G^\ddagger/RT)/C^\circ$；再乘假设的 GSH 硫负离子比例，得到表观二阶速率。5 mM 缓冲 GSH 的**不可逆或仅正向捕获**极限下，

$$I(t)=I_0e^{-k_{GSH}[GSH]t},\qquad t_{1/2}=\frac{\ln2}{k_{GSH}[GSH]}.$$

实际 W03/W07/W08 场景包含 GSH 加合物逆反应。令 $a=k_{GSH}[GSH]$、$b=k_{reverse,GSH}$、$p=a/(a+b)$，初始全为游离药物时，加合物浓度满足 $C_{SG}(t)/I_0=p[1-e^{-(a+b)t}]$，$I(t)/I_0=1-C_{SG}(t)/I_0$。正向半衰期 $\ln2/a$ 不同于松弛半衰期 $\ln2/(a+b)$；实际转化 50% 的时间仅在 $p>0.5$ 时存在，为 $-\ln(1-0.5/p)/(a+b)$。$p\le0.5$ 时有限时间内达不到 50%。Pareto 图使用明确标注的正向半衰期，而不是可能不存在的实际半转化时间。

比值 $(k_{inact}/K_D)/k_{GSH}$ 虽然无量纲，比较的是两种反应效率，不是临床安全窗。当前 Pareto 集为 **W02、W06、W08**；用户指定的 30 min 反应性筛查阈值标记 **W03、W04、W05**。这些分类依赖未校准速率和能垒。GSH 结合可能体现解毒；缓冲 GSH 体系既没有模拟细胞内 GSH 耗竭，也不能预测特异质 DILI。[Flanagan 等，GSH 本征反应性](https://doi.org/10.1021/jm501412a)。

**参数来源与下一步决策。** 轨道量来自实际半经验计算，结合速率与能垒基线来自假设。优先测含抑制剂耗竭记录的时间依赖抑制、完整蛋白加合物、结合蛋白周转的洗脱恢复以及相同 pH 下的 GSH 动力学。采用真实观测窗口与参数不确定性后，才能讨论驻留与细胞效力的关系。特定激酶体系中已经证实可逆共价化学能够调节驻留，但这不能直接赋予本项目假设分子同样性质。[Bradshaw 等](https://doi.org/10.1038/nchembio.1817)。

[弹头结果](projects/task03_covalent_kinetics/data_task3/warhead_results.csv) · [洗脱轨迹](projects/task03_covalent_kinetics/data_task3/washout_trajectories.csv) · [项目完整报告](projects/task03_covalent_kinetics/COVALENT_DRUG_KINETICS_REPORT_ZH.md)

![任务3：表观失活、洗脱恢复与条件性反应性Pareto前沿](figures_omnibus/fig_task3_covalent_kinetics.png)

<a id="task4"></a>
## 4. PBPK、IVIVE 与暴露方案筛选

### 清除率与质量守恒的推导

微粒体清除率通过明确结合校正和单位转换放大：

$$CL_{int,u}=\frac{CL_{int,mic}}{f_{u,mic}}\,MPPGL\,m_{liver}\frac{60}{10^6}\quad[\mathrm{L/h}].$$

肝稳态下，流入减流出等于代谢，即 $Q_H(C_{in}-C_{out})=f_{u,b}CL_{int,u}C_{out}$。解得

$$CL_{H,b}=\frac{Q_Hf_{u,b}CL_{int,u}}{Q_H+f_{u,b}CL_{int,u}},\quad f_{u,b}=f_{u,p}/R_b,\quad F_H=1-CL_{H,b}/Q_H.$$

已实现系统包含血液、肝、GI 药物库、肾、脑、肺与剩余组织七个**质量状态**。GI 是管腔药物库，不是灌流肠壁。定义组织 $C_i=M_i/V_i$、静脉血浓度 $C_{v,i}=R_bC_i/K_{p,i}$、血浆观测 $C_p=M_B/(V_BR_b)$、动脉浓度 $C_a=C_{v,lung}$：

$$\dot M_g=-k_aM_g,\qquad\dot M_{lung}=Q_c(R_bC_p-C_a),$$
$$\dot M_i=Q_i(C_a-C_{v,i})\quad(i=kidney,brain,rest),$$
$$\dot M_L=Q_H(C_a-C_{v,L})+k_aF_aF_gM_g-H,$$
$$\dot M_B=\sum_{i=L,kidney,brain,rest}Q_iC_{v,i}-Q_cR_bC_p-CL_{renal,p}C_p.$$

其中 $H=CL_{int,u}f_{u,p}C_L/K_{p,L}$，$CL_{renal,p}=GFRf_{u,p}$。再积分肝消除、肾消除及 $(1-F_aF_g)k_aM_g$ 三个损失项，就得到“体内储存＋累计消除＝累计给药”。肝首过已经通过肝代谢产生，不能再给肝输入重复乘一次 $F_H$。

分配系数采用 **Poulin–Theil 组成形式与近似 logD 替代**，并非完整 Rodgers–Rowland 模型。令 $A(f)=D(f_{nl}+0.3f_{ph})+f_w+0.7f_{ph}$，则 $K_{p,t}=A(t)f_{u,p}/[A(p)f_{u,t}]$。相同 pH 的被动分配、组织组成和结合均包含近似；没有 BBB 转运/通透机制。[Poulin 与 Theil](https://doi.org/10.1002/jps.10005)。

### 稳态必须求解，不能由“第七天”代替

线性情况下，给药间隔内 $\dot M=AM$，静注或口服给药是精确状态跳变 $d$。周期给药后状态满足

$$M_{ss,+}=e^{A\tau}M_{ss,+}+d,\qquad M_{ss,+}=(I-e^{A\tau})^{-1}d.$$

完整暴露为 $AUC_{0,\infty}=e_B^T(-A)^{-1}M_0/(V_BR_b)$。这些矩阵解独立核对 ODE、长尾与周期稳态。可选饱和代谢使用 $H=V_{max}C_{u,L}/(K_m+C_{u,L})$，需要额外给定 $K_m$，此时不能继续假设剂量线性缩放。

合成默认条件给出 $CL_{int,u}=38.88$ L/h、$CL_{H,b}=3.72699$ L/h，同剂量口服/静注 AUC 比为 **86.273%**。100 mg 口服峰浓度为 **0.393085 mg/L**。终末半衰期 **62.5233 h**，口服 AUC 有 **57.95%** 位于 48 h 之后。第七天给药前谷浓度 **1.25124 mg/L**，周期稳态对应值 **1.48126 mg/L**；第七天状态差异 **15.53%**，**未达到 1% 稳态准则**。

在假设剂量网格和暴露约束下，研究筛选选择 **15 mg BID**，游离浓度覆盖率 **90.443%**。临床推荐剂量字段仍为空。假设毒性阈值和 IC90 到人体药效的映射都没有得到验证。

**参数来源与下一步决策。** 化合物明确为合成示例；生理、分配与毒性阈值均为情景输入。需测结合、微粒体/肝细胞清除、转运和有充分长尾采样的口服/静注 PK，再针对明确用途确认模型。模型报告需要说明使用情境和支撑证据，软件测试通过不能代替这一步。[FDA PBPK 报告指南](https://www.fda.gov/regulatory-information/search-fda-guidance-documents/physiologically-based-pharmacokinetic-analyses-format-and-content-guidance-industry)。

[实际输入](projects/task04_pbpk/parameters_used.json) · [结果汇总](projects/task04_pbpk/results_summary.json) · [项目完整报告](projects/task04_pbpk/PBPK_DOSE_PREDICTION_REPORT_ZH.md)

![任务4：给药途径、组织分布和七天累积；第七天不等于稳态](figures_omnibus/fig_task4_pbpk_pharmacokinetics.png)

<a id="task5"></a>
## 5. CYP 抑制、酶恢复与药物相互作用

### 从酶损失连接到动态清除

酶丰度归一化为基线一，方程为

$$\dot E=k_{deg}(1-E)-k_{obs}(t)E,\qquad k_{obs}(t)=\frac{k_{inact}I_u(t)}{K_I+I_u(t)}.$$

恒定抑制剂下，$E_{ss}=k_{deg}/(k_{deg}+k_{obs})$，$E(t)=E_{ss}+[E(0)-E_{ss}]e^{-(k_{deg}+k_{obs})t}$。可逆抑制进一步把功能活性缩放为 $a=E/(1+I_u/K_i)$。完全移除抑制剂后才有 $E(t)=1-[1-E(0)]e^{-k_{deg}t}$；体内残留抑制剂会使这一简单恢复式不再适用。

归档静态筛查实现

$$R_1=1+C_{max,u}/K_i,\qquad R_2=1+\frac{k_{inact}5C_{max,u}}{k_{deg}(K_I+5C_{max,u})}.$$

它们是筛查量，不是 AUCR。实现的肝筛查触发值为 $R_1\ge1.02$、$R_2\ge1.25$，口服肠道筛查另行计算。M12 对游离浓度情境及 TDI 的因子 5 有明确说明。本项目使用模拟的第 14 天周期峰值，但未证明已经达到周期稳态，因此保留暴露输入的局限。[FDA/ICH M12，第 2.1.2 与 7.5 节](https://www.fda.gov/media/161199/download)。

仅含抑制的简化静态模型为

$$AUCR=\frac{1}{a_g(1-F_g)+F_g}\frac{1}{a_hf_m+(1-f_m)}.$$

酶活性下降应使该比值上升。更一般的实现根据受抑制内在清除重新计算肝提取和肠逃逸，显式保留肾清除。内在清除权重不自动等于系统性代谢分数。恒定抑制情况下，解析表达式与 PBPK 积分进行独立对照。

### 已实现耦合与计算发现

Task 5 调用 Task 4 的导数函数，添加肝/肠 CYP3A4、CYP2D6、CYP2C9 六个酶池，同时修改匹配的代谢与损失账本。28 次 BID 给药覆盖 14 天，最后一次给药后再观察十天恢复。口服咪达唑仑 2 mg 和美托洛尔 50 mg 分开模拟，并非实测临床混合探针实验。

六个药名标记的是构造场景。酮康唑、克拉霉素、利托那韦、氟康唑、奎尼丁和近零抑制阿莫西林场景的第 14 天咪达唑仑动态 AUCR 分别是 **11.6746、6.4009、16.9622、2.1647、1.4253、1.0000**，不能作为这些真实药物的临床实测倍数引用。克拉霉素场景的肝/肠 90% 酶恢复时间为 **4.9976/3.3666 天**；利托那韦两个 CYP3A4 酶池十天仍未恢复，记作**右删失**。

来源表收录克拉霉素已发表模型估计：$K_I=5.3$ µM，肝/肠 $k_{inact}=0.4/4$ h⁻¹。这三个数对应的 246 行恒定暴露传递函数补充计算独立于主要假设 PK 场景。[Quinney 等](https://pmc.ncbi.nlm.nih.gov/articles/PMC2812061/)。

**参数来源与下一步决策。** 多数 PK/抑制参数是假设；三个文献数值是模型估计，不是本次测得的游离动力学常数。仅凭活性随时间下降，无法区分血红素配位、血红素损伤或蛋白共价修饰。需要 NADPH/时间/浓度对照、耗竭与游离比例、移除后恢复、代谢物鉴定、诱导以及外部探针 PK。本模型省略诱导、自抑制、抑制性代谢物、转运与遗传表型，适合辅助实验与采样设计，不能给出禁忌或患者停药间隔。

[参数来源](projects/task05_cyp_ddi/results/parameter_provenance.csv) · [动态与静态比较](projects/task05_cyp_ddi/results/dynamic_static_comparison.csv) · [项目完整报告](projects/task05_cyp_ddi/CYP_DDI_KINETICS_REPORT_ZH.md)

![任务5：合成预孵育、动态探针PK与依赖假设的筛查矩阵](figures_omnibus/fig_task5_cyp_ddi_mbi.png)

<a id="task6"></a>
## 6. ASD 热力学、玻璃化转变与过饱和

### 相平衡与自由能曲率

药物质量分数 $w$ 转为体积分数 $\phi=(w/\rho_d)/[(w/\rho_d)+(1-w)/\rho_p]$。在 $N_d=1$ 与假设聚合物链段数 $N_p$ 下，无量纲格点混合自由能为

$$f(\phi)=\frac{\phi}{N_d}\ln\phi+\frac{1-\phi}{N_p}\ln(1-\phi)+\chi\phi(1-\phi).$$

加权 Hansen 近似为 $\chi=V_d[(\Delta\delta_D)^2+0.25(\Delta\delta_P)^2+0.25(\Delta\delta_H)^2]/RT$。因为 MPa·cm³ = J，指定单位下 $\chi$ 无量纲。

二阶导数为 $f''=1/(N_d\phi)+1/[N_p(1-\phi)]-2\chi$，旋节线满足 $f''=0$。双节线则满足公切线条件 $f'(\phi_1)=f'(\phi_2)=[f(\phi_2)-f(\phi_1)]/(\phi_2-\phi_1)$，等价于两相化学势相等。联立 $f''=f'''=0$ 得到

$$\phi_{d,c}=\frac{\sqrt{N_p}}{\sqrt{N_d}+\sqrt{N_p}},\quad\chi_c=\tfrac12(N_d^{-1/2}+N_p^{-1/2})^2.$$

当 $N_p=100$，临界**药物**体积分数是 0.909。这里讨论无定形—无定形相分离；在该模型中单相不等于不会结晶。

### 分子运动、有限剂量与正确的成核指数

Gordon–Taylor 必须使用绝对温度：$T_g=(w_dT_{g,d}+Kw_pT_{g,p})/(w_d+Kw_p)$。相对湿度先经声明的吸湿关系转成水质量比例，RH 不能直接作为水质量分数。VFT 松弛模型 $\tau=\tau_0\exp[B/(T-T_0)]$ 只在 $T>T_0$ 计算。**32 个储存条件有 12 个位于适用数学域之外**；即使域内外推出极长时间，也不能当作货架期。

未溶药量 $U$、溶解药量 $L$、新晶体 $X$ 与吸收汇 $A$ 满足

$$\dot U=-J_d,\quad\dot L=J_d-J_p+J_r-J_a,\quad\dot X=J_p-J_r,\quad\dot A=J_a.$$

因此 $U+L+X+A$ 守恒。令 $C=L/V$，Noyes–Whitney 质量通量为 $J_d=(D/h)S_0(U/U_0)^{2/3}(C_{source}-C)_+$；若写成浓度导数，必须除以容器体积 $V$。聚合物一方面可降低药物活度和源溶解度，另一方面延迟析晶，这两种作用可能竞争。

球形晶核的自由能为 $\Delta G(r)=4\pi r^2\gamma-(4\pi/3)r^3\Delta\mu/v$。令一阶导数为零，得 $r_*=2\gamma v/\Delta\mu$、$\Delta G_*=16\pi\gamma^3v^2/(3\Delta\mu^2)$。再代入 $\Delta\mu=k_BT\ln S$：

$$\frac{\Delta G_*}{k_BT}=\frac{16\pi\gamma^3v^2}{3(k_BT)^3(\ln S)^2},\quad v=V_m/N_A,\quad S>1.$$

其中 $\gamma$ 使用 J/m²，$v$ 使用每分子 m³；输入的 cm³/mol 摩尔体积必须先乘 $10^{-6}$ 转为 m³/mol，再除以 $N_A$。进入指数的是这个无量纲量；原任务含额外有量纲能垒的写法没有被沿用。$S\le1$ 时关闭成核。

第五个状态 $z$ 是无量纲有效生长位点激活度。定义 $(x)_+=\max(x,0)$，固定聚合物浓度 $c_p$、初始剂量 $D_0$、能垒 $B=\Delta G_*/k_BT$，实现的闭合关系为

$$\dot z=k_{nuc}(S-1)_+e^{-B}-k_{loss}(1-S)_+z,$$
$$J_p=\frac{k_{growth}}{1+\beta c_p}\left\{z+s_0(X/D_0)^{2/3}\right\}(C-C_{crys})_+V,$$
$$J_r=K_{crys}(X/D_0)^{2/3}(C_{crys}-C)_+,\qquad J_a=k_aL.$$

$K_{crys}=DA_{crys}/h$ 为体积/时间，秒需转为小时；$k_{nuc},k_{loss},k_{growth},k_a$ 的单位为 h⁻¹，$\beta$ 为 mL/mg。代码在通量计算中使用有界指数和非负状态值，但保留原始求解器微小负偏差用于诊断。位点不是实测颗粒数，诱导期定义为析出量达到剂量 1%，不表示第一个分子晶核出现。方程对应 [ASD 驱动](projects/task06_asd/run_task6_asd_formulation_supersaturation_kinetics.py) 中的 `simulation(...).rhs`。

### 已计算的制剂权衡

假设 100 mg 药物置于 900 mL 容器，20% HPMC-AS 场景的封闭容器 AUC 为 **0.319124 mg·h/mL**，是晶体对照的 **7.997 倍**；$S>1.05$ 持续 **5.921 h**。六小时吸收汇质量比 **8.795**，不等于人体生物利用度。若 $J_a=k_aL$，则 $A=k_aV\int Cdt$，同一吸收汇系统中的吸收量和 AUC 存在数学依赖，不能互相充当独立验证。

30% 载药量对应更高封闭 AUC，**0.370222**，因此没有证明 20% 最优。27 条件敏感性网格给出吸收比 **5.961–8.803**，是人为选定情景范围，不是置信区间。仅标注 pH 6.8，而不规定胆盐、胶体及游离药物测量，不能称为经过验证的 FaSSIF 实验。

**参数来源与下一步决策。** API 是假想分子，大部分聚合物、输运和成核参数是假设。两种聚合物典型 Tg 有厂家来源，但不代表特定批次。应测 DSC/XRPD、吸湿、相行为、总溶解与游离分子浓度、颗粒演化，再联合校准溶出和渗透后连接 PBPK。产业意义是设计能被实验否证的制剂比较，识别活度—分子运动—析晶之间的权衡。

[输入](projects/task06_asd/inputs/parameters.json) · [结果](projects/task06_asd/results/summary.json) · [完整报告及原始来源](projects/task06_asd/ASD_FORMULATION_KINETICS_REPORT_ZH.md)

![任务6：两相共存、玻璃化转变与有限剂量过饱和](figures_omnibus/fig_task6_asd_supersaturation.png)

<a id="task7"></a>
## 7. 肿瘤免疫 QSP 与配对虚拟队列

### 延迟损伤与免疫反馈

系统使用天、mg、mg/L 和归一化 CTL 密度。Simeoni 类型增长函数为

$$G(w_0)=\frac{\lambda_0w_0}{[1+(\lambda_0w_0/\lambda_1)^\psi]^{1/\psi}}.$$

小质量时 $G\approx\lambda_0w_0$，大质量时 $G\to\lambda_1$；后者单位为 mg/day，不是承载量。小分子损伤进入三个不增殖中转池：

$$\dot w_0=G-k_2C_sw_0-k_{kill}Ew_0,\quad\dot w_1=k_2C_sw_0-w_1/\tau,$$
$$\dot w_2=(w_1-w_2)/\tau,\quad\dot w_3=(w_2-w_3)/\tau.$$

求和得 $\dot W=G-w_3/\tau-k_{kill}Ew_0$，说明进入损伤态并不立即减少总质量。三个指数驻留阶段给出平均总延迟 $3\tau$，本例为 **4.5 天**。[Simeoni 等](https://pubmed.ncbi.nlm.nih.gov/14871843/)。

令抗原通量代理 $J_{ag}=w_3/\tau$，抗体占据 $RO=C_m/(K_D+C_m)$：

$$\dot E=s_{basal}+\alpha\frac{J_{ag}}{K_{ag}+J_{ag}}-[\mu_E+k_{exh}PDL1(1-RO)]E.$$

取 $s_{basal}=(\mu_E+k_{exh}PDL1)E_{baseline}$，保证未经治疗时非零免疫平衡。抗原来源未包含直接免疫杀伤，这是模型简化。口服 PK 叠加 $FDk_a[e^{-k_et}-e^{-k_at}]/[V(k_a-k_e)]$；抗体采用 0.5 小时有限输注，而非瞬时推注。积分在每个给药变化点分段。

### 协同与事件终点

同一个配对虚拟个体内，$f_A=1-W_A/W_V$，Bliss 参考为 $f_A+f_B-f_Af_B$；组合抑制减去这一参考就是 Excess over Bliss。Loewe 则要求可求逆的单药剂量—效应曲线：$CI=d_A/D_A(f_{AB})+d_B/D_B(f_{AB})$。若单药在支持剂量范围内达不到组合效应，逆剂量就应缺失，不能外推一个数补齐矩阵。

清除率、倍增时间和基础 CTL 采用独立对数正态变异。算术均值 $m$、CV 为 $c$ 时，$\sigma^2=\ln(1+c^2)$、$\mu=\ln m-\sigma^2/2$。同一 50 个体接受四种反事实方案，因此是 **50 个独立的合成参数向量**，不是 200 位独立患者。

第 60 天载体对照（Vehicle）、靶向、抗 PD-1、联合的平均肿瘤质量分别为 **347.151、140.959、84.198、7.744 mg**。配对平均 Bliss 超额为 **0.06825**。25 个内部剂量格点中仅 **9 个**有受支持的 Loewe 指数，**16 个保留缺失**。

事件定义为首次模型肿瘤质量超过基线 120%，60 天未发生者行政删失。Kaplan–Meier 为 $\hat S(t)=\prod_{t_j\le t}(1-d_j/n_j)$。这是**模型质量进展，不是 RECIST 或临床 PFS**。四组事件数为 50、43、28、14，早期进展不会因晚期缩小而取消。RECIST 的病灶测量与其他规则未在这个质量模型中实现。[RECIST 工作组](https://recist.eortc.org/recist-1-1/)。

限制平均无事件时间对 $\hat S$ 积分到 60 天。联合相对靶向、相对抗体的差值分别为 **31.493、15.073 天**，均是条件性模拟结果。按个体成对自助抽样和按个体聚类的 Cox 不确定性保留了共享个体设计。合成 p 值、HR 和 SEM 仅衡量选定模拟变异，不能证明临床疗效或比例风险成立。

**参数来源与下一步决策。** 全部效力/PK 数值是假设；100 mg 肿瘤与 70 kg 剂量换算的组合不代表经过验证的跨物种转化。需要单药剂量范围、匹配的组合时序、肿瘤暴露、CTL/抗原观测和大小—质量观测模型。产业用途应是实验设计与机制假说比较，再用保留的纵向结果验证预测并传播不确定性。

[虚拟队列参数](projects/task07_qsp/virtual_patients.csv) · [结果](projects/task07_qsp/results_summary.json) · [项目完整报告](projects/task07_qsp/QSP_IMMUNO_ONCOLOGY_REPORT_ZH.md)

![任务7：肿瘤轨迹、Bliss曲面与模型质量无进展曲线](figures_omnibus/fig_task7_qsp_immuno_oncology.png)

<a id="task8"></a>
## 8. ADC 制造分布、载荷释放与组织输运

### 有限试剂的随机偶联

四个链间二硫键独立还原，$b\sim Binomial(4,p)$，其中 $p=1-e^{-k_{red}t_{red}}$，可反应容量 $m=2b$。联合状态 $f_{m,n}$ 同时记录容量和实际 DAR，$n\le m$。25 个状态保留奇数 DAR。游离连接子当量为 $L$ 时，

$$J_{m,n}=k_0(m-n)e^{-\alpha n}Lf_{m,n},\quad\dot f_{m,n}=J_{m,n-1}-J_{m,n},\quad\dot L=-\sum_{m,n}J_{m,n}.$$

边界外通量为零。求和可证明抗体比例及 $L+\sum n f_{m,n}$ 守恒；Gillespie 事件在有限分子群上实现相同化学计量。2、4、6、10 当量时平均 DAR 为 **2.000000、3.999320、5.912516、7.523958**。四当量条件下，24 次每次 1,000 抗体的随机重复给平均 **3.999500**，符合有限群体抽样，不要求与确定性极限完全相同。

HIC 场景叠加归一化 Voigt 峰，$y(t)=\sum_nf_nV(t-t_n;\sigma,\gamma)$，假设 $t_n=2+2n^{1.15}$ min。DAR0–8 全部参与。假设 DAR4 收集窗的面积纯度 **96.9147%**、回收率 **94.6290%** 仅适用于这个无噪声模型。已知峰形的 NNLS 恢复是算术检查，不是实测方法学验证。下游材料仍是全部四当量产物，没有偷偷替换为纯化 DAR4。

### 细胞与空间守恒

结合、内吞与分选把 ADC 送入溶酶体。有效处理通量 $v=V_{max}L_y/(K_m+L_y)$ 释放 $\overline{DAR}\,v$ 载荷。表面/回收受体、抗体及载荷各自有守恒式。连接子切割、自消除和抗体降解在此被合并，不把有效蛋白酶唯一指定为组织蛋白酶 B；特定 Val-Cit 系统的 CatB 缺失实验支持这一限制。[Caculitan 等](https://doi.org/10.1158/0008-5472.CAN-17-2391)。

十个状态 $(X,R,B,E,L_y,R_i,Z,P,Q,M)$ 依次为胞外 ADC、游离受体、结合 ADC、内体 ADC、溶酶体 ADC、回收受体、已降解 ADC、胞质载荷、已输出载荷和已代谢载荷。定义

$$J_b=k_{on}C_{ext}R,\ J_u=k_{off}B,\ J_i=k_{int}B,\ J_r=k_{rec}R_i,\ J_s=k_{sort}E,\ J_c=v,\ J_x=k_{perm}P,\ J_m=k_{met}P.$$

分子数以每细胞计，$C_{ext}[\mathrm{nM}]=X/(N_A10^{-9}V_{bath}[\mathrm L])$，$k_{on}$ 为 nM⁻¹h⁻¹。完整细胞 ODE 为

$$\dot X=-J_b+J_u,\quad\dot R=-J_b+J_u+J_r,\quad\dot B=J_b-J_u-J_i,$$
$$\dot E=J_i-J_s,\quad\dot L_y=J_s-J_c,\quad\dot R_i=J_i-J_r,\quad\dot Z=J_c,$$
$$\dot P=\overline{DAR}J_c-J_x-J_m,\quad\dot Q=J_x,\quad\dot M=J_m.$$

对应 [ADC 驱动](projects/task08_adc/run_task8_adc_dar_cleavage_bystander_dynamics.py) 的 `cell_rhs`。直接求和得到 $X+B+E+L_y+Z=X_0$、$R+B+R_i=R_0$、$\overline{DAR}(X+B+E+L_y)+P+Q+M=\overline{DAR}X_0$。受体在内吞时已进入独立回收通路，所以抗体降解释放载荷不会再次消耗受体。

球对称组织的胞外浓度满足扩散算子 $r^{-2}\partial_r(r^2D\partial_rC_e)$ 加释放、清除和交换项。实现中积分的是**球壳药量**，相邻面的通量为

$$J_{i+1/2}=D\epsilon\frac{4\pi r_{i+1/2}^2}{\Delta r}(C_{e,i}-C_{e,i+1}).$$

同一面两侧增减相消。邻细胞交换为 $k_{perm}(V_iC_{e,i}-N_i)$，胞外项取相反符号。中心零通量对应对称条件，外边界吸收作为明确统计的汇。隐式 Euler 求解 $(I-\Delta tA)N_{k+1}=N_k+\Delta tb_k$，在 80 个球壳上运行 72 小时。

半径 20 µm 的源细胞核心位于半径 200 µm 的球体内。条件性存活观测为 $S=\exp[-\int k_{max}C_i/(EC50+C_i)\,dt]$，死亡并不反过来改变几何或释放。抗原阳性源模型单向供应载荷，没有包含组织载荷向这些源细胞的反扩散。

假设高/低通透率为 0.2/0.002 h⁻¹，抗原阴性平均存活为 **0.958800/0.999848**。高通透场景的胞外 50 nM 等值线达到 **23.75 µm** 球壳中心，这不是经过验证的杀伤半径。结果显示该几何下平均旁观者效应有限，并非广泛清除肿瘤。约 $10^{-14}$ 的质量误差说明账本精度，不能衡量生物预测准确度。

**参数来源与下一步决策。** 动力学、HIC 峰形、受体丰度、通透及反应参数是假设，仅有机制文献动机。应测分 DAR 化学与检测响应、LC-MS 释放物、血浆稳定性、内吞、扩散/外排及空间共培养存活。产品 DAR 优化还需要暴露与可制造性；历史 DAR4 实例或本模型都不能证明普遍最优 DAR。

[参数来源](projects/task08_adc/parameter_provenance.csv) · [结果](projects/task08_adc/summary.json) · [完整报告及来源](projects/task08_adc/ADC_TRANSLATIONAL_ENGINEERING_REPORT_ZH.md)

![任务8：DAR、假设HIC、溶酶体供给与径向载荷暴露](figures_omnibus/fig_task8_adc_multiscale.png)

<a id="task9"></a>
## 9. 合成结构异质性、空腔与力学响应

### 公开结构不能赋予合成轨迹真实时间

结构锚点为 **T4 溶菌酶 L99A 的 X 射线结构 4W51 与 4W59**，不是癌症激酶或冷冻电镜数据。匹配后有 1,278 个重原子、164 个 Cα。核心 Kabsch 对齐后，$x(q)=(1-q)x_c+qx_b+\delta x(q)$ 加入声明的环区扰动。500 个样本形成可重复的几何基准，不是 MD、能量极小路径或 apo 平衡布居。[Merski 等](https://doi.org/10.1073/pnas.1500806112)。

核 PCA 对对齐坐标的 RBF 相似矩阵中心化并对角化，得到两条几何潜变量轴。原子序数加权点经 Gaussian 平滑产生密度代理，体素 1 Å、FWHM 2.5 Å。没有粒子图像、CTF、取向推断、FSC 或 cryoDRGN 训练；模拟平滑宽度不是实验分辨率。

### 可逆状态与布居自由能

相邻转移计数 $C$ 被转为指定邻接上的对称通量 $F=(C+C^T)/2+0.5S$，随后

$$T_{ij}=F_{ij}/\sum_jF_{ij},\qquad\pi_i=\frac{\sum_jF_{ij}}{\sum_{ij}F_{ij}}.$$

$F$ 对称直接保证 $\pi_iT_{ij}=\pi_jT_{ji}$。这是构造性可逆估计器，不是极大似然拟合。状态布居差为

$$G_i-G_0=-RT\ln(\pi_i/\pi_0).$$

该式不包含过渡态势垒。把开放态设为吸收态，首步关系 $m_i=\tau+\sum_jT_{ij}m_j$ 给出 $(I-T_{NN})m=\tau\mathbf1$。赋予的 $\tau=10$ ns 是**合成时钟**，无序实验快照不能辨识这个时间尺度。

五个状态计数为 **[225,129,79,61,6]**。估计开放—闭合布居自由能 **2.057351 kcal/mol**，条件性参数自助区间 **1.101320–3.334307**；生成器设定的差值为 **1.6 kcal/mol**。MFPT **2630.064 ns** 依赖指定时钟和状态模型。200 条自助轨迹有 22 条缺少某状态，揭示采样不足；区间不包含真实几何与物理动力学的不确定性。

### 空腔截断与奇异 Hessian

局部探针/遮挡算法统计连通体素体积，不是经过验证的全局可结合性检测器。原始无界分数保留并使用无量纲体积对数，另输出 sigmoid 作为未经校准的有界启发式量。记录的局部分量为 **76.78–209.67 Å³**，但 **500 帧中 135 帧接触 ROI 边界**，完整体积和联合门槛状态应为未知；其余 365 帧未达到指定联合门槛。扩大 ROI 还会改变连通性，数值变大不自动表示同一空腔已收敛。

弹性网络每条接触的能量为 $\tfrac12k[(\delta r_i-\delta r_j)\cdot u_{ij}]^2$。二次求导得到对角 $kuu^T$ 与非对角 $-kuu^T$ Hessian 块。六个刚体零模使 $H$ 奇异，去掉零模后使用 $\delta r=H^+F$，不能直接求普通逆。对残基 $j$ 的各向同性单位力，

$$M_{ij}=\frac13\|H^+_{ij}\|_F^2,$$

来自 $\mathbb E(ff^T)=I/3$。随机力估计随样本增加趋于这个解析期望。归一化耦合支持接触图上的路线假说，不能证明能量流或因果传播。所选 **106→11 路线仅 12.624658 Å**，没有获得原指定的 >30 Å 结果。

**参数来源与下一步决策。** 坐标来自实验，状态能量、时钟、插值、评分和弹簧模型均为声明的假设。产业接口是检查结构支持程度和几何失败原因的可审计筛选。发现新口袋前，需要目标特异的结构系综、动力学测量、侧链/溶剂精修、可结合性基准与扰动实验。[ANM](https://doi.org/10.1016/S0006-3495(01)76033-X) 和 [PRS](https://doi.org/10.1371/journal.pcbi.1000544) 提供方法依据，不能验证这条具体路线。

[结果](projects/task09_cryoem_allostery/results/summary.json) · [ROI 敏感性](projects/task09_cryoem_allostery/results/roi_sensitivity.json) · [项目完整报告](projects/task09_cryoem_allostery/CRYOEM_CRYPTIC_POCKET_REPORT_ZH.md)

![任务9：合成潜变量、状态布居自由能与机械PRS响应](figures_omnibus/fig_task9_cryoem_allostery.png)

<a id="task10"></a>
## 10. RNA 系综、配体热力学与 U1 招募

### 配分函数及实验坐标几何

ViennaRNA 在声明的最近邻模型下计算 $Z=\sum_se^{-G_s/RT}$。系综自由能为 $-RT\ln Z$，$P_{ij}=Z^{-1}\sum_{s:(i,j)\in s}e^{-G_s/RT}$。每个核苷酸还包含未配对态 $u_i=1-\sum_jP_{ij}$，因此

$$H_i=-\sum_jP_{ij}\log_2P_{ij}-u_i\log_2u_i,$$

其中 $0\log0=0$。不含未配对项的量不是完整位置状态分布的熵。26 nt 序列 `AUACUUACCUGUUCGGGAGUAAGUCU` 是**人为连接的发夹**：将假尿苷替换为 U，并用人工 UUCG 连接两条链，不是天然 SMN2 前体 mRNA。

真实执行的 14 次折叠覆盖 WT/A18C 和七个温度。37°C WT 系综能量为 **−4.840479 kcal/mol**，位置熵最大 **0.925397 bits**，该条件下没有位置超过 1.2 bits。这个阴性结果不能推广到全部温度：60°C WT 有 7 个位置超过门槛。未配对概率或位置熵均不能直接给出三级结构碱基翻转能垒。

独立几何分析使用 **6HMI/6HMO** 各 20 个 NMR 模型，保留已沉积修饰。该研究的配体是 **SMN-C5，并非利司扑兰**。NMR 模型差异不是时间序列，也不是等权 Boltzmann 布居。[Campagne 等](https://doi.org/10.1038/s41589-019-0384-5)。

560 行几何结果区分跨链磷酸中心弦长、局部空球直径、水探针净空、投影深度及粗粒度屏蔽电势。它们是操作性描述符，不是标准沟槽宽度或结合自由能。Debye–Hückel 采用均匀介电与单价盐近似，计算屏蔽长度 **8.011 Å**；电势超过约 26.7 mV 热电压时，线性响应的定量解释受到限制。

### 精确结合循环

令 $q=e^{-\Delta G_{conf}/RT}$ 表示未结合时可结合态相对于闭合态的权重，$K_O=C^\circ e^{\Delta G_{bind,O}/RT}$，$K_C$ 为闭合态解离常数。将结合与游离权重求和，得到

$$K_{D,eff}=\frac{1+q}{1/K_C+q/K_O}.$$

只有闭合态结合可忽略时，才约化为 $K_O(1+e^{\Delta G_{conf}/RT})$。结合态构象差满足 $\Delta G_{conf,bound}=\Delta G_{conf}+\Delta G_{bind,O}-\Delta G_{bind,C}$，两条路径因而热力学闭合。单凭平衡循环不能在动力学上区分诱导契合与构象选择。

四种骨架能量分解是**明确假说**，不是 docking、FEP 或真实药物亲和力。平面骨架的打开代价 5.5 kcal/mol、可结合态结合能 −12 kcal/mol 给出完整有效 $K_D$ **26.2809 µM**。配体特异的打开代价对应不同可结合亚态，而不是同一个普适转变的四个矛盾自由能。

### 有限配体与 U1 总量

每个目标/脱靶 RNA 池包含 $C,O,CL,OL,CU,OU,CLU,OLU$ 八态。令 $u=U/K_U$、$\alpha=e^{-\Delta\Delta G_{splice}/RT}$，权重为

$$[1,q,L/K_C,qL/K_O,u,qu,Lu/K_C,qLu\alpha/K_O].$$

归一化后得到结合概率，两种总量分别要求

$$L_{tot}=L+\sum_rR_{r,tot}P_r(L\text{ bound}),\qquad U_{tot}=U+\sum_rR_{r,tot}P_r(U\text{ bound}).$$

数值求根确定游离 $L,U$，无需预先强加 Hill 曲线。默认目标 RNA 2 nM、合并脱靶 RNA 30 nM、U1 100 nM。目标 U1 占据率乘 100 定义为**包含率代理**，不是实测外显子包含率。真实剪接包含此平衡模型未表达的动力学及 ATP 依赖步骤。branaplam 的 U1-C 研究不能普遍套用于利司扑兰，更不能混淆 U1-C 与 U1-A。[White 等](https://www.nature.com/articles/s41467-024-53124-5)。

在 **976** 个主要条件与 **120** 个敏感性条件中，平面 WT 代理由 **4.72469%** 上升到解析饱和极限 **93.27945%**，模型总配体 EC50 为 **1.78816 µM**。只对有支持的单调增强反应给出 EC50；非单调情况保留交点及方向，不强行压成单一 EC50。目标/脱靶合格带可以不连续或右删失；平面场景上限被剂量网格截断，不能解释为安全剂量上限。

**参数来源与下一步决策。** 折叠来自模型热力学，NMR 坐标来自公开实验，配体能量与招募耦合来自未经校准假设。A18C 二级结构对照和结合模型中的反事实“突变”是独立构造。应优先获得相同盐条件下的结合、SHAPE/DMS/NMR 约束、剪接报告浓度反应、蛋白依赖和全转录组选择性。在连接 PBPK 或药物化学排序前，还需真实化合物身份、渗透与暴露。

[结果](projects/task10_rna_splicing/results/results_summary.json) · [热力学假设](projects/task10_rna_splicing/results/thermodynamic_cycles.csv) · [项目完整报告](projects/task10_rna_splicing/RNA_TARGETED_CADD_REPORT_ZH.md)

![任务10：二级结构配对、假设热力学分量与U1占据反应](figures_omnibus/fig_task10_rna_targeted_cadd.png)

<a id="translation"></a>
## 转化路线、本地计算策略与验收门槛

### 目前本地计算能增加什么证据

本套件使用小型代数系统、ODE、有限体积输运、半经验描述符与中等规模结构矩阵，适合有资源约束的本地 CPU 运行。多个 Agent 可以独立审查方程、补充参数扫描、制图和核对来源，数值进程仍应限制线程与峰值内存。更多 Agent 不会产生更多真实实验重复、患者或独立生物验证。

下一步最有价值的本地工作应围绕决策：比较不同模型结构，判断哪些参数会改变实验选择，并保留负结果。对响应 $y(\theta)$，无量纲局部敏感性 $S_j=(\theta_j/y)\partial y/\partial\theta_j$ 只在导数和尺度有意义时适用。人为扫描不是概率分布；置信区间需要明确抽样/推断模型；结构不确定性需要比较不同模型，不能仅在同一模型内改变参数。

实验拟合应先定义观测模型 $y_k=h[x(t_k;\theta)]+\epsilon_k$，包括真实噪声和删失机制。曲线拟合良好仍可能存在参数不可辨识。应结合剖面似然、参数相关性、保留实验条件与预测检查，再判断是否识别了唯一机制。不能用同一条合成曲线同时完成“校准”和“生物验证”。

### 跨模块传递需要明确契约

| 拟连接方向 | 最低传递要求 | 当前状态 |
|---|---|---|
| 1 → 4 | 同一化合物/微观态、实测结合清除、单位换算 | 能导出描述符，尚非校准暴露链 |
| 6 → 4 | 游离溶解浓度、渗透、转运时间与首过模型 | 吸收汇尚不是已验证吸收输入 |
| 4 → 5 | 共同质量账本、受抑内在清除及匹配损失项 | Task 5 已实现 |
| 2/3 → 7 | 细胞游离暴露、靶点作用、周转与效应映射 | 未建立校准传递 |
| 8 → 7 | 肿瘤载荷暴露、空间细胞响应和兼容肿瘤观测 | 未建立校准传递 |
| 9 → 2/3 | 实验支持的结合几何与构象布居 | 目前仅合成结构假说 |
| 10 → 4/7 | 明确化合物、真实剪接效应、暴露与疾病映射 | 未建立校准传递 |

统一化合物标识和参数来源是前提。把 MPO 当清除率、GSH 比值当毒性、MSM 布居差当结合能垒，或把 U1 占据当临床外显子包含率，会制造没有证据的精确感。

### 面向产业用途的分阶段推进

1. **冻结身份与观测量。** 明确化合物、批次、介质、温度、游离/总浓度、时间单位和实际实验终点；保存原始数据及转换。
2. **校准最小可辨识模型。** 选择有信息量的浓度/时间与阴性对照；区分未知参数和实测输入，记录不确定性及缺失。
3. **检验前瞻预测。** 保留与用途相关的浓度、时序、化合物或批次，与更简单基线比较，并记录失败预测。
4. **连接有证据的接口。** 在耦合边界验证质量与单位守恒并传播不确定性，说明连接后是否改变可实验检验的决策。
5. **确认一个窄而明确的用途。** 制剂初筛、采样设计、分析方法原型或机制假说排序，比直接声称端到端临床给药预测更符合当前证据。

### 复现与证据保留

在仓库根目录使用统一入口，将本次结果写入新的工作目录：

```bash
python run_ai4pharm_omnibus_suite.py --out work/omnibus
python -m unittest discover -s tests -v
python tools/validate_repository_layout.py --out work/omnibus_layout_validation.json
```

运行模式与资源控制以根 README 和入口 `--help` 为准。各项目现有 manifest 保留其运行输入、版本和哈希。新结果应写入新的 `work/` 目录，不能为迎合新依赖版本悄悄覆盖归档数据。本次 omnibus 图是可检查数值的新多面板表达；300 DPI 是制图输出规格，不能证明已经达到可发表科学水平。

验证应包括独立解析极限、守恒、物理要求下的非负性、求解器/网格加密、输入结构与来源哈希，以及实际图形检查。这些工作处理软件与数值错误。外部实验、PK 和结构证据处理模型是否适合真实生物体系；二者不能互相替代。

当前套件的贡献是跨十类开发问题透明、可复现地比较假设及其后果。论文强度取决于明确的新问题、适当真实数据、合理基线与前瞻验证。归档中的缺失几何、VFT 无效域、Loewe 逆解不支持、ROI 截断、较短力学路径，以及 **37°C WT** 未出现高熵 RNA 位点，都是需要保留的结果，不应被隐藏或改成预设成功。
