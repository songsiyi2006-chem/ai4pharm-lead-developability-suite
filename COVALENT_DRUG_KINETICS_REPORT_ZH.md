# 任务3：共价抑制剂动力学、靶点驻留与谷胱甘肽反应性

## 研究范围与证据等级

本工作是可复现的计算演示，并非经过实验校准的候选药物预测。八种假想分子共享
芳基氨基嘧啶核心；未提供明确靶蛋白、结合构象或实测动力学数据。这些结构不是
索托拉西布、奥希替尼或伊布替尼。轨道描述符来自实际低层级量子计算，默认结合
速率与化学反应速率则是演示参数。能垒使用未经校准的描述符替代模型，没有执行
过渡态搜索。本结果不确立临床疗效、治疗窗、安全性或药物性肝损伤概率。

## 3A：亲电性与硫负离子反应能垒情景

RDKit 检查结构，以原子映射定位反应碳，生成 ETKDGv3 构象，并从收敛的 MMFF
构象中选择最低能者。默认 EHT/YAeHMOP 扩展休克尔法是真实半经验轨道计算，
不是 DFT；绝对轨道能量不可与 xTB 数值混用。可选 GFN2-xTB/ALPB 水溶剂模式
先优化中性分子，再在相同几何上计算中性和阴离子单点，保存原始日志与 JSON。

采用 I≈−EHOMO、A≈−ELUMO 的轨道近似：

$$\mu=(E_H+E_L)/2,\quad\eta=(E_L-E_H)/2,\quad
\omega=\mu^2/(2\eta),\quad S=1/(2\eta).$$

这里明确使用上述硬度与软度约定。EHT 的局部 f+ 是原子分辨的 LUMO 布居，
属于冻结轨道近似；xTB 的 f+=q(N)−q(N+1) 来自垂直有限差分。局部软度为 S f+。
程序比较迈克尔受体 β 碳与 α 碳的局部值，并导出反应碳在全部重原子中的排名。
未支持 β 位局域化时如实报告；负布居不截断、不修饰。氯乙酰胺在连氯碳发生
SN2 取代，不适用迈克尔 β 碳描述。中性微观状态、单一选定构象及阴离子束缚
状态的不确定性限制了解释；二甲氨基的质子化需要另行建模。

能垒情景为：

$$\Delta G^\ddagger=B_{warhead}+\mathrm{clip}[2(E_L-\overline{E_L}),-4,4]
\quad\mathrm{kJ/mol}.$$

系数为 2 kJ mol⁻¹ eV⁻¹，各亲电反应基团的基线在参数 JSON 中公开，未由实验拟合。
采用 1 M 标准态、传递系数为 1 的 Eyring 换算：

$$k_{thiolate}=(k_BT/h)/C^\circ\exp[-\Delta G^\ddagger/(RT)],$$
$$f_{thiolate}=1/(1+10^{pK_a-pH}),\qquad k_{GSH}=f_{thiolate}k_{thiolate}.$$

显式输入 kGSH 可覆盖该模型，此时能垒列仍是独立演示情景，并非从实测速率推断。
真正的加成或 SN2 能垒验证需统一溶剂与微观状态、过渡态优化、单一恰当虚频、
IRC 路径核验、自由能修正及实验校准；本脚本不会将替代值冒充这些计算。

## 3B：两步共价失活与非线性回归

记游离靶点为 F、非共价复合物为 B、共价复合物为 C：

$$\dot F=-k_{on}IF+k_{off}B,$$
$$\dot B=k_{on}IF-(k_{off}+k_{inact})B+k_{rev}C,$$
$$\dot C=k_{inact}B-k_{rev}C.$$

不可逆极限 krev=0，且结合达到快速预平衡时，B/(F+B)=I/(KD+I)，其中
KD=koff/kon。于是未共价修饰靶点 U=F+B 满足 Udot=−kinact I U/(KD+I)，
得到题目中的饱和双曲线。该近似要求游离药物恒定、耗竭可忽略、结合平衡快于
化学失活。KD 一般不等于回归所得 KI；准稳态 KI,QSSA=(koff+kinact)/kon。

恒定游离浓度下的精确慢衰减特征值为：

$$a=k_{on}I,\quad S_r=a+k_{off}+k_{inact},\quad
\lambda_{slow}=\frac{2ak_{inact}}{S_r+\sqrt{S_r^2-4ak_{inact}}}.$$

低浓度斜率为 kinact/KI,QSSA，高浓度趋向 kinact。完整瞬态为双指数，慢特征值
通常不是严格双曲线。程序在 1 nM–100 µM 范围积分全质量作用方程，跳过 12 个
快模态时间常数，再以非线性最小二乘拟合 A exp(−kobs t)。第二阶段拟合正参数
双曲线，使用对数速率残差，使各浓度近似按相对误差等权。表格分别列出 KD、
KI,QSSA、拟合 KI、拟合 kinact、各效率定义及双曲线近似误差。

这些轨迹没有加入实验噪声，因而不将回归协方差冒充实验置信区间。模拟采样时间
随慢速率调整，可能远超实际实验时长；真实实验需要限定窗口、重复测量及参数
可辨识性检验。图1统一使用 krev=0 比较正向效率；细胞洗脱和动态浓度模拟恢复
各化合物的实际设定逆反应。可逆共价分子的正向效率不能替代平衡占有率或驻留。
瞬时抑制为 B+C，未共价修饰比例则为 F+B，两者不是同一实验读数。

另有有限初始剂量模型，同时包括药物耗竭、清除、GSH 捕获和逆向释放，并核验
I+B+C+SG+已清除药物的质量守恒。该六小时模型不含蛋白周转；周转由独立洗脱
模块处理。本任务未包含底物竞争。题设 <10³、10³–10⁴、10⁴–10⁶、>10⁶ M⁻¹s⁻¹
分级仅作筛选标签，不能将其中某一范围称为普遍适用的临床最优标准。

## 3C：靶点驻留、蛋白周转与快速洗脱

亲和力参与识别与浓度依赖的占有过程，驻留时间也不能独自决定临床疗效；游离
暴露、靶点合成、通路对占有率的敏感性、组织分布和选择性同样重要。共价占有
可在游离药物清除后持续，Bradshaw 等的可逆共价激酶抑制剂实验提供了实例。

全部靶点状态以 kdeg=ln2/t1/2,protein 降解，新蛋白以
ksyn=kdeg Etotal,baseline 合成到游离池：

$$\dot F=k_{syn}-k_{deg}F-k_{on}IF+k_{off}B,$$
$$\dot B=k_{on}IF-(k_{off}+k_{inact}+k_{deg})B+k_{rev}C,$$
$$\dot C=k_{inact}B-(k_{rev}+k_{deg})C.$$

先恒定游离药物处理两小时，随后瞬时将游离药物降为零。已结合药物保留；其后
释放的药物立即移除，不考虑细胞药物库与外源再结合。不过 C→B 后仍可在同一
复合物内部再次形成共价键。脚本也单独积分题设简化模型
Fdot=ksyn−(kdeg+kobs)F，洗脱后解析解为
F(t)=Etotal−[Etotal−F(0)]exp(−kdeg t)，并核验数值解。

单态非共价抑制剂本征驻留时间为 1/koff，含周转时的细胞结合寿命为
1/(koff+kdeg)。理想不可逆键的化学寿命无限；1/kdeg 是蛋白/加合物的平均寿命，
ln2/kdeg 才是周转控制的恢复半衰期。可逆共价分子的 1/krev 仅为键打开的平均
等待时间，并非完整解离驻留时间。程序通过结合态矩阵 Q 计算洗脱时初始占有
加权的平均首次离开时间 tau=1ᵀ(−Q)⁻¹p_bound(0)，包含解离、断键、再次成键和
蛋白降解。恢复半时指抑制降至洗脱初值的一半，未在窗口内达到则记为 null。

## 3D：GSH 共轭孪生模型与反应性权衡

缓冲恒定 GSH 默认 5 mM。不可逆拟一级动力学为
I(t)=I0 exp(−kGSH[GSH]t)，t1/2=ln2/(kGSH[GSH])。数值积分与解析结果互相核验。
可逆 GSH 加合物的比例为 p[1−exp(−(a+b)t)]，其中
a=kGSH[GSH]、b=kreverse,GSH、p=a/(a+b)。输出正向半衰期、弛豫半衰期、
平衡加合物比例及达到总量 50% 加合的时间；p≤0.5 时最后一个指标不存在。

题设安全比 (kinact/KD)/kGSH 无量纲，同时给出 QSSA 版本。该比值没有包含暴露、
蛋白组选择性或临床结局，因此并非真实治疗指数。图3以靶向效率和 GSH 正向
半衰期同时越大越好绘制 Pareto 前沿，阴影标记题设筛选区。虽然文件名保留 radar，
内容按要求采用散点图。正向半衰期低于 30 分钟触发高反应性筛选标记，而非经过
验证的特异质性 DILI 预测。GSH 共轭可能有解毒作用；耗竭与蛋白共价修饰取决于
暴露、再生、代谢及细胞区室。恒定 GSH 试验本身并未模拟细胞 GSH 耗竭。
Flanagan 等提供 GSH 本征反应性测量方法；这类测量并不直接给出临床安全结论。

## 可复现性、使用及限制

随机种子仅控制构象生成，动力学轨迹未加随机噪声。软件版本、文件 SHA-256
及验证结果均保存；检查包括浓度非负性、质量守恒、解析/数值一致性、速率极限、
蛋白恢复和图片 300 DPI。同一软件构建下应可复现，但不承诺跨平台逐字节一致。
图4是混合机理、独立演示速率与计算 LUMO 的探索性线性回归，无论相关性强弱均
不能证明因果或预测能力。xTB 与 EHT 均为近似方法，默认能垒和速率未作实验验证。

运行 `python run_task3_covalent_kinetics_residence_time.py --help` 查看参数。
可修改 `data_task3/parameters_template.json` 后通过 `--parameters` 导入带来源的
速率；仅覆盖部分字段时，其余仍为演示值，来源说明应写清这一点。kGSH 覆盖值
按当前实验条件的表观速率使用，不自动作 pH 修正。`--self-test-only` 执行核心
数学检查。`--git-sync` 要求已有 main 分支检出及 origin 远端，仅暂存生成文件，
正常提交推送，不强制推送或切换分支。重复运行幂等更新 README 区块并覆盖专用
输出路径；请勿同时向同一目录运行多个实例。

## 本次结果

- 电子结构方法：EHT/YAeHMOP; frozen-orbital LUMO population proxy。
- 温度 310.15 K；pH 7.4；GSH pKa 假设 8.7；GSH 5 mM。
- Pareto 候选：W02, W06, W08。
- 高反应性标记：W03, W04, W05。
- β 位局域化未获支持或贡献很弱的分子：W05, W06。
- 探索性 LFER R²：0.2992。
- ODE 拟合 kobs 相对特征值的最大误差：7.81e-08。

| ID | Warhead | LUMO (eV) | kinact/KD (M^-1 s^-1) | GSH forward t1/2 (min) | Ratio | Flag <30 min |
|---|---|---:|---:|---:|---:|---|
| W01 | Acrylamide | -9.479 | 3e+04 | 68.73 | 8.92e+05 | no |
| W02 | Dimethylaminomethyl acrylamide | -9.364 | 2e+04 | 240.32 | 2.08e+06 | no |
| W03 | Alpha-cyanoacrylamide | -9.780 | 6e+04 | 3.60 | 9.36e+04 | yes |
| W04 | Chloroacetamide | -9.329 | 4.5e+04 | 3.47 | 6.76e+04 | yes |
| W05 | Vinyl sulfone | -9.330 | 1.07e+04 | 11.09 | 5.12e+04 | yes |
| W06 | Methacrylamide | -9.329 | 400 | 790.61 | 1.37e+05 | no |
| W07 | Beta-methyl cyanoacrylamide | -9.583 | 9e+04 | 43.02 | 1.68e+06 | no |
| W08 | Beta-isopropyl cyanoacrylamide | -9.509 | 1.2e+05 | 145.82 | 7.57e+06 | no |

### 各化合物参数来源

- W01: kinetics = Illustrative scenario; no experimental calibration; GSH = Eyring/LFER scenario; uncalibrated
- W02: kinetics = Illustrative scenario; no experimental calibration; GSH = Eyring/LFER scenario; uncalibrated
- W03: kinetics = Illustrative scenario; no experimental calibration; GSH = Eyring/LFER scenario; uncalibrated
- W04: kinetics = Illustrative scenario; no experimental calibration; GSH = Eyring/LFER scenario; uncalibrated
- W05: kinetics = Illustrative scenario; no experimental calibration; GSH = Eyring/LFER scenario; uncalibrated
- W06: kinetics = Illustrative scenario; no experimental calibration; GSH = Eyring/LFER scenario; uncalibrated
- W07: kinetics = Illustrative scenario; no experimental calibration; GSH = Eyring/LFER scenario; uncalibrated
- W08: kinetics = Illustrative scenario; no experimental calibration; GSH = Eyring/LFER scenario; uncalibrated

## 参考资料

- [Flanagan et al. (2014), experimental GSH reactivity and computational methods](https://doi.org/10.1021/jm501412a)
- [Bradshaw et al. (2015), reversible covalent inhibitors and tunable residence](https://doi.org/10.1038/nchembio.1817)
- [RDKit: rdEHTTools / YAeHMOP API](https://www.rdkit.org/docs/source/rdkit.Chem.rdEHTTools.html)
- [xTB: orbital properties and machine-readable JSON output](https://xtb-docs.readthedocs.io/en/latest/properties.html)
- [Assay Guidance Manual: mechanism-of-action assays](https://www.ncbi.nlm.nih.gov/books/NBK92001/)
