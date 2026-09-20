# 任务四：PBPK 与基于暴露的剂量筛选

## 摘要
化合物：**SYNTHETIC_LEAD_DEMO**；合成示例输入：**True**。
本交付为可复现的简化七状态 PBPK 科研模型，不是经过临床验证的预测平台。未提供指定先导化合物的实测清除率、组织分配实测值、毒理学数据或人体药代观测，因此示例不能产生可靠的人体推荐剂量。毒性阈值为演示假设，不代表已确立的安全上限。将 synthetic_inputs 改成 false 只改变标签，不构成模型验证。

筛选得到的**科研情景方案**为：每次 15 mg，BID（30 mg/日），覆盖率 90.443%，峰浓度 0.2802 mg/L。机器可读结果中的 clinical_recommendation 保持为空。QD 表示每 24 小时一次，BID 表示每 12 小时一次；剂量按活性化合物质量计，不自动进行盐型质量换算。

## 药学工程意义与范围
PBPK 将分子性质、器官容量、血流和清除机制连接起来，用于解释给药剂量如何形成动态体内暴露。把游离暴露与药效阈值比较，可用于提出制剂或给药频率的研究假设，并确定下一步应优先测量的参数。这里把 IC90 当作游离浓度阈值；若体外数值是培养体系名义浓度，必须另行校正蛋白结合和测定条件。超过 IC90 的时间比例仅是暴露替代指标，不能等同于受体占有率或临床疗效。

七个药物量状态为血液/血浆、肝、胃肠吸收库、肾、脑、肺、其余组织。中央状态存储全血药物量，再通过 Rb 换算混合静脉血浆浓度。胃肠为腔内药物库，不是灌流肠壁；肠壁、门静脉与动脉血均未单独建室。肺处于体循环串联位置；肝血流汇总门静脉及肝动脉供应。按照任务中的简化肾清除写法，肾排泄从中央室扣除，肾组织室只承担分布，避免重复清除。其余组织是有效成分混合室，默认用作靶组织代理，并非经过验证的独立靶器官。

## 参数与组织分配
完整输入与单位见 parameters_used.json，可编辑模板见 compound_example.json。MW=450 g/mol，cLogP=3，电离类型=base，酸/碱 pKa=None/7.8，fu,p=0.1，Rb=1，CLint,mic=8 μL/min/mg，fu,mic=1，ka=1 h^-1，Fa=0.9，Fg=1。MW 用于浓度单位换算，不在该分配公式中独立决定 Kp。中央血容量=5 L；器官体积、血流和组成属于固定成人情景假设，不代表个体化生理参数。肝质量 1,800 g 与分布体积 1.8 L 分开使用。

| 组织 | 体积 L | 血流 L/h | 组织:血浆 Kp | 中性脂质、磷脂、水体积分数 |
|---|---:|---:|---:|---|
| liver | 1.8 | 90 | 2.74764 | 0.0138,0.0303,0.705 |
| kidney | 0.31 | 66 | 2.38782 | 0.0121,0.024,0.783 |
| brain | 1.4 | 45 | 6.25609 | 0.0391,0.0533,0.77 |
| lung | 0.5 | 360 | 0.925614 | 0.003,0.009,0.811 |
| rest | 60.99 | 159 | 6.13669 | 0.05,0.015,0.65 |

组织分配采用 Poulin–Theil 组成公式[1]，并明确以 Henderson–Hasselbalch 估算的 logD7.4 处理电离；组织蛋白结合采用可配置的间质/血浆蛋白比例。这是等 pH、被动平衡下的近似，**不是完整 Rodgers–Rowland 或 Berezhkovskiy 实现**。未描述酸性磷脂结合、溶酶体捕获、主动转运或血脑屏障通透性。组织组成数值为明确列出的演示假设，脑和混合剩余组织尤其需要实测验证。支持按组织覆盖 Kp；两性化合物仅使用简化电离近似。

## IVIVE 与质量守恒方程
微粒体清除率按基于孵育总浓度的表观清除率解释，除以 fu,mic 后外推；若输入已经是游离清除率，应设 fu,mic=1。乘以 60/10^6 完成 μL/min 到 L/h 的转换。本情景 CLint,u=38.88 L/h，肝血清除率 CLH,b=3.72699 L/h，低浓度肝逃逸率 FH=0.958589。FH 仅是肝逃逸率，不是总体口服利用度；肝内代谢已经实现首过，吸收输入不能再次乘 FH。

以下代码块列出完整公式，M 为 mg，浓度为 mg/L，流量及清除率为 L/h，时间为 h。血液流出浓度为 Rb*C_i/Kp_i；这种换算避免将血流与血浆浓度直接相乘。

```text

Let M_B be mixed-venous blood drug mass; Cp=M_B/(V_B Rb).
For tissue i, C_i=M_i/V_i and C_vi=Rb*C_i/Kp_i (venous blood).
C_a=C_v,lung; Qc=sum_i Q_i, i=liver,kidney,brain,rest.

    dM_g/dt = -ka*M_g
    dM_lung/dt = Qc*(Rb*Cp - C_a)
    dM_i/dt = Qi*(C_a-C_vi)                      [kidney, brain, rest]
    dM_L/dt = QH*(C_a-C_vL) + ka*Fa*Fg*M_g - H
    dM_B/dt = sum_i(Qi*C_vi) - Qc*Rb*Cp - CLrenal,p*Cp
    H_linear = CLint,u * fu,p * C_L/Kp,L
    H_saturable = Vmax*Cu,L/(Km+Cu,L), Cu,L=fu,p*C_L/Kp,L
    Vmax = CLint,u*Km; CLrenal,p = GFR*fu,p

Bookkeeping: dE_H/dt=H; dE_R/dt=CLrenal,p*Cp;
dE_pre/dt=(1-Fa*Fg)*ka*M_g; dAUC/dt=Cp.
M_B+M_g+sum(M_tissues)+E_H+E_R+E_pre = cumulative administered dose.

    CLint,u [L/h] = (CLint,mic / fu,mic) * MPPGL * liver_mass_g * 60/1e6
    fu,b = fu,p/Rb
    CLH,b = QH*fu,b*CLint,u / (QH+fu,b*CLint,u)
    EH = CLH,b/QH; FH = 1-EH; F_linear = Fa*Fg*FH
    Cfree [nM] = Cp [mg/L] * fu,p * 1e6 / MW [g/mol]

Kp approximation (common pH 7.4):
    D = 10^cLogP / (1 + acid_term + base_term)
    acid_term = 10^(7.4-pKa_acid); base_term = 10^(pKa_base-7.4)
    [include only terms appropriate to the selected ionization class]
    fu,t = 1/(1 + protein_ratio*(1/fu,p-1))
    A(f) = D*(f_nl+0.3*f_ph)+f_water+0.7*f_ph
    Kp,t = A(tissue)/A(plasma) * fu,p/fu,t

Linear model dM/dt=A*M:
    AUC_0_inf = e_B^T*(-A)^(-1)*M0/(V_B*Rb)
    M_ss,post = (I-exp(A*tau))^(-1)*dose_vector

```

默认模型为线性。设置 km_unbound_mg_l 后启用 Michaelis–Menten 饱和代谢，并令 Vmax=CLint,u*Km。仅凭微粒体 CLint 不能推定 Km，必须另行提供或明确假设。此时良好搅拌模型指标只代表低浓度极限；同剂量口服/静脉 AUC 比值应解释为表观暴露比，不一定等于进入循环的绝对剂量比例，也可能超过 100%。非线性剂量筛选逐剂量计算，不使用线性叠加。

## 数值计算与结果
采用 Radau 或 BDF 分段积分，每次给药以精确状态跳跃加入胃肠药物库。0、12、…、156 小时共给药 14 次，168 小时记录给药前值。七个药物量状态之外，另有肝消除、肾消除、首过前损失和 AUC 四个记账积分，不把它们计入生理房室。

线性模型通过矩阵方法独立核验 ODE、计算完整 AUC 尾部，并直接求周期稳态。非线性 AUC 积分至残余质量低于剂量的 10^-9 后，补入低浓度尾部近似；超出迭代或时间上限会报错。半衰期为模型低浓度渐近特征值对应的半衰期，口服指标包含吸收限速的可能，不把 48 小时图形回归值伪装成终末半衰期。单次给药的 Cmax 在 0–48 小时内搜索。

| 指标 | 结果 |
|---|---|
| IV Cmax / Tmax | 20 mg/L / 0 h |
| Oral Cmax / Tmax | 0.393085 mg/L / 0.389014 h |
| IV AUC 0-inf | 22.3364 mg h/L |
| Oral AUC 0-inf | 19.2703 mg h/L |
| Oral apparent terminal half-life | 62.5233 h |
| AUC-ratio oral F | 86.273% |
| Oral AUC after 48 h | 57.946342612108964% |
| Day 7 pre-dose trough | 1.25124 mg/L |
| True SS peak / pre-dose trough | 1.86809 / 1.48126 mg/L |
| True SS within-interval minimum | 1.48126 mg/L |
| True SS 24h coverage | 100% |
| Cmax accumulation ratio | 4.752378630997845 |
| Day 7 vs SS state relative error | 0.155287 |
| Day 7 reached SS (1% state/peak/trough criterion) | False |

第七天结果与独立周期稳态比较，不预设七天足够达到稳态。同时报告给药前谷浓度和给药间隔内实际最低浓度，因为吸收延续可能使两者不同。蓄积比使用真正稳态峰浓度除以首次给药后 0–12 小时峰浓度。完整数值、质量守恒误差和周期边界残差见 results_summary.json。

## 剂量筛选、局限与复现
对 QD、BID 分别在每次 10–800 mg、步长 1 mg 的网格上求最小剂量，要求每 24 小时游离浓度超过 IC90 的时间不少于 90%，且**总血浆**峰浓度严格低于 5 mg/L。默认 QD 日剂量范围为 10–800 mg，BID 为 20–1600 mg。阈值交点使用求根精化后按持续时间计算，避免仅统计采样点。两种频率按总日剂量比较，相同日剂量时选择峰值较低者。图中的阴影根据绘图网格展示可行区域，最终候选剂量另在搜索网格中验证。

sensitivity_analysis.csv 提供 fu,p、微粒体清除率、吸收及全部 Kp 的单因素变化结果，用于显示情景敏感性；它不是置信区间、群体变异预测或全局敏感性分析。模型缺少制剂溶出、肠壁代谢实测依据、转运、器官疾病、个体差异、毒理和人体观测。需要外部实验与人体数据检验预测表现，才能支持临床剂量决策。FDA 报告指南[3]强调清楚说明模型用途及支持证据；数值测试通过不等于临床验证。

运行 `python projects/task04_pbpk/run_task4_pbpk_pharmacokinetics_dose_prediction.py --self-test --out work/task4_reproduction` 可复现交付。使用 `--write-example work/task4_compound.json` 导出模板，编辑后通过 `--config work/task4_compound.json` 读取。verification.json 保存检查结果，manifest.json 保存版本和文件校验值，四张图均为 300 DPI。脚本自包含，可重新生成所有表格和中英文报告；Task 1 描述符可通过 JSON 模板导入，但不能替代蛋白结合、清除率或毒理学实测参数。

## 参考资料
1. Poulin & Theil (2002), tissue-composition partitioning: https://doi.org/10.1002/jps.10005
2. Poulin & Theil (2002), generic PBPK models: https://doi.org/10.1002/jps.10128
3. FDA (2018), PBPK report format and content: https://www.fda.gov/regulatory-information/search-fda-guidance-documents/physiologically-based-pharmacokinetic-analyses-format-and-content-guidance-industry
