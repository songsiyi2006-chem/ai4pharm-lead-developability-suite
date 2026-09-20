# Task 5–7：计算实现与交付记录

[返回七任务首页](../README.md) · [原始任务与公式修正](SCIENTIFIC_NOTES.md)

计算执行日期：2026-09-20，对应提交 `3785df8`。本页导航已更新为新的项目目录；历史验收记录中的旧路径和哈希保留原始含义。当前迁移检查见 [目录整改记录](reorganization/README.md)。

三个任务均提供可运行 Python 脚本、参数记录、计算数据、四张 300 DPI 图和中英文报告。此次没有新增湿实验、患者数据或临床模型标定。默认计算是明确标注的机制模拟；文献方法和少量已核实的材料信息不会自动使其成为药物或制剂的真实预测。

| 项目 | 实际完成的计算 | 项目入口 |
|---|---|---|
| Task 5：CYP DDI/TDI | 6 个药物名称对应的假设参数情景、3 个 CYP 亚型、14 天 BID 给药及停药后 10 天；2 个探针在首日和第 14 天给药，共 24 组动态/静态比较；酶恢复、周转敏感性和较严求解容差复核 | [代码与结果](../projects/task05_cyp_ddi/) · [中文](../projects/task05_cyp_ddi/CYP_DDI_KINETICS_REPORT_ZH.md) · [EN](../projects/task05_cyp_ddi/CYP_DDI_KINETICS_REPORT_EN.md) |
| Task 6：ASD | 4 个聚合物情景、324 个热力学网格点、32 个储存条件/载药量组合、69 次溶出/吸收及敏感性 ODE 运行；求解共同切线与旋节线 | [代码与结果](../projects/task06_asd/) · [中文](../projects/task06_asd/ASD_FORMULATION_KINETICS_REPORT_ZH.md) · [EN](../projects/task06_asd/ASD_FORMULATION_KINETICS_REPORT_EN.md) |
| Task 7：肿瘤免疫 QSP | 同一组 50 名虚拟患者的四方案配对模拟，60 天、48,200 行轨迹；6×6 剂量矩阵与用于单药求逆的 78 个剂量条件；Bliss、Loewe、KM、配对置换和描述性 Cox 分析 | [代码与结果](../projects/task07_qsp/) · [中文](../projects/task07_qsp/QSP_IMMUNO_ONCOLOGY_REPORT_ZH.md) · [EN](../projects/task07_qsp/QSP_IMMUNO_ONCOLOGY_REPORT_EN.md) |

## 需要怎样理解结果

Task 5 实际导入并调用仓库原有 Task 4 的 PBPK 实现，再加入酶耗竭/再合成和探针清除的变化。基础筛查比值、官方静态简化式、包含肾清除的模型极限和动态 AUC 分别保存。默认药物名称用于演示不同作用模式，数值参数未拟合真实临床 PK。恢复时间按最后一次给药后酶丰度持续回到基线 90% 定义，超过观察窗的结果保留为删失；这不等于可逆抑制已经消失。独立文献附录还保存 Quinney 2010 的三个已发表模型参数及 246 行响应计算，明确区分文献模型估计和新测定值。

Task 6 区分质量分数与体积分数，保存成核指数的量纲和有限药物质量守恒。两个聚合物典型 Tg 有厂商来源，其他未测参数明确属于假设。32 个储存情景中有 12 个超出所设 VFT 公式的适用域，未强行计算货架寿命。20% HPMC-AS 是预设候选配方；示例中溶出 AUC 约为晶型的 8.00 倍，吸收量代理比约为 8.80，均不是人体生物利用度增益，也没有证明工业最优配方。

Task 7 用相同随机个体在四个方案中作配对比较。模型中组合方案表现出正的 Bliss 超额效应，但 25 个内部组合点中仅 9 个能在单药支持范围内估计 Loewe 指数，其余 16 个保留缺失及原因。模拟质量进展不是 RECIST 1.1 或临床 PFS；置信区间仅反映指定虚拟人群和模型下的变异，不包括未知生物机制、参数失配或临床外推误差。

## 如何复跑

从仓库根目录运行，使用已具备依赖的 Python。输出目录必须尚不存在，防止覆盖冻结结果：

```powershell
python projects/task05_cyp_ddi/run_task5_cyp_ddi_mechanism_based_inhibition.py --out work/task5_reproduce
python projects/task06_asd/run_task6_asd_formulation_supersaturation_kinetics.py --out work/task6_reproduce
python projects/task07_qsp/run_task7_qsp_tumor_immune_pkpd_synergy.py --out work/task7_reproduce
python -m unittest discover -s tests -v
```

三项脚本均提供 `--self-test`；输入模板和自定义参数的格式见各项目 README。Task 5 需要仓库中的 Task 4 模块，已记录其来源哈希。主程序不会自动提交或推送 Git，也不需要联网下载模型。

Windows 本次复用了已有 `chem-ai4s` Python 3.12 环境执行新计算；完整仓库测试使用含已有 seaborn 的 Miniconda base Python 3.14 环境。没有为此安装或替换环境。原有 Task 1–4 科学实现与结果保留，首页、导航和共享 README 哈希按新索引更新。

数值检查、独立复跑、文件哈希、图像元数据和仓库完整性证据见 [原始验收记录](archive/task5_7_delivery_2026-09-20/validation/summary.json)。这是一份软件与示例计算交付，不代表临床或工业验收通过。
