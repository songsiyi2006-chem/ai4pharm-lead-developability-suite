# Task 7 — 肿瘤免疫联合用药 QSP / Synthetic tumor–immune QSP

[返回仓库首页](../../README.md) · [项目列表](../README.md)

本模块已执行 50 名虚拟个体、4 个配对治疗臂、60 天模拟，以及 6×6 剂量矩阵和 Loewe 单药反演。它是可复现的研究原型；全部药物、肿瘤、免疫参数和个体差异均为明确标注的假设，未使用患者数据或完成临床验证。

[TASK.md](TASK.md) 保留未经改写的用户原始需求。原文中的目标与预期结果不代表已经达到的科学结论；以下报告、数据和验证记录说明实际完成情况及证据边界。

| 定位 | 文件 |
|---|---|
| 完整独立驱动脚本 | [run_task7_qsp_tumor_immune_pkpd_synergy.py](run_task7_qsp_tumor_immune_pkpd_synergy.py) |
| 双语技术报告 | [中文](QSP_IMMUNO_ONCOLOGY_REPORT_ZH.md) / [English](QSP_IMMUNO_ONCOLOGY_REPORT_EN.md) |
| 参数与文献来源 | [parameters_used.json](parameters_used.json) / [sources.json](sources.json) |
| 50 名虚拟个体 | [virtual_patients.csv](virtual_patients.csv) |
| 48,200 行状态轨迹 | [trajectories.csv](trajectories.csv) |
| 结果与配对生存统计 | [results_summary.json](results_summary.json) / [survival_statistics.json](survival_statistics.json) |
| 剂量与动态协同 | [dose_synergy_matrix.csv](dose_synergy_matrix.csv) / [bliss_timecourse.csv](bliss_timecourse.csv) |
| 单药支持域与 PK | [loewe_monotherapy_support.csv](loewe_monotherapy_support.csv) / [pk_exposure.csv](pk_exposure.csv) |
| 验证、运行日志、SHA256 | [verification.json](verification.json) / [run_log.json](run_log.json) / [manifest.json](manifest.json) |
| 四张 300 DPI 图 | [figures_task7](figures_task7/) |

从仓库根目录运行，要求 Python ≥3.10、NumPy、SciPy、Matplotlib；测试另需 Pillow。普通运行没有网络或 Git 副作用。

```bash
python projects/task07_qsp/run_task7_qsp_tumor_immune_pkpd_synergy.py --self-test
python projects/task07_qsp/run_task7_qsp_tumor_immune_pkpd_synergy.py --out work/task7_reproduction
python projects/task07_qsp/run_task7_qsp_tumor_immune_pkpd_synergy.py --write-example work/qsp_parameters.json
python projects/task07_qsp/run_task7_qsp_tumor_immune_pkpd_synergy.py --config work/qsp_parameters.json --out work/task7_custom
python -m unittest discover -s tests -p test_task7.py -v
```

`--out` 必须是新目录，以保护已存证结果。交付目录中的生成文件来自同一次完整运行；报告和本 README 是人工编写的分析说明。`manifest.json` 记录驱动脚本与生成文件哈希；根目录交付清单另覆盖报告。耗时日志与软件版本随环境变化，CSV/数值结果在相同软件环境下可重现。

自定义 `n_patients` 和 `horizon_day` 会同时更新生成数据、图题和删失终点。通用字段使用 `endpoint_day`、`eob_endpoint_mean`、`endpoint_mass_mean_mg` 等名称；静态双语报告描述本次 50 人、60 天默认计算。10 人、30 天的完整运行回归检查见 [自定义运行测试日志](validation_custom_run.txt)。

第 60 天配对个体 Bliss 均值为 **0.06825 ±0.00686 SEM**；25 个内部剂量格点中仅 **9 个**有支持域内 Loewe 指数。质量超过基线 120% 的终点是模型定义，**不是 RECIST 1.1**。模拟 HR、置换检验和均值区间仅描述假设模型及有限虚拟样本，不能解释为真实治疗获益。
