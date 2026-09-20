# Task 5 · CYP inhibition and recovery

[返回仓库首页](../../README.md) · [项目列表](../README.md)

任务五把 CYP3A4、CYP2D6、CYP2C9 的可逆抑制与酶失活/再合成，接入仓库中**实际的 [Task 4 七状态药物分布模型](../task04_pbpk/README.md)**。六个药名用于情景定位；主要 PK、Ki、KI、kinact 数值均为明确列出的假设，不能用于临床给药或停药决策。

[TASK.md](TASK.md) 保存用户原始需求，未经改写；其中的目标、公式和预期结果不等于已经验证的结论。实际实现、必要修正与完成情况以下列报告和结果为准。

| 入口 | 内容 |
|---|---|
| [中文技术报告](CYP_DDI_KINETICS_REPORT_ZH.md) / [English report](CYP_DDI_KINETICS_REPORT_EN.md) | 机制、数学修正、实际计算、证据边界与优化实验 |
| [执行脚本](run_task5_cyp_ddi_mechanism_based_inhibition.py) | 实际导入 Task 4；离线、确定性、拒绝覆盖已有输出 |
| [可编辑六药面板](compound_panel_example.json) | `--config` 输入；完整保留假设标签 |
| [结果摘要](results/results_summary.json) | 24 个探针对照、36 个酶恢复条目与数值检验 |
| [参数证据台账](results/parameter_provenance.csv) | 每项数字的假设/用户指定来源 |
| [文献参数](results/literature_parameter_reference.csv) / [响应计算](results/literature_parameter_response.csv) | Quinney 2010 三个已发表模型估计；独立于六药 PBPK |
| [动态与静态比较](results/dynamic_static_comparison.csv) | 按相同观察时窗计算 AUCR |
| [酶恢复](results/washout_recovery.csv) / [周转敏感性](results/turnover_sensitivity.csv) | 从末次给药起算，保留残余抑制剂 PK |
| [四张 300 DPI 图](results/figures_task5/) | 预孵育、双探针 PK、条件筛选图和恢复曲线 |
| [来源](results/sources.json) / [SHA256 清单](results/manifest.json) | 一手来源与可追溯结果 |

从仓库根目录运行，使用已有的 NumPy、SciPy、Matplotlib：

```sh
python projects/task05_cyp_ddi/run_task5_cyp_ddi_mechanism_based_inhibition.py --self-test
python projects/task05_cyp_ddi/run_task5_cyp_ddi_mechanism_based_inhibition.py --config projects/task05_cyp_ddi/compound_panel_example.json --out work/task5_reproduction
python projects/task05_cyp_ddi/run_task5_cyp_ddi_mechanism_based_inhibition.py --write-example work/editable_panel.json
```

输出目录必须不存在，避免覆盖冻结的计算。配置只允许 `compound_panel` 键，必须包含六个原有药名和所有字段，允许修改数值；此版本不把用户修改的参数自动升格为实测证据。运行需要仓库中的 Task 4 脚本，但不修改它；不是脱离仓库即可运行的独立单文件模型。`--self-test` 不写入结果。

四幅图的位置是 `results/figures_task5/`，保持任务内代码、报告和结果集中。FDA/ICH 筛选阈值只用于条件性模型演示，低于阈值不标注为“安全”，AUCR ≥ 5 也不自动判定禁忌。三个 CYP 是用户指定子集，不能替代完整监管面板。
