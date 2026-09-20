# 十任务整合交付 · 2026-09-21

[仓库首页](../../README.md) · [中文总报告](../../AI4PHARM_DECADE_TREATISE_ZH.md) · [English](../../AI4PHARM_DECADE_TREATISE_EN.md) · [原始要求](TASK.md)

本次完成统一运行入口、10 张多面板 300 DPI 图、两份十章学术总报告，以及支持筛选、搜索、放大和来源导航的离线 HTML 图谱。三个 Agent 分别完成前五图、后五图和双语报告，又交叉审阅数值、公式、来源与执行器；主任务负责集成、本机计算和发布。

## 两条证据链

1. **冻结基准与综合展示。** 根目录 `figures_omnibus/` 由原 `projects/` 中的归档数值重新绘制，未改写原项目。38 个绘图输入的 SHA256 记录在 [source_hashes.json](../../data_omnibus/source_hashes.json)。报告中的数字明确对应这些归档计算。
2. **本次独立复算。** 通过统一入口在本机顺序执行全部 10 个原始驱动，数值库限制单线程，无 GPU、HPC 或新增实验。10 项退出码均为 0，总子任务用时 **342.125 秒**。这只是本次机器和环境的记录，不是通用性能保证。复算之后也实际重绘了全部 10 图，检验新输出的目录适配，没有旧数值回填。

精选复算输入按任务存放在 [data_omnibus/recompute/projects/](../../data_omnibus/recompute/projects/)，共 38 个 CSV/JSON/NPZ 文件；完整本次工作目录仍保存在本地 `work/omnibus_recompute_20260921/`。精选数据包含全部综合图所需数值，不宣称包含每个原驱动的所有中间文件。

[逐任务执行记录](../../data_omnibus/recompute/execution.json) · [原始运行日志](../../data_omnibus/recompute/logs/) · [复算输入哈希](../../data_omnibus/recompute/selected_source_hashes.json) · [运行环境](../../data_omnibus/recompute/run_summary.json)

| Task | 本次子进程时间 / s | 状态 |
|---|---:|---|
| 1 · MPO / ADMET | 170.205 | 完成；构象部分收敛按原结果保留 |
| 2 · TPD | 50.549 | 完成 |
| 3 · 共价动力学 | 20.432 | 完成；EHT 路径 |
| 4 · PBPK | 8.949 | 完成 |
| 5 · CYP DDI | 32.699 | 完成 |
| 6 · ASD | 11.016 | 完成 |
| 7 · QSP | 8.443 | 完成 |
| 8 · ADC | 5.997 | 完成 |
| 9 · 合成构象 / PRS | 28.745 | 完成 |
| 10 · RNA | 5.091 | 完成；ViennaRNA 2.7.2 |

环境：Python 3.14.6、NumPy 2.5.2、SciPy 1.18.1、Matplotlib 3.11.1、RDKit 2026.3.5、pandas 3.0.5；ViennaRNA 使用现有隔离目录。没有安装或替换已有科研环境。

## 数值比较与新增数据

[逐列比较报告](recompute_comparison.json) 覆盖综合图引用的 **29 个 CSV、1,401,826 对数值**：16 个文件字节相同，13 个在 `rtol=1e-7, atol=1e-12` 内，0 个超差、0 个缺失、0 个审计故障。比较同时检查表头、行数、文本、缺失与列类型，不把缺失当零。容差仅用于跨运行数值一致性，不能替代各模型专门的守恒、量纲和临床适用性判据。其余 JSON/NPZ 文件记录哈希和来源，不被计入该 CSV 容差结论。

综合图新增的可检查派生表：

- [Task 3 QSSA / GSH Pareto](../../data_omnibus/task3_pareto_qssa.csv)：8 个情景，非支配集合 W02/W06/W08；前向 GSH 半时不代表临床安全。
- [Task 6 绘图热力学网格](../../data_omnibus/task6_thermodynamic_plot_grid.csv)：4 种聚合物 × 301 个载药点。
- [Task 7 肿瘤轨迹分位数](../../data_omnibus/task7_tumor_quantiles.csv)：964 行四组中位数与 10–90% 分位数。
- [Task 8 溶酶体释放通量](../../data_omnibus/task8_lysosomal_payload_flux.csv)：1,442 行，从原 MM 速率式与状态量计算。

独立审阅重新计算了 Task 6 的 FH/Tg（最大差异约 2.8×10⁻¹⁷）、Task 7 的四组 KM（约 1.1×10⁻¹⁶），并核对 Task 8 通量及 Task 10 概率/单位。这些是实现一致性检查，不能替代生物学验证。

## 验证与复现

[108 项测试输出](unittest.txt) · [综合产物检查](validation.json) · [目录与原项目完整性](layout_validation.json) · [视觉与交互检查](visual_review.json)

```bash
# 重新生成归档综合图和页面
python run_ai4pharm_omnibus_suite.py --out work/omnibus

# 完整本机复算并绘图
python run_ai4pharm_omnibus_suite.py --mode recompute --out work/omnibus_fresh

# 读取本次已发布的精选复算数值重新绘图，不启动科学计算
python run_ai4pharm_omnibus_suite.py --source-root data_omnibus/recompute --out work/omnibus_selected

# 独立重做本次 CSV 比较；输出必须为新文件
python tools/compare_omnibus_recompute.py --fresh-root data_omnibus/recompute --execution data_omnibus/recompute/execution.json --out work/omnibus_comparison.json

python tools/validate_omnibus.py --out work/omnibus_validation.json
```

入口拒绝覆盖已有输出、原项目和证据目录；失败会留下状态和日志。最终执行器对超时作进程树终止处理，避免遗留 Task 1 的构象 worker；该故障路径以隔离单元测试检查，本次正常复算未发生超时。原项目的科学驱动未修改。

## 科学解释修正

原始要求作为附件原样保留，不被当成必须“算出来”的结论。两份报告同步补充了可逆 GSH 分支、ASD 第五状态与析晶/再溶解闭合、ADC 十状态与守恒、CNT 的 SI 单位换算。Task 4 第 7 天尚未稳态；Task 9 的状态自由能不是活化势垒，135 个 ROI 截断结果仍不确定，12.625 Å 路径未达到 30 Å；Task 10 在 **37°C WT 基准**中没有超过 1.2 bit 的位置，不能推广至全部温度。

现阶段仅 Task 4 → 5 实际调用同一模型，其余是分别参数化的计算研究。公开结构、执行过的算法、假设参数、代理终点和实验验证需求均分别表述。本次未新增临床预测、实验验证或跨尺度统一校准。
