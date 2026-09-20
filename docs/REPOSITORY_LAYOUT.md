# 仓库目录与文件定位

[首页](../README.md) · [项目列表](../projects/) · [整改记录](reorganization/README.md)

## 统一规则

1. 一个研究任务对应一个 `projects/taskNN_topic/` 目录，编号补齐两位，便于排序。
2. 项目内 `README.md` 说明研究目的，并列出代码、输入、数据、图表、双语报告与复现命令。
3. 已有数据目录名和结果文件名保留，便于对照原技术报告。项目内部的组织差异在各 README 中明确解释。
4. 共享依赖位于根目录 `requirements.txt`；项目锁定的历史依赖快照随项目保存。
5. `tests/`、`tools/`、`.github/` 为仓库公共设施。`docs/archive/` 保存历史证据，`work/` 只用于本地复算。

## 十个项目

| 当前目录 | 原位置 | 主要文件 |
|---|---|---|
| [task01_lead_developability](../projects/task01_lead_developability/) | 根目录 Task 1 文件和 `data/`、`results_task1/`、`figures_task1/` | MPO/ADMET 代码、结构来源、模型卡、报告和验证说明 |
| [task02_tpd](../projects/task02_tpd/) | 根目录 Task 2 文件和 `data_task2/`、`figures_task2/` | TPD 代码、构象与动力学数据、双语报告、执行日志 |
| [task03_covalent_kinetics](../projects/task03_covalent_kinetics/) | 根目录 Task 3 文件和 `data_task3/`、`figures_task3/` | 共价动力学代码、参数、结构、结果和双语报告 |
| [task04_pbpk](../projects/task04_pbpk/) | `task4_pbpk/` | PBPK 代码、输入模板、时序结果与报告 |
| [task05_cyp_ddi](../projects/task05_cyp_ddi/) | `task5_cyp_ddi/` | CYP 代码、六药输入、`results/`、双语报告与 `TASK.md` |
| [task06_asd](../projects/task06_asd/) | `task6_asd/` | ASD 代码、`inputs/`、`results/`、图表、双语报告与 `TASK.md` |
| [task07_qsp](../projects/task07_qsp/) | `task7_qsp/` | QSP 代码、虚拟队列、轨迹、图表、双语报告与 `TASK.md` |
| [task08_adc](../projects/task08_adc/) | 本次新增 | ADC 代码、DAR/HIC/细胞释放/扩散结果、图表及双语报告 |
| [task09_cryoem_allostery](../projects/task09_cryoem_allostery/) | 本次新增 | 公开蛋白结构、合成系综与密度、MSM/几何/PRS、双语报告 |
| [task10_rna_splicing](../projects/task10_rna_splicing/) | 本次新增 | 公开 RNA 结构、配分函数/几何/热力学结果、图表及双语报告 |

完整的逐文件旧→新映射保存在 [path_map.json](reorganization/path_map.json)。该历史映射覆盖迁移时的 Task 1–7；Task 8–10 是随后新增项目。旧路径不再作为运行入口，不保留重复的源码副本。

## 运行路径

命令统一从仓库根目录执行。各项目 README 的示例为新的计算显式指定 `work/` 下的输出目录，避免覆盖冻结结果。Task 1、Task 2、Task 3 的 `--output-dir` 默认值是脚本所属项目目录；建议复算时仍显式指定新目录。Task 5–10 要求新的输出位置。

```bash
python projects/task01_lead_developability/run_task1_mpo_admet_developability.py --conformers 0 --output-dir work/task01_quick
python projects/task03_covalent_kinetics/run_task3_covalent_kinetics_residence_time.py --self-test-only
python projects/task05_cyp_ddi/run_task5_cyp_ddi_mechanism_based_inhibition.py --self-test
python -m unittest discover -s tests -v
python tools/validate_repository_layout.py --out work/layout_validation.json
```

Task 5 的导入路径已指向 `projects/task04_pbpk/`。从单个项目目录执行时，使用该目录内的脚本文件名即可；Task 5 仍需保留相邻 Task 4 目录。

## 历史与当前校验

科学结果和图表按原始字节迁移，本次没有用新假设替换已有计算。旧首页、原始清单和验收日志保存在 [archive](archive/README.md)，里面的路径、代码哈希和测试计数描述当时的提交。

Task 1–7 项目清单中的 `layout_migration` 记录迁移基线、当前代码哈希和发生变更的文档/代码原哈希。其原 `script_sha256` 仍表示生成科学结果时的代码，不能被当作本次重新计算的证明。当前文件校验和已适配新目录。

Task 8–10 直接按新目录组织；清单中的 `script_sha256` 对应各自实际执行的驱动脚本，`files_sha256` 记录项目内产物。它们不使用旧迁移记录伪装为历史项目。
