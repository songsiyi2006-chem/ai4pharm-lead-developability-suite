# Task 09 · Synthetic structural heterogeneity and cryptic-pocket audit

[Return to projects](../README.md) · [Original instructions](TASK.md) · [中文报告](CRYOEM_CRYPTIC_POCKET_REPORT_ZH.md) · [English report](CRYOEM_CRYPTIC_POCKET_REPORT_EN.md)

使用公开 T4 溶菌酶 L99A 的 4W51/4W59 晶体坐标，执行 500 个合成构象、非线性嵌入、高斯密度体素、可逆五态 MSM、网格口袋/SASA 和 ANM/PRS。该项目是可复现方法基准：没有真实 cryo-EM 颗粒、cryoDRGN 训练、MD、配体结合或癌症靶点验证。

| 文件夹/文件 | 内容 |
|---|---|
| [单文件驱动](run_task9_cryoem_cryptic_pocket_allostery.py) | 四模块计算、4图、双语报告、数值检查 |
| [inputs](inputs/) | 原始 PDB、来源/哈希、可编辑参数 |
| [results](results/) | 500坐标/latent、三维density、MSM bootstrap、逐帧几何、ROI/网格敏感性、PRS |
| [figures_task9](figures_task9/) | 4张300 DPI计算图 |
| [verification.json](verification.json) | 数值恒等式与收敛检查 |
| [manifest.json](manifest.json) | 代码及产物SHA256 |

默认结果的开态只出现6帧；人口自由能估计为2.057 kcal/mol，95%参数自助区间1.10–3.33，不能当作真实动力学势垒。ROI内空腔为76.8–209.7 Å³；135帧触边，完整体积与阈值分类未知，未证明成药窗口，也未排除其存在。最小接触耦合路径106→11，实际端距12.625 Å，未达到30 Å。

复算需要已安装的 Python、NumPy、SciPy、Matplotlib；不下载模型、不在线取数。**输出目录必须是新目录**，CLI不会覆盖此处归档。

```powershell
python projects/task09_cryoem_allostery/run_task9_cryoem_cryptic_pocket_allostery.py --out work/task9_reproduce
python projects/task09_cryoem_allostery/run_task9_cryoem_cryptic_pocket_allostery.py --write-example work/task9_config.json
python projects/task09_cryoem_allostery/run_task9_cryoem_cryptic_pocket_allostery.py --config work/task9_config.json --out work/task9_custom
python projects/task09_cryoem_allostery/run_task9_cryoem_cryptic_pocket_allostery.py --self-test
```

参数调整不等于实验校准。若某状态没有采样，密度计算会明确失败，不会自动换随机种子得到目标结果。原始PDB采用A链alternate A，配体只提供ROI中心并从蛋白体积/密度计算中移除。新增项目与既有PBPK、QSP之间尚无已校准的直接耦合。
