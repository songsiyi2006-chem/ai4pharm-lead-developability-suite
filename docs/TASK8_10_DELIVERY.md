# Task 8–10：ADC、构象变构与 RNA 剪接

[返回仓库首页](../README.md) · [十个项目](../projects/) · [原始需求解释](SCIENTIFIC_NOTES.md)

执行日期：2026-09-20。此次工作继续使用公开数据和本机 CPU；没有新增湿实验、临床样本或 HPC 作业。每项任务分别保存驱动脚本、原始任务、输入与参数、数值结果、四张 300 DPI 图和中英文技术报告。

| 任务 | 计算内容 | 报告与代码 |
|---|---|---|
| 08 · ADC 工程 | 二硫键还原与 DAR 0–8 偶联、HIC 分离、受体转运与溶酶体释放、球形组织扩散和摄取 | [项目](../projects/task08_adc/) · [中文](../projects/task08_adc/ADC_TRANSLATIONAL_ENGINEERING_REPORT_ZH.md) · [EN](../projects/task08_adc/ADC_TRANSLATIONAL_ENGINEERING_REPORT_EN.md) |
| 09 · 构象与变构 | 公开蛋白结构起点、合成构象和模拟三维密度、五态可逆 MSM、口袋几何及 ANM/PRS 接触网络 | [项目](../projects/task09_cryoem_allostery/) · [中文](../projects/task09_cryoem_allostery/CRYOEM_CRYPTIC_POCKET_REPORT_ZH.md) · [EN](../projects/task09_cryoem_allostery/CRYOEM_CRYPTIC_POCKET_REPORT_EN.md) |
| 10 · RNA 剪接 | 最近邻配分函数和配对概率、公开 NMR 结构几何、四个假设配体的耦合热力学及 U1 占有率 | [项目](../projects/task10_rna_splicing/) · [中文](../projects/task10_rna_splicing/RNA_TARGETED_CADD_REPORT_ZH.md) · [EN](../projects/task10_rna_splicing/RNA_TARGETED_CADD_REPORT_EN.md) |

## 本次运行规模和主要结果

- Task 8：4 个偶联条件、24 组各 1,000 个抗体的 Gillespie 模拟、2 个细胞/组织通透性情景、2 次细化及 4 次运输参数敏感性计算，共 23,200 行空间数据。默认平均 DAR 为 3.99932；假设 HIC 收集窗的 DAR 4 纯度为 96.915%。72 h 邻近抗原阴性细胞的平均存活代理为 0.958800 / 0.999848；这些数值描述高/低通透性的预设模型情景。原始数表见 [Task 8 数据](../projects/task08_adc/data/)。
- Task 9：500 个合成构象、200 次 MSM 参数自助抽样、500 帧口袋/SASA、6 次网格及 9 次区域大小敏感性计算；ANM 为 492×492，PRS 为 164×164。开态仅采到 6 帧；135 帧口袋触边，完整体积与阈值分类保持未知。106→11 的接触耦合路径端距为 12.625 Å，没有达到原任务要求的 30 Å。详见 [Task 9 数据与限制](../projects/task09_cryoem_allostery/results/summary.json)。
- Task 10：14 次 ViennaRNA 配分函数、40 个公开 NMR 构象、560 行沟槽朝向几何描述符、976 个剂量平衡点及 120 个假设敏感性点。37 °C 默认序列的位置熵最大值为 0.925 bit，没有达到预设的 1.2 bit 门槛。非单调响应保留半动态范围交点，单一增强型 EC50 不适用；选择性代理按连续合格网格段保存，并标记触及扫描上限的截断。详见 [Task 10 结果](../projects/task10_rna_splicing/results/results_summary.json)。

## 数据身份与结论范围

Task 8 的化学和运输速率属于机制情景。DAR 4 是分离方法分析对象，不是所有 ADC 的通用最优值；膜通透性比较不等于真实 DXd/MMAF 制剂之间的定量药效比较。组织浓度越过假设阈值不能独立证明旁观者细胞已死亡。

Task 9 的公开起点是 T4 溶菌酶 L99A 的 [4W51](https://www.rcsb.org/structure/4W51) 和 [4W59](https://www.rcsb.org/structure/4W59) 晶体结构。它们是几何锚点；中间构象、模拟密度、状态布居和时间标尺需要独立标注。这里没有真实 cryo-EM 粒子重构或 cryoDRGN 训练；也没有因为任务预设了数值范围就强行给出对应口袋体积、开放代价或通路距离。

Task 10 的 [6HMI](https://www.rcsb.org/structure/6HMI) 与 [6HMO](https://www.rcsb.org/structure/6HMO) 是 U1 snRNA/SMN2 剪接位点的 NMR RNA 结构，其中结合配体是 SMN-C5。NMR 模型序号不是时间轴，也不是可直接统计热力学布居的轨迹。ViennaRNA 计算使用明确说明的序列模型；三级翻转自由能和配体参数另作假设。平衡 U1 占有率不等于经实验校准的外显子包含率，不能从中获得临床安全剂量窗口。

## 复现与维护

从仓库根目录进入各项目 README，按其命令将新计算写入 `work/` 内的新目录。公开结构输入已随项目保存，常规复算无需再次下载。Task 10 使用官方 ViennaRNA 2.7.2，已列入共享依赖；本次 Windows 运行仅在被 Git 忽略的工作目录增加独立模块，没有替换既有化学环境。

```bash
python -m unittest discover -s tests -v
python tools/validate_repository_layout.py --out work/layout_validation.json
python tools/validate_task8_10_delivery.py --out work/task8_10_validation.json
```

第二条命令同时核对历史 Task 1–7 迁移产物的字节保留；第三条核对本次精确要求的脚本、双语报告、12 张图及项目清单。软件测试、图像检查和哈希核对不代替实验或临床验证。

本次完整仓库测试 **96 项通过**，其中新增 Task 8/9/10 分别为 12/14/14 项。三个驱动另在独立输出目录复算，按文件核对结果；原有 Task 1–7 的 170 个项目文件与上次提交逐字节一致。详细证据、排除的运行耗时元数据及图像检查范围见 [验收记录](task8_10_validation/README.md)。
