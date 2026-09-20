# Task 8–10 交付核对

[交付说明](../TASK8_10_DELIVERY.md) · [仓库首页](../../README.md)

这些记录说明本机实际执行、复现和文件核对的范围，不代表实验或临床模型认证。

| 记录 | 内容 |
|---|---|
| [summary.json](summary.json) | 测试数量、运行环境、图像检查及证据限制 |
| [unittest.txt](unittest.txt) | 完整仓库 96 项测试的终端结果 |
| [reproductions.json](reproductions.json) | 三个冻结驱动独立复算后的逐文件比较；仅排除运行耗时元数据 |
| [layout.json](layout.json) | 十个项目入口、Markdown 链接、清单及历史迁移核对 |
| [delivery.json](delivery.json) | 三份驱动、六份报告、12 张 300 DPI 图的文件与哈希检查 |
| [prior_artifacts.json](prior_artifacts.json) | 原有七个项目的 170 个文件与基线提交的逐字节比较 |
| [staged_integrity.json](staged_integrity.json) | 实际待提交 Git blob 与项目清单哈希一致性 |

独立复算使用 Python 3.12 化学环境；完整测试使用已有 Python 3.14 环境。ViennaRNA 2.7.2 以官方 wheel 安装于被忽略的工作目录，未替换既有环境。运行时间和 Python 补丁版本只描述本次机器；新环境依赖见仓库 requirements。

所有 12 张图由各任务执行者检查，根任务另复核 ADC 图 1/4、构象任务图 1/2/4 和 RNA 全部四图，并核对所发现文字遮挡/裁切的修复。图像 DPI 与文件可读性由交付工具另行检查。数值不变量、视觉检查及已执行输入都不保证所有未知输入均无缺陷。

验收期间修复了：ADC 自定义终点标签；口袋触边帧的未知分类和区域敏感性；PRS 图裁切和密度共用色标；RNA PDB 作者编号、配分函数常数一致性、有限 U1 饱和极限、非单调响应交点、分段选择性网格与无耦合对照。
