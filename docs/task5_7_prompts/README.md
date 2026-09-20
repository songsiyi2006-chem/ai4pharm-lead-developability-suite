# Task 5–7 原始任务与实现说明

此处保存用户提供的三份任务附件。它们记录原始需求，不作为已校验的公式、临床证据或项目结果。

- [Task 5 原始任务](TASK5_ORIGINAL.md)：CYP DDI/TDI。
- [Task 6 原始任务](TASK6_ORIGINAL.md)：ASD 相行为与过饱和。
- [Task 7 原始任务](TASK7_ORIGINAL.md)：肿瘤免疫 QSP 与组合效应。

实现对应任务目录中的脚本与双语报告。需要修正的原始表述如下：

| 原始任务中的问题 | 本次实现采用的解释 |
|---|---|
| 把 DDI 基础 R 比值与 AUCR 混用；给出的净效应公式可能随抑制增强而降低 AUCR | 分别输出基础筛查和机制模型预测，使用有明确清除与肠道逃逸假设的公式，并检验无抑制和强抑制极限 |
| 把强抑制剂类别直接等同于禁忌；将二维图绿色称为安全 | 模拟类别与筛查结果不产生临床禁忌或安全建议 |
| 溶出通量未明确除以溶出体积；CNT 指数的能量量纲不完整 | 用有限固体质量守恒模型，并将成核势垒转成无量纲指数 |
| 用相对湿度直接推出水含量和货架寿命 | 明示水吸附假设、VFT 适用域和未标定的诱导时间；不报告经过验证的货架寿命 |
| 将 Simeoni 指数到线性增长称为承载容量模型 | 区分渐近线性增长和 logistic 承载容量；每级驻留时间为 τ 时三段总均值为 3τ |
| 将 Bliss 正值等同于真实抗癌协同 | 仅称相对于所选模型零假设的超额效应；Loewe 仅在单药曲线支持的效应区间求逆 |
| 将肿瘤质量超过基线 20% 称为 RECIST 1.1 | 保留明确命名的模拟质量进展终点；它不等于临床影像学 RECIST 或临床 PFS |

DDI 筛查依据核对了 [FDA ICH M12 最终指南](https://www.fda.gov/regulatory-information/search-fda-guidance-documents/m12-drug-interaction-studies)及 [FDA 抑制剂分类说明](https://www.fda.gov/drugs/drug-interactions-labeling/healthcare-professionals-fdas-examples-drugs-interact-cyp-enzymes-and-transporter-systems)。真实 RECIST 的目标病灶进展涉及相对研究期最小径线和的增加及绝对增加条件，参见 [EORTC RECIST 1.1](https://recist.eortc.org/recist-1-1/)。更完整的模型公式、来源、参数和限制在各任务报告中逐项说明。

原始任务中的“无 bug”“临床预测”“优化制剂”等目标不被当作已经成立的验收结论。交付记录只报告实际执行的计算、测试、图表检查和当前证据边界。
