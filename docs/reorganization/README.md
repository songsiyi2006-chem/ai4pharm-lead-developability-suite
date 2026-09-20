# 目录整改记录 · 2026-09-20

[返回首页](../../README.md) · [当前目录说明](../REPOSITORY_LAYOUT.md)

迁移基线为 `3785df88d244bf276de7f9b6e4d1d74995b530b9`。该基线的 190 个版本化文件均在 [path_map.json](path_map.json) 中逐项记录；其中 177 个移动到对应项目或历史档案目录。新增的项目 README 和公共导航单独纳入本次提交。

本次工作调整目录、代码导入、默认输出位置、文档链接与校验路径。已有科学数据、结构文件、输入、图表和科学汇总保持原始字节。没有新增科学实验或借目录调整更新数值结果。

| 检查 | 记录 |
|---|---|
| 全仓库测试及必要 CLI 检查 | [validation.json](validation.json) |
| 全仓库测试终端输出 | [unittest.txt](unittest.txt) |
| 项目入口、文件迁移、数据保留、链接和哈希 | [layout_validation.json](layout_validation.json) |

`path_map.json.previous_sha256` 来源于基线 Git 对象。当前源代码哈希保存在各项目 manifest 的 `layout_migration.current_sources_sha256`，原始执行哈希及完整旧 manifest 均保留。已改变的代码/文档哈希刷新属于元数据维护，不代表再次运行完整科学计算。

目录检查工具为 [validate_repository_layout.py](../../tools/validate_repository_layout.py)，同时接入 CI。历史计算比较工具 [validate_task5_7_delivery.py](../../tools/validate_task5_7_delivery.py) 已适配新任务位置。
