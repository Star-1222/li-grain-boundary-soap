# SOAP历史归档

`stages/`完整保留早期六阶段脚本、输入、结果与报告。它们反映当时研究路径，不是当前任务清单；历史文件中的“最终”描述当时版本。

- [六阶段汇总报告](stages/00_汇总报告/Ni_SOAP_最终报告.md)
- [文件索引](stages/00_汇总报告/项目文件索引_脚本输出报告_按阶段.md)
- Stage01—06按原目录保留，科学脚本与数据内容未改变。

若确需复现这套历史流程，从仓库根目录运行：

```bash
python maintenance/tools/prepare_soap_workspace.py --dest runs/ni-soap-history-001
```

该命令准备的是旧六阶段流程。当前参数研究从[当前工作入口](../current/README.md)进入；当前结论统一看[研究状态](../../maintenance/docs/RESEARCH_STATUS.md)。

路径迁移见[迁移记录](../../maintenance/docs/PATH_MIGRATION.md)，逐文件原新路径与哈希见[迁移清单](../../maintenance/records/LAYOUT_MIGRATION_2026-11-15.json)。
