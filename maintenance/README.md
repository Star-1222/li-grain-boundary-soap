# 仓库维护与资料

这里集中放置跨主线文档、台账、工具和工具测试。研究内容从[晶界构造](../grain_boundary/README.md)和[SOAP 方法](../soap/README.md)进入。

| 内容 | 入口 |
|---|---|
| 工作流与版本 | [工作流](docs/WORKFLOWS.md) · [版本关系](docs/VERSIONS.md) |
| 学习与评估 | [学习路线](docs/STUDY_GUIDE.md) · [项目评估](docs/PROJECT_REVIEW.md) · [研究记录](docs/RESEARCH_STATUS.md) |
| 文件和缺失数据清单 | [产物台账](records/ARTIFACT_LEDGER.md) · [数据恢复清单](records/DATA_MANIFEST.md) |
| 维护工具 | [tools](tools/README.md) |
| 提交与台账操作 | [维护说明](docs/MAINTENANCE.md) |

在仓库根目录运行：

```bash
python maintenance/tools/prepare_soap_workspace.py --dest runs/ni-soap-001
python maintenance/tools/artifact_ledger.py
python -m unittest discover -s maintenance/tests -v
```

数据恢复清单里的路径均相对于仓库根目录。迁移前后的对应关系见[目录迁移表](docs/PATH_MIGRATION.md)。
