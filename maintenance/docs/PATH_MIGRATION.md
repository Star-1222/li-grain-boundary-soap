# 目录迁移

2026-10-31 按两条研究主线收拢根目录。旧路径来自提交 `4c3a2d8`，下表的新路径均相对仓库根目录。

| 原路径 | 新路径 |
|---|---|
| `Lithium_GrainBoundary_Project` | `grain_boundary` |
| `Li_Sigma5_310_2x2_阶段一验收_2026-09-05_new` | `grain_boundary/experiments/phase1_2026-09` |
| `Li_Sigma5_控压方案对照_2026-10-15` | `grain_boundary/experiments/barostat_2026-10` |
| `Ni_SOAP_Method_Exploration_by_Stage` | `soap` |
| `_archive` | `grain_boundary/archive` |
| `docs` | `maintenance/docs` |
| `tools` | `maintenance/tools` |
| `tests` | `maintenance/tests` |
| `ARTIFACT_LEDGER.md` | `maintenance/records/ARTIFACT_LEDGER.md` |
| `ARTIFACT_LEDGER.csv` | `maintenance/records/ARTIFACT_LEDGER.csv` |
| `DATA_MANIFEST.md` | `maintenance/records/DATA_MANIFEST.md` |
| `DATA_MANIFEST.csv` | `maintenance/records/DATA_MANIFEST.csv` |
| `refresh_ledger.cmd` | `maintenance/refresh_ledger.cmd` |

科研脚本、输入、结构、日志与检查点按原始内容迁移。历史报告内的原机器绝对路径保留，代表当时环境；运行前以现行工作流为准。DATA_MANIFEST 只更新路径，保留原记录的大小与校验和。

## 标准流程收敛（2026-10-31）

现行入口为 `grain_boundary/current/`。旧目录迁入历史区：

| 原路径 | 新路径 |
|---|---|
| `grain_boundary/01_GB_Construction` | `grain_boundary/history/construction_baseline_3600` |
| `grain_boundary/Li_Sigma5_310_metastable_2x2` | `grain_boundary/history/prototype_2x2` |
| `grain_boundary/experiments` | `grain_boundary/history/experiments` |
| `grain_boundary/Reports` | `grain_boundary/history/reports` |
| `grain_boundary/07_Documentation` | `grain_boundary/history/presentations` |
| `grain_boundary/archive` | `grain_boundary/history/legacy` |

## 当前研究与历史归档分离（2026-11-15）

| 原路径 | 当前路径 |
|---|---|
| `soap/Stage01_*`至`soap/Stage06_*`及阶段检索表 | `soap/history/stages/`下同名目录与文件 |
| `soap/00_汇总报告` | `soap/history/stages/00_汇总报告` |
| `soap/06_环境与实现` | `soap/history/stages/06_环境与实现` |
| `soap/90_历史归档` | `soap/history/stages/90_历史归档` |
| 未入库文献原路径`soap/05_文献` | 恢复时使用`soap/history/stages/05_文献` |

当前SOAP研究入口是 `soap/current/`；原 `maintenance/tools/prepare_soap_workspace.py`转为明确的历史六阶段准备入口。科学脚本、输入、数值结果和图像按字节保留；未入库数据清单仅迁移路径，保留大小和SHA256。

逐文件旧路径、新路径和迁移前后哈希见[迁移清单](../records/LAYOUT_MIGRATION_2026-11-15.json)。
