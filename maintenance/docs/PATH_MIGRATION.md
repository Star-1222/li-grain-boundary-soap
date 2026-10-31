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
