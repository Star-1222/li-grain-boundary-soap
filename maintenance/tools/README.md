# 仓库工具

所有命令从仓库根目录运行；Python 3.10+。

| 工具 | 用途 | 依赖 |
|---|---|---|
| `prepare_soap_workspace.py` | 将 Ni 六阶段脚本和 CFG 复制到独立运行目录，生成执行顺序与来源清单 | 标准库 |
| `artifact_ledger.py` | 生成 `ARTIFACT_LEDGER.md/csv`，记录项目文件、用途、SHA256 与内容分析 | NumPy |
| `build_artifact_catalog.py` | 台账分类与扫描的底层模块；日常更新使用上面的入口 | 标准库 |
| `compute_gb_energy.py` | 从历史 LAMMPS 日志复算晶界过剩能 | 标准库 |
| `install_git_hook.ps1` | 可选：提交前自动刷新并暂存台账 | PowerShell、Git、Python、NumPy |

## 准备 SOAP 运行副本

```bash
python maintenance/tools/prepare_soap_workspace.py --dest runs/ni-soap-001
```

目标目录必须不存在。工具只复制脚本、两份原始 CFG 和依赖列表，不复制旧结果、不启动计算。新目录中的 `SOURCE_MANIFEST.json` 记录来源与 SHA256，`RUN_ORDER.md` 给出逐阶段命令。计算环境与执行限制见[工作流](../docs/WORKFLOWS.md)。

## 台账与能量复算

```bash
python maintenance/tools/artifact_ledger.py
python maintenance/tools/compute_gb_energy.py
```

Windows 也可以运行 `maintenance/refresh_ledger.cmd`。台账是派生物，生成前应先暂存源文件；大小与哈希取暂存区 blob，内容分析取工作区。生成后把两份台账一起暂存。完整操作见[维护说明](../docs/MAINTENANCE.md)。

`DATA_MANIFEST.md/csv` 保存未入库数据的恢复清单，只有在数据完整的原机器上才用 `--data-manifest` 更新。全新克隆保留现有清单。

工具默认使用其所在仓库，旧台账工具支持用 `RESEARCH_ROOT` 指向其他根目录；`LEDGER_MIRROR` 为空时不生成镜像副本。

## 可选 Git 钩子

在仓库根目录执行一次：

```powershell
powershell -ExecutionPolicy Bypass -File maintenance\tools\install_git_hook.ps1
```

该脚本会写入 `.git/hooks/pre-commit`，安装前自行保留已有自定义钩子。提交时它刷新并自动暂存台账；生成失败仍放行，所以不能代替手工检查。钩子不随克隆自动安装。

## 验证工具

```bash
python -m unittest discover -s maintenance/tests -v
```

测试验证目录覆盖、运行副本与来源完整性、覆盖保护和文档链接。台账的覆盖与稳定性测试另需 NumPy。它们不认证历史 MD / SOAP 结果。

现行晶界运行与起点保存从 [grain_boundary/current](../../grain_boundary/current/README.md) 进入；能量复算工具默认查阅历史基准日志，不代表现行候选起点已经验收。
