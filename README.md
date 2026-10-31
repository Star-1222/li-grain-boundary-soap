# 锂晶界与 SOAP 方法研究

研究 BCC Li Σ5(310)[001] 晶界的构建与热处理，并用 Ni 双晶探索 SOAP 局域结构分析方法。研究结论、实验记录与运行代码分别保存在对应项目中。

## 从这里开始

| 要做什么 | 入口 |
|---|---|
| 看懂两条研究线、选择运行版本 | [工作流与输入输出](docs/WORKFLOWS.md) |
| 区分基准、验收包、对照实验和旧版本 | [版本关系](docs/VERSIONS.md) |
| 从零学习项目 | [学习路线](docs/STUDY_GUIDE.md) |
| 查看结果与待解决问题 | [研究记录](docs/RESEARCH_STATUS.md) · [项目评估](docs/PROJECT_REVIEW.md) |

## 项目入口

| 目录 | 用途 |
|---|---|
| [Lithium_GrainBoundary_Project](Lithium_GrainBoundary_Project/README.md) | Li 构建基准与早期 2×2 工作流 |
| [Li 阶段一验收包](Li_Sigma5_310_2x2_阶段一验收_2026-09-05_new/README.md) | 14400 原子实验的阶段一结果与续算检查点 |
| [Li 控压对照](Li_Sigma5_控压方案对照_2026-10-15/README.md) | 同一淬火态的 x-only / aniso 对照 |
| [Ni SOAP 方法探索](Ni_SOAP_Method_Exploration_by_Stage/README.md) | Stage01–06 脚本、报告与输出 |
| [_archive](_archive/README.md) | 历史副本，仅供追溯 |

## 常用操作

在仓库根目录运行，工具要求 Python 3.10+：

```bash
python tools/prepare_soap_workspace.py --dest runs/ni-soap-001
python tools/artifact_ledger.py
python -m unittest discover -s tests -v
```

第一条只准备独立的 SOAP 运行副本，具体执行顺序见[工作流](docs/WORKFLOWS.md)。台账工具另需 NumPy。

完整产物见 [ARTIFACT_LEDGER.md](ARTIFACT_LEDGER.md)；克隆中缺失的大文件及校验和见 [DATA_MANIFEST.md](DATA_MANIFEST.md)。环境、台账更新和提交约定见[维护说明](docs/MAINTENANCE.md)。
