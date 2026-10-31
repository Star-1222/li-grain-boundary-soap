# 版本关系

整理基线：`450e2ee`。截至该提交，远程只有 `main`，共 16 次线性提交，没有 tag 或 GitHub Release。目录名中的 `v2`、日期和 `new` 是历史命名，不能用于判断整个仓库的最新版本。

## 科学版本按实验区分

| 标识 | 规范位置 | 与其他版本的关系 |
|---|---|---|
| Li 3600 基准 | `Lithium_GrainBoundary_Project/01_GB_Construction/Sigma5_310_Metastable_GB/v2_production/` | 已验证的构建基准；其中 `src/generator/` 是该基准代码入口 |
| Li 早期 2×2 | `Lithium_GrainBoundary_Project/Li_Sigma5_310_metastable_2x2/` | 面内扩大后的早期脚手架；不能用其 NOT RUN 报告判断后续实验状态 |
| Li 阶段一，2026-09 | `Li_Sigma5_310_2x2_阶段一验收_2026-09-05_new/` | 独立验收记录；以包内报告、日志和检查点为准 |
| Li 控压对照，2026-10 | `Li_Sigma5_控压方案对照_2026-10-15/` | 后续敏感性实验；与阶段一并存，协议不能混用 |
| Ni SOAP Stage01–06 | `Ni_SOAP_Method_Exploration_by_Stage/` | 一条方法探索链；汇总报告为规范入口，Stage 报告副本供查阅 |
| 历史副本 | `_archive/`、Ni 的 `90_历史归档/` | 仅用于追溯，不作为新运行入口 |

## Git 历史中的定位点

| 提交 | 内容 |
|---|---|
| `4f3f800` | 初始研究基线，已含 Li / Ni 多套历史材料 |
| `6eb5bd6`–`7d9c4d9` | 工具可移植性、行尾与台账哈希修正 |
| `6f9992a` | 将续算检查点纳入 Git |
| `87b10a6`、`3b8d90f` | 新增首页与学习路线 |
| `2d9345b`、`9e1e730` | SOAP 参数报告及归位 |
| `263edc8`、`468aeb5` | 新增 Ni 分阶段输入输出检索表 |
| `450e2ee` | 新增 10 月控压对照实验 |

本次整理以新增提交记录，不改写以上历史。没有把尚未冻结的方法包装成稳定发布版本。

## 后续提交约定

- `main` 保存可查阅、可追溯的研究状态；每次改动在独立分支完成，通过 PR 审查后合并。
- 一次提交只表达一个目的，例如 `docs: simplify project navigation`、`fix: include barostat experiment in ledger`。
- 新实验按目的和日期命名，例如 `exp/li-barostat-2026-11-01`。先记录输入来源、协议、环境，再提交结果。
- 成熟的里程碑需要 tag 时，采用描述性的实验标识，如 `li-phase1-2026-09`；先核对对应提交和数据完整性，再创建 tag。
- 已完成实验目录保留原始证据。修正结论时注明依据，重新运行时另建目录，避免产生 `final_new_v3` 一类名称。

日常操作见[维护说明](MAINTENANCE.md)。
