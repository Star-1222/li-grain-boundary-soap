# 版本关系

整理基线：`450e2ee`。截至该提交，远程只有 `main`，共 16 次线性提交，没有 tag 或 GitHub Release。目录名中的 `v2`、日期和 `new` 是历史命名，不能用于判断整个仓库的最新版本。

当前目录按两条研究主线组织，旧路径对照见[目录迁移表](PATH_MIGRATION.md)。

## 当前与历史

- **现行协议与代码：** [grain_boundary/current](../../grain_boundary/current/README.md)，标准采用淬火后 10 K aniso 松弛、三项稳定性审查、保存起点。
- **标准起点状态：** [STATUS.json](../../grain_boundary/current/STATUS.json)，尚未登记通过新审查程序的结构。
- **历史证据：** [grain_boundary/history](../../grain_boundary/history/README.md)，保留 3600 基准、早期 2×2、阶段一验收、控压对照和旧报告；均不作为现行默认版本。
- **SOAP：** [soap](../../soap/README.md) 保持自己的六阶段方法探索流程。

历史参数表和旧报告只描述对应运行。后续使用现行 model.json 和 PROTOCOL.md，不能按目录日期或 final/new 字样选择流程。

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
