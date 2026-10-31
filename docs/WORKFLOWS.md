# 工作流与输入输出

Li 物理实验与 Ni SOAP 方法验证是两条研究线。Ni Stage01–06 与 Li 的 Stage0–4 编号彼此独立；Ni 结果不能直接当作 Li 的校准结果。

## Li：按实验目的选择入口

| 单元 | 输入 → 处理 → 输出 | 当前用途 |
|---|---|---|
| 3600 原子基准 | 体相晶格常数 → 构建、检查、极小化、300 K 退火、淬火、最终极小化 → `data/final/stage4_final.data` | 构建与几何回归基准 |
| 早期 2×2 脚手架 | `stage0_build_validate/` → `lammps/` Stage0–4 | 14400 原子初始构建与旧协议参考；报告的 NOT RUN 是当时快照 |
| 9 月阶段一验收 | 初始构型 → 极小化、363 K 退火 25 ps、10 K/ps 淬火、10 K 保持 5 ps → `p1_held.data` / `.restart` | 已完成实验的报告、日志与续算起点 |
| 10 月控压对照 | 同一 `stage3b_quenched.restart` → x-only / aniso 各保持 5 ps → 热力学与结构分析 | 检验控压方式；不能直接覆盖前一实验的验收口径 |

- 基准代码：[v2_production](../Lithium_GrainBoundary_Project/01_GB_Construction/Sigma5_310_Metastable_GB/v2_production)。生成器在 `src/generator/`，LAMMPS 输入在 `lammps/`，已有结果在 `data/`。
- 早期脚手架：[Li_Sigma5_310_metastable_2x2](../Lithium_GrainBoundary_Project/Li_Sigma5_310_metastable_2x2)。它与验收包的协议、目录命名不同，不能拼成一次运行。
- 已完成阶段一：[验收包入口](../Li_Sigma5_310_2x2_阶段一验收_2026-09-05_new/README.md)。续算必须按该包的报告选择检查点。
- 控压对照：[实验说明](../Li_Sigma5_控压方案对照_2026-10-15/README.md)。复现说明依赖原 WSL 运行目录；全新克隆还需恢复上游输入，不能直接照抄本机绝对路径。

构建基准的 Python 调用以当前目录结构为准，历史双语报告中的 `outputs/`、`generator/` 路径是重排前布局：

```bash
cd Lithium_GrainBoundary_Project/01_GB_Construction/Sigma5_310_Metastable_GB/v2_production
python -m pip install -r requirements.txt
# Bash / WSL
export PYTHONPATH="$PWD/src"
python -m generator.generate --a 3.506186475312 --outdir ../../../../runs/li-build-001
```

这条命令只构建新模型。完整 MD 需要 LAMMPS + OpenKIM、MPI 和协议对应的输入；本次整理没有认证旧 runner 的断点恢复逻辑。新计算应使用新的输出及日志目录，避免已有 `.done` 标记跳过计算。

## Ni：先准备运行副本，再按依赖执行

归档目录按 Stage 拆分了脚本和输出，部分脚本仍按“同目录模块 + `outputs/`”读写。准备工具把六阶段脚本和两份 CFG 复制到一个新目录，原始归档留在原位：

```bash
python tools/prepare_soap_workspace.py --dest runs/ni-soap-001
python -m pip install -r Ni_SOAP_Method_Exploration_by_Stage/06_环境与实现/requirements_stages.txt
cd runs/ni-soap-001
```

目标目录必须不存在。工具会记录每个复制文件的原路径与 SHA256，并生成 `RUN_ORDER.md`。它不复制历史输出，也不启动计算；后续结果写到新目录下的 `outputs/`。

| 顺序 | 脚本，按从左到右运行 | 主要依赖与结果 |
|---|---|---|
| 01 | `inspect_structures.py` → `soap_test.py` → `stage1_qc.py` | 两份 Ni CFG → 数据检查与 SOAP 初测 |
| 02 | `stage2_profiles.py` → `stage2_soap.py` → `stage2_5_audit.py` → `stage2_75_full.py` | CFG、代表原子 → SOAP 矩阵、原子元数据、距离审计 |
| 03 | `stage3a_lae_threshold.py` → `stage3b_lae_ler.py` | Stage02 SOAP 矩阵与元数据 → 阈值稳定性、LAE/LER |
| 04 | `stage4a_continuous_disorder.py` → `stage4b_disorder_energy.py` | 上游矩阵与分类 → 连续无序度、能量关联 |
| 05 | `stage5a_physical_spatial_residual.py` → `stage5b_soap_direction.py` | 上游逐原子结果与 CFG → 控制物理量后的残差与方向分析 |
| 06 | `stage6_final_analysis.py` | Stage02、04、05 输出 → 最终评估 |

每条脚本命令为 `python <脚本名>`。先确认 Stage01 检查通过，再进入计算量更大的阶段。准备副本解决归档后的路径关系；完整六阶段尚未在本次整理中重跑。依赖清单按导入模块整理，未锁定包版本，复现实验还需记录实际环境与参数。

已有报告、数据表和逐文件索引从 [Ni 项目入口](../Ni_SOAP_Method_Exploration_by_Stage/README.md) 查阅。若从中间阶段开始，需恢复该阶段的全部上游数据，缺失项见 [DATA_MANIFEST](../DATA_MANIFEST.md)。
