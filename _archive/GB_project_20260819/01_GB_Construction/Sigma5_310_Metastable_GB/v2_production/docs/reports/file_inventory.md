# 文件清单与分类 file_inventory.md

> 生成日期：2026-08-12
> 扫描范围（3 个根目录）：
> 1. 工作区 `C:\Users\1\Documents\Codex\2026-08-01\bcc-li-5-001-sigma5-310`（153 个文件）
> 2. `C:\Users\1\Documents\Li_Project`（55 个文件，不含 venv；venv 另有 9215 个文件）
> 3. `C:\Users\1\Documents\Li_GB_Project`（52 个文件）
>
> **目标根目录假设（待确认）**：新 `Li_Project/` = `C:\Users\1\Documents\Li_Project\`
> 分类标准：A=FINAL（生产），B=ARCHIVE（旧版但有价值），C=TEMP（可日后删除）
> 移动前必须经您确认；本文件只做登记，不删除任何科研文件。

---

## 0. 判断依据（关键事实）

| 事实 | 说明 |
|---|---|
| v2_final 是权威生产包 | 63 项测试通过；Stage 0–4 全部 `.done`；3600 原子；末速度全零；fmax=3.55e-7 eV/Å；Pxx=3.6e-4 bar |
| 旧生成器存在 3 份完全相同的拷贝 | `legacy_results/`、`Li_Project/Li_Bicrystal_Generator/`、`Li_GB_Project/Li_Bicrystal_Generator/` 的 builder/config/geometry/exporter/main.py 的 SHA256 完全一致 |
| 旧流程几何有缺陷 | 源码 `SUPERCELL=(20,12,12)` + 全局 `OVERLAP_CUTOFF=2.3`；实际输出 3224 原子（3600→3224，10.44%）；三斜盒 85.977×60.795×14.892 |
| 新流程 | 镜像法构建；3720→3600（仅删 GB 面重合原子 120 个，3.2258%）；正交盒 66.525×33.263×35.062（a=3.506186475312 Å） |
| 势函数无冲突 | 三处全部使用 `SNAP_ZuoChenLi_2019_Li__MO_732106099012_000` |
| 晶格常数有差异（有意） | 旧：a=3.51（初值）；新：a=3.506186475312（SNAP 平衡值） |

---

## 1. 工作区 `v2_final/` → `production/`（FINAL，81 个文件）

### 1.1 根文件（6）

| 文件 | 作用 | 类别 | 移动目的地 | 理由 |
|---|---|---|---|---|
| README.md | 权威工作流说明 | FINAL | production/README.md | v2_final 主文档 |
| FINAL_STATUS.md | 最终验证状态 | FINAL | production/FINAL_STATUS.md | 记录 Stage 0–4 收敛结果 |
| ENVIRONMENT.txt | 运行环境记录 | FINAL | production/ENVIRONMENT.txt | 复现所需环境信息 |
| SHA256SUMS.txt | 文件完整性校验 | FINAL | production/SHA256SUMS.txt | 校验 v2_final 文件指纹 |
| pytest.ini | pytest 配置 | FINAL | production/pytest.ini | 测试配置 |
| requirements.txt | 依赖清单 | FINAL | production/requirements.txt | 依赖声明 |

### 1.2 `generator/`（10 个，全部 FINAL → production/generator/）

| 文件 | 作用 | 理由 |
|---|---|---|
| config.py | 参数：a、θ、KIM 模型 | 新架构核心 |
| geometry.py | 晶体学核心（旋转矩阵、立方群、逆格点、CellList） | 已修复 24 元立方群 |
| builder.py | 对称镜像法双晶构建 + overlap 删除 | 新几何核心 |
| exporters.py | data/XYZ/CIF/grain_map/删除记录导出 | 导出器 |
| verify.py | 11 项自动几何验证 | 验证器 |
| lammps_io.py | data/restart 解析 | 重导入验证依赖 |
| generate.py | CLI：构建+验证+导出 | 入口 |
| bulk.py | bulk 参考结构 | 参考能计算 |
| reimport_check.py | ASE 重新读入交叉验证 | 导出验证 |
| __init__.py | 包初始化 | 包结构 |

### 1.3 `lammps/`（11 个，全部 FINAL → production/lammps/）

| 文件 | 作用 | 理由 |
|---|---|---|
| in.0_check | Stage 0：势+结构 sanity（run 0） | Stage 0–4 标准流程 |
| in.1_minimize | Stage 1：初始最小化（松弛 x，固定 y/z） | 同上 |
| in.2_anneal | Stage 2：300 K NVT 20 ps | 同上 |
| in.2_anneal_resume | Stage 2 断点续算（run … upto） | 真 restart 续算 |
| in.3_quench | Stage 3：300→1 K 淬火 50 ps | 目标 1 K（禁 0 K NVT） |
| in.3_quench_resume | Stage 3 断点续算 | 真 restart 续算 |
| in.4_final | Stage 4：0 K 最终最小化（fmax 判据） | 力收敛 + Pxx≈0 |
| in.bulk_energy | bulk 参考能 | 参考能脚本 |
| in.dump2data | dump→data 重建（box yes add yes） | 已修复空盒问题 |
| run_all.sh | 全流程驱动（.done 检查点） | 安全中断/续跑 |
| resume.sh | 续跑入口 | 同上 |

### 1.4 `tests/`（14 个，全部 FINAL → production/tests/）

conftest.py、test_builder.py、test_conventional_bcc.py、test_grain_id.py、test_legacy_regression.py、test_no_nvt_zero.py、test_outputs_consistency.py、test_overlap_removal.py、test_periodicity.py、test_sigma5_geometry.py、test_velocities_zero.py、test_verify.py、test_workflow_scripts.py —— 覆盖用户要求的全部测试点（conventional BCC、Σ5 几何、周期性、grain ID、overlap 局部性、输出一致性、末速度零、禁 0K NVT）。

### 1.5 `outputs/`（24 个，全部 FINAL → production/outputs/）

| 文件 | 作用 |
|---|---|
| Sigma5_310.data / .xyz / .cif / .grain_map / _grain_id.xyz | 初始双晶结构（3600 原子） |
| removal_log.csv | 120 条删除记录 |
| params.json | 构建参数（a=3.506186475312） |
| build_report.md / verification_report.md | 构建与验证报告 |
| stage1_minimized.data / .restart | Stage 1 输出 |
| stage2_anneal.dump / stage2_annealed.data / .restart | Stage 2 输出 |
| stage3_quench.dump / stage3_quenched.data / .restart | Stage 3 输出 |
| stage4_final.data / .restart | Stage 4 输出（末速度全零） |
| bulk/bulk_Li.data / .xyz / .cif / bulk_Li_minimized.data / .restart | bulk 参考结构 |

### 1.6 `logs/`（16 个，FINAL → production/logs/）

bulk_energy.lammps、dump2data_validation.lammps、kim.log、stage0_check.done、stage0_check_a0_mpi12.lammps、stage1_minimize.done、stage1_minimize_checkpoint_mpi12.lammps、stage2_anneal.done、stage2_anneal_mpi12.lammps、stage2_anneal_mpi12.screen、stage3_quench.done、stage3_quench_mpi12.lammps、stage3_quench_mpi12.screen、stage4_final.done、stage4_final_mpi12.lammps、stage4_final_mpi12.screen —— 各阶段独立日志与 `.done` 检查点。

### 1.7 `validation/`（1 个）

| 文件 | 作用 | 类别 | 移动目的地 |
|---|---|---|---|
| rebuilt_stage2.data | dump2data 重建验证产物 | FINAL | production/validation/rebuilt_stage2.data |

---

## 2. 工作区 `legacy_results/` → `archive/`（ARCHIVE，30 个文件）

### 2.1 旧生成器源码（9 → archive/legacy_generator/）

| 文件 | 作用 | 类别 | 理由 |
|---|---|---|---|
| builder.py / config.py / geometry.py / exporter.py / main.py | 旧生成器（primitive cell 当 conventional、固定盒旋转、全局删原子） | ARCHIVE | 保留以追溯历史错误 |
| README_generator.md | 旧生成器文档 | ARCHIVE | 同上 |
| README_project.md | 旧项目说明（与 Li_GB_Project/README.md 同内容） | ARCHIVE | 同上 |
| requirements.txt / SETUP.md | 旧环境配置 | ARCHIVE | 同上 |

### 2.2 旧 LAMMPS 工作流（17 → archive/old_lammps/）

annealed.data、anneal_save.in、continue.in、dump2data.in、final_gb.data、final_gb.restart、in.bicrystal_relax、in.quench、in.resume、kim.log、log.lammps、run_relax.sh、Sigma5_310.lmp、Sigma5_310_minimized.data、Sigma5_310_minimized.restart、traj_anneal.dump、traj_quench.dump —— 旧流程结果与脚本（含 `fix nvt temp 300.0 0.0`、`in.resume` 与主脚本相同等问题），作为历史/错误分析保留。

### 2.3 旧输出（4 → archive/old_outputs/）

Sigma5_310.cif、Sigma5_310.lmp、Sigma5_310.xyz、Sigma5_310_grainID.xyz（旧 3224 原子结构）。

---

## 3. 工作区 `v2_obsolete_20260802_150439/` → `archive/debugging_history/`（ARCHIVE，41 个文件）

DO_NOT_USE.txt、README.md、pytest.ini、requirements.txt、generator/（10）、lammps/（10，含 in.resume 与“x 固定/y,z 松弛”的旧策略）、reports/（4）、tests/（13）。
理由：是 v2_final 之前的中间版本（几何核心已修正，但工作流策略与脚本尚未定型）；目录自带 DO_NOT_USE 标记，仅作历史追溯，不作科研数据。

---

## 4. 工作区 `work/`（1 个文件）

| 文件 | 作用 | 类别 | 移动目的地 | 理由 |
|---|---|---|---|---|
| cs5_310_ref_result.json | Node 独立参考验证结果（34 KB） | ARCHIVE | archive/debugging_history/cs5_310_ref_result.json | 验证溯源；建议同时留一份到 production/validation/ |

---

## 5. 工作区根目录 `outputs/`（0 个文件）

空目录。类别 TEMP（空目录，无文件，不移动内容；是否保留目录待确认）。

---

## 6. `C:\Users\1\Documents\Li_Project`（55 个文件，不含 venv）

### 6.1 `Li_Bicrystal_Generator/`（旧生成器完整拷贝，ARCHIVE）

与 Li_GB_Project/Li_Bicrystal_Generator 完全重复（核心 5 个 .py 的 SHA256 相同；另含 utils.py、__init__.py、__main__.py、examples/、__pycache__/、lammps_workflow/、output/）。
→ 与 `archive/legacy_generator/` 去重：**保留一份**（建议保留 Li_GB_Project 那份或工作区 legacy_results 那份），其余拷贝移入 `temp/duplicates/`（不删除）。

### 6.2 `output/`（旧输出，4 个文件）

Sigma5_310.cif/.lmp/.xyz/_grainID.xyz —— 与 legacy_results/output/ 相同。ARCHIVE → archive/old_outputs/（去重，多余拷贝进 temp/duplicates/）。

### 6.3 根目录早期 LAMMPS 环境验证（6 个文件）

final.data、in.test_li、kim.log、log.lammps、README.md、run.sh —— 与 Li_GB_Project/lammps_test/ 重复。ARCHIVE → archive/old_lammps/early_env_test/（或按去重规则）。

### 6.4 `markdown file/test 1.md`

Markdown 入门教程（与项目无关）。TEMP → temp/unrelated/。

### 6.5 `Li_Bicrysmkdir/`（空目录，拼写错误）

TEMP（空目录，无可移动文件）。

### 6.6 `venv/`（Python 3.14 虚拟环境，9215 个文件）

TEMP（环境产物，可再生）。**建议原地保留、不移动**（venv 内含绝对路径，移动会失效）；如确认后想清理，请手动删除并重建，或改用 requirements.txt。

---

## 7. `C:\Users\1\Documents\Li_GB_Project`（52 个文件）

已组织良好的旧项目根：README.md、SETUP.md、requirements.txt、setup_env.sh、.gitignore、docs/（PDF+PPTX 开题材料）、literature/（8 篇文献 + soap_method/）、lammps_test/、Li_Bicrystal_Generator/、notes/（空）。
**建议：保持不动**（不在本次 Li_Project 整理范围内；如需要可后续并入 archive/ 或保留为独立文献目录）。`notes/` 为空目录。

---

## 8. 汇总

| 来源 | 文件数 | FINAL | ARCHIVE | TEMP |
|---|---|---|---|---|
| 工作区 v2_final | 81 | 81 | 0 | 0 |
| 工作区 legacy_results | 30 | 0 | 30 | 0 |
| 工作区 v2_obsolete | 41 | 0 | 41 | 0 |
| 工作区 work | 1 | 0 | 1 | 0 |
| 工作区 outputs | 0 | 0 | 0 | 0（空目录） |
| Li_Project（除 venv） | 55 | 0 | 约 27 | 约 28（含重复拷贝） |
| Li_Project venv | 9215 | 0 | 0 | 9215（原地保留） |
| Li_GB_Project | 52 | 0 | 52（保持原位） | 0 |

> 说明：Li_Project 中“重复拷贝”按去重原则移入 temp/duplicates/，不删除任何文件。
