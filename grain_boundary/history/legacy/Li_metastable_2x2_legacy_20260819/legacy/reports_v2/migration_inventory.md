# 迁移清单 migration_inventory.md

> 生成日期：2026-08-12
> 目标位置：`D:\Research\Li_Project`（D:\Research 目前不存在，确认后创建）
> 来源：工作区 `C:\Users\1\Documents\Codex\2026-08-01\bcc-li-5-001-sigma5-310`、`C:\Users\1\Documents\Li_Project`、`C:\Users\1\Documents\Li_GB_Project`
> 原则：不删除任何文件；只标记删除候选；不修改代码/LAMMPS 参数/Python 文件；只做整理与路径迁移。

---

## 0. 迁移映射总览

| 目标目录 | 来源 |
|---|---|
| src/generator/ | v2_final/generator（完整包，见冲突 C1） |
| lammps/stage0_test … stage4_relax/ | v2_final/lammps/in.*（见冲突 C2） |
| data/initial, minimized, annealed, quenched, final/ | v2_final/outputs 各阶段结构 |
| tests/ | v2_final/tests |
| logs/ | v2_final/logs |
| docs/README.md, docs/reports/ | v2_final/README.md + 验证/状态/校验文件 |
| archive/old_versions/ | v2_obsolete + 旧生成器 1 份 + 旧 LAMMPS 工作流 1 份 + 旧输出 1 份 |
| archive/debugging/ | work/cs5_310_ref_result.json + 早期环境测试 + 调试日志 |

---

## 1. v2_final → 生产目录（81 个文件，全部迁移）

### 1.1 → src/generator/（10 个）

| 原路径（v2_final/generator/） | 类型 | 大小 | 是否为空 | 迁移 | 目标 |
|---|---|---|---|---|---|
| __init__.py | py | 86 B | 否 | 是 | src/generator/__init__.py |
| config.py | py | 2123 B | 否 | 是 | src/generator/config.py |
| geometry.py | py | 9860 B | 否 | 是 | src/generator/geometry.py |
| builder.py | py | 8990 B | 否 | 是 | src/generator/builder.py |
| exporters.py | py | 5732 B | 否 | 是 | src/generator/exporters.py |
| verify.py | py | 17145 B | 否 | 是 | src/generator/verify.py |
| lammps_io.py | py | 3485 B | 否 | 是 | src/generator/lammps_io.py |
| generate.py | py | 10320 B | 否 | 是 | src/generator/generate.py |
| bulk.py | py | 2345 B | 否 | 是 | src/generator/bulk.py |
| reimport_check.py | py | 4313 B | 否 | 是 | src/generator/reimport_check.py |

> **冲突 C1**：目标结构要求的 src/geometry、src/exporters、src/verification 会拆散 Python 包（`from .geometry import …` 相对导入），违反“不修改 Python 文件”。建议保留完整包于 src/generator/；src/geometry、src/exporters、src/verification 可各放一个 README 说明模块位置（或按您确认后另议）。

### 1.2 → lammps/（11 个，拆分到 stage 子目录）

| 原路径（v2_final/lammps/） | 类型 | 大小 | 迁移 | 目标 |
|---|---|---|---|---|
| in.0_check | in | 662 B | 是 | lammps/stage0_test/in.0_check |
| in.1_minimize | in | 1549 B | 是 | lammps/stage1_minimize/in.1_minimize |
| in.2_anneal | in | 953 B | 是 | lammps/stage2_anneal/in.2_anneal |
| in.2_anneal_resume | in | 1095 B | 是 | lammps/stage2_anneal/in.2_anneal_resume |
| in.3_quench | in | 964 B | 是 | lammps/stage3_quench/in.3_quench |
| in.3_quench_resume | in | 1069 B | 是 | lammps/stage3_quench/in.3_quench_resume |
| in.4_final | in | 1682 B | 是 | lammps/stage4_relax/in.4_final |
| in.bulk_energy | in | 990 B | 是 | lammps/stage4_relax/in.bulk_energy（或 lammps/） |
| in.dump2data | in | 546 B | 是 | lammps/ |
| run_all.sh | sh | 1226 B | 是 | lammps/run_all.sh（需更新脚本路径引用，见 C2） |
| resume.sh | sh | 1943 B | 是 | lammps/resume.sh（同上） |

> **冲突 C2**：run_all.sh / resume.sh 用相对文件名引用 in.*。拆分到 stage 子目录后，这两个脚本需做“仅路径字符串”修改（如 `in.0_check` → `stage0_test/in.0_check`），不涉及任何 LAMMPS 参数。确认后才改。

### 1.3 → data/（24 个，按阶段归类）

| 原路径（v2_final/outputs/） | 大小 | 迁移 | 目标 |
|---|---|---|---|
| Sigma5_310.data / .xyz / .cif / .grain_map / _grain_id.xyz | 各数十万 B | 是 | data/initial/ |
| removal_log.csv | 8467 B | 是 | data/initial/ |
| params.json | 557 B | 是 | data/initial/ |
| build_report.md / verification_report.md | ~2 KB | 是 | docs/reports/ |
| stage1_minimized.data / .restart | 284566/317785 B | 是 | data/minimized/ |
| stage2_anneal.dump | 8283528 B | 是 | data/annealed/ |
| stage2_annealed.data / .restart | 467158/317785 B | 是 | data/annealed/ |
| stage3_quench.dump | 10373208 B | 是 | data/quenched/ |
| stage3_quenched.data / .restart | 284332/317785 B | 是 | data/quenched/ |
| stage4_final.data / .restart | 284352/317785 B | 是 | data/final/ |
| bulk/bulk_Li.cif/.data/.xyz + bulk_Li_minimized.data/.restart | 各 ~10 KB | 是 | data/bulk/（目标结构未列，新增建议） |

### 1.4 → tests/（14 个）

v2_final/tests/ 全部（conftest.py + 13 个 test_*.py）→ `tests/`。

### 1.5 → logs/（16 个）

v2_final/logs/ 全部（含 5 个 .done、各阶段 .lammps/.screen、kim.log、bulk_energy.lammps、dump2data_validation.lammps）→ `logs/`。

### 1.6 → docs/（根文件 6 个）

| 原路径 | 目标 |
|---|---|
| v2_final/README.md | docs/README.md |
| v2_final/FINAL_STATUS.md | docs/reports/FINAL_STATUS.md |
| v2_final/ENVIRONMENT.txt | docs/reports/ENVIRONMENT.txt |
| v2_final/SHA256SUMS.txt | docs/reports/SHA256SUMS.txt |
| v2_final/pytest.ini | tests/pytest.ini |
| v2_final/requirements.txt | requirements.txt（项目根，或 docs/） |

### 1.7 → data/ 或 docs/（validation 1 个）

v2_final/validation/rebuilt_stage2.data → `data/validation/rebuilt_stage2.data`（或 data/initial/）。

---

## 2. legacy_results → archive/（30 个文件，去重后迁移）

| 原路径 | 数量 | 去重 | 目标 |
|---|---|---|---|
| legacy_results/*.py（builder/config/geometry/exporter/main）+ README_generator/README_project/requirements/SETUP | 9 | 与 Li_Project、Li_GB_Project 完全重复 → 保留 1 份 | archive/old_versions/legacy_generator/ |
| legacy_results/lammps_workflow/*（17 个：in.*、data、restart、dump、log） | 17 | 与 Li_GB_Project 完全重复 → 保留 1 份 | archive/old_versions/old_lammps/ |
| legacy_results/output/*（4 个） | 4 | 与 Li_Project/output 完全重复 → 保留 1 份 | archive/old_versions/old_outputs/ |

---

## 3. v2_obsolete_20260802_150439 → archive/old_versions/（41 个文件，全部迁移）

DO_NOT_USE.txt、README.md、pytest.ini、requirements.txt、generator/(10)、lammps/(10)、reports/(4)、tests/(13) → `archive/old_versions/v2_obsolete_20260802_150439/`（整目录）。
> 含空目录 logs/、outputs/（空目录不迁移，见第 7 节）。

---

## 4. work → archive/debugging/（1 个文件）

cs5_310_ref_result.json（34 KB）→ `archive/debugging/cs5_310_ref_result.json`。

---

## 5. C:\Users\1\Documents\Li_Project（55 个文件，不含 venv）

| 路径 | 类型 | 大小 | 是否为空 | 迁移 | 目标 |
|---|---|---|---|---|---|
| Li_Bicrystal_Generator/（源码 10 + examples + lammps_workflow 17 + output 空 + __pycache__ 17） | py/in/data/sh/log/dump | 约 10 MB | 否 | 去重：已由 legacy_results 覆盖（缺 utils/__init__/__main__/examples → 以 Li_GB_Project 版本为完整源） | archive/old_versions/legacy_generator/（完整版 1 份） |
| output/Σ5_310.*（4 个） | 结构文件 | 各 ~150-287 KB | 否 | 重复 → 候选删除 | — |
| README.md / run.sh / final.data / in.test_li / kim.log / log.lammps | 早期环境测试 | 各 0.5-9 KB | 否 | 与 lammps_test 重复 → 保留 1 份 | archive/debugging/early_env_test/ |
| markdown file/test 1.md | md 教程 | 2448 B | 否 | 无关文件 → 候选删除 | — |
| venv/ | Python 环境 | 9215 文件，337 MB | 否 | 不迁移（绝对路径失效），候选删除（手动） | — |
| Li_Bicrysmkdir/ | 空目录（拼写错误） | — | 是 | 不迁移，候选删除 | — |
| Li_Bicrystal_Generator/output/ | 空目录 | — | 是 | 不迁移，候选删除 | — |

---

## 6. C:\Users\1\Documents\Li_GB_Project（52 个文件）

| 路径 | 建议 |
|---|---|
| README.md、SETUP.md、requirements.txt、setup_env.sh、.gitignore | 可选：迁移到 archive/old_versions/legacy_generator/（与 legacy_results 同内容） |
| docs/（PDF+PPTX 开题材料，2 个） | 可选：迁移到 docs/presentations/（目标结构未列，需确认） |
| literature/（8 篇 PDF + soap_method/） | 可选：迁移到 docs/literature/（目标结构未列，需确认） |
| lammps_test/（6 个） | 可选：迁移到 archive/debugging/early_env_test/（与 Li_Project 重复） |
| Li_Bicrystal_Generator/ | 归档完整源（含 utils/__init__/__main__/examples） |
| notes/（空目录） | 不迁移，候选删除（空目录） |

> 默认建议：**Li_GB_Project 保持原位**（已组织良好），仅当其内容需并入 D:\ 时才迁移；见待确认项。

---

## 7. 空目录清单（不迁移，标记候选删除）

| 空目录 | 说明 |
|---|---|
| 工作区 outputs/ | 空 |
| v2_obsolete/logs/、v2_obsolete/outputs/ | 空 |
| Li_Project/Li_Bicrysmkdir/ | 拼写错误目录 |
| Li_Project/Li_Bicrystal_Generator/output/ | 空 |
| Li_GB_Project/notes/ | 空 |

---

## 8. 删除候选清单（只标记，不删除）

| 项目 | 数量/大小 | 理由 |
|---|---|---|
| 重复的旧生成器拷贝 | 2 份 | 保留 1 份归档 |
| 重复的旧输出（legacy_results/output vs Li_Project/output） | 1 份 | 保留 1 份归档 |
| 重复的旧 lammps_workflow | 1 份 | 保留 1 份归档 |
| 重复的早期环境测试 | 1 份 | 保留 1 份归档 |
| __pycache__ / *.pyc | 17 个 | Python 缓存 |
| venv/ | 9215 文件 / 337 MB | 环境产物，可再生 |
| markdown file/test 1.md | 1 个 | 无关教程 |
| 空目录 | 5 个 | 见第 7 节 |
| v2_obsolete（旧版） | 41 个 | 已确认无科研价值，仅归档（不删除，移入 archive） |

---

## 9. 特殊检查结果

### 9.1 重复版本（已核验，SHA256）

- 旧生成器核心 5 .py：3 份相同（legacy_results / Li_Project / Li_GB_Project）。
- 旧输出 4 文件：legacy_results/output == Li_Project/output。
- 旧 lammps_workflow 数据：legacy_results == Li_GB_Project。
- 早期环境测试 5 文件：Li_Project 根 == Li_GB_Project/lammps_test。
- legacy_results/README_project.md == Li_GB_Project/README.md。

### 9.2 硬编码路径 / 旧目录引用（不改内容，仅记录）

| 文件 | 内容 | 处理 |
|---|---|---|
| v2_obsolete/README.md、reports/reference_energy_bulk.md | `/mnt/c/Users/1/Documents/Codex/2026-08-01/bcc-li-5-001-sigma5-310/v2` | 归档，保留原文（含过期路径） |
| Li_GB_Project/SETUP.md、setup_env.sh | `cd ~/Li_GB_Project`、`PROJECT_DIR="$HOME/Li_GB_Project"` | 归档/原位，保留原文 |
| v2_final 全部 .in/.sh | 无绝对路径（${outdir}/${input_data}/${input_restart}/$HERE） | 干净，无动作 |

### 9.3 LAMMPS 脚本迁移后注意

- in.1_minimize 用 `${input_data}`、in.4_final 用 `${input_restart}`、其余用 `${outdir}`：迁移后用法不变，只需在命令行把 outdir 指向 `D:\Research\Li_Project\data` 等。
- run_all.sh 默认 `OUTDIR=$HERE/../outputs`：确认后改为 `../data`，且脚本路径改到 stage 子目录（见 C2）。

---

## 10. 待确认事项（确认后再执行移动）

- [ ] 1. 确认目标根 `D:\Research\Li_Project`（将创建 D:\Research）。
- [ ] 2. src 结构：接受“src/generator 保持完整包、geometry/exporters/verification 仅作模块说明”的方案吗？（拆分会破坏 Python 导入）
- [ ] 3. lammps 按 stage 子目录拆分，并仅修改 run_all.sh / resume.sh 的脚本路径引用，可以吗？
- [ ] 4. 重复文件“保留 1 份入 archive、其余标记删除候选（不删除）”可以吗？
- [ ] 5. venv（337 MB）不迁移、标记候选删除，可以吗？
- [ ] 6. Li_GB_Project 的 docs/literature 是否并入 D:\Research\Li_Project\docs/？
