# 版本冲突报告 conflict_report.md

> 生成日期：2026-08-12
> 规则：发现冲突不自动覆盖；逐条记录「文件A / 文件B / 区别 / 推荐保留 / 理由」。

## C1. 多个 generator（5 份）

| 项 | 值 |
|---|---|
| 文件A | `v2_final/generator/`（10 个 .py） |
| 文件B | `legacy_results/`（builder/config/geometry/exporter/main.py） |
| 文件C | `Li_Project/Li_Bicrystal_Generator/`（完整旧包） |
| 文件D | `Li_GB_Project/Li_Bicrystal_Generator/`（与 C 完全重复） |
| 文件E | `v2_obsolete_20260802_150439/generator/`（v2 中间版） |
| 区别 | A=新架构（镜像法、3600 原子、已验证）；B/C/D=旧生成器（primitive 当 conventional、固定盒旋转、全局删原子、3224 原子），三份 SHA256 完全一致；E=几何已修正但工作流未定型 |
| 推荐保留 | A（production/generator/）；B/C/D 去重保留 1 份（archive/legacy_generator/）；E 入 archive/debugging_history/ |
| 理由 | A 是唯一通过 63 项测试并产出 v2_final 结果的版本；B/C/D 为同一历史版本，保留一份即可；E 作为历史中间态 |

## C2. 多个 config.py

| 项 | 值 |
|---|---|
| 文件A | 旧 `config.py`（legacy_results/Li_Project/Li_GB_Project 三份相同）：LATTICE_CONSTANT=3.51、SUPERCELL=(20,12,12)、OVERLAP_CUTOFF=2.3、Σ5/13/17/29 GB 数据库 |
| 文件B | `v2_final/generator/config.py`：DEFAULT_A=3.51（构建初值）、KIM_MODEL=SNAP_ZuoChenLi_2019_Li__MO_732106099012_000 |
| 文件C | `v2_obsolete_.../generator/config.py`：与 B 基本一致（旧策略） |
| 区别 | 旧版含全局删除参数（错误根源）；新版为镜像法参数；v2_final 实际用 a=3.506186475312 生成初始结构 |
| 推荐保留 | B（production/generator/config.py） |
| 理由 | 新参数对应已验证的 3600 原子正交盒；旧参数是历史错误来源，仅归档 |

## C3. 多个 in.* LAMMPS 输入

| 项 | 值 |
|---|---|
| 文件A | 旧 `lammps_workflow/in.*`（in.bicrystal_relax、in.quench、in.resume、anneal_save.in、continue.in、dump2data.in） |
| 文件B | `v2_final/lammps/in.*`（in.0_check … in.4_final + resume 脚本） |
| 文件C | `v2_obsolete_.../lammps/in.*`（in.0_check…in.4_final + in.resume） |
| 区别 | 旧：`fix nvt temp 300.0 0.0`、`in.resume` 与主脚本相同、`dump2data.in` 空盒 add no、minimize 1e-6/1e-8；新：淬火目标 1 K、真 resume、box yes add yes、minimize 0.0 1e-6 |
| 推荐保留 | B（production/lammps/） |
| 理由 | B 是实际跑完 Stage 0–4 并全部 `.done` 的版本 |

## C4. 用户要求的关键文件名 vs 实际文件名（命名冲突，不自动改名）

| 项 | 值 |
|---|---|
| 期望文件名 | `stage0_test.in`、`stage1_minimize.in`、`stage2_anneal.in`、`stage3_quench.in`、`stage4_relax.in` |
| 实际文件名 | `in.0_check`、`in.1_minimize`、`in.2_anneal`、`in.3_quench`、`in.4_final` |
| 区别 | 命名风格不同；Stage 4 实际是「0 K 最小化」而非「relax」 |
| 推荐保留 | 保留实际已验证文件名，并在 production/README.md 中建立映射表 |
| 理由 | 改名会破坏 run_all.sh、resume 脚本与文档引用；若您要求严格按期望名，可另建别名副本 |

## C5. 晶格常数差异（有意为之）

| 项 | 值 |
|---|---|
| 文件A | 旧结构：a=3.51 Å（初值） |
| 文件B | v2_final：a=3.506186475312 Å（SNAP 平衡晶格常数） |
| 区别 | 旧流程未做势函数平衡晶格标定；v2_final 以平衡 a 构建初始结构 |
| 推荐保留 | v2_final 的 a=3.506186475312 |
| 理由 | 与 SNAP 势函数平衡态一致（params.json、README 均有记录） |

## C6. 原子数量差异（反映几何错误 vs 修复）

| 项 | 值 |
|---|---|
| 文件A | 旧：3600→3224（删 10.44%），盒 85.977×60.795×14.892 三斜 |
| 文件B | 新：3720→3600（删 3.2258%，仅 GB 面重合原子），盒 66.525×33.263×35.062 正交 |
| 区别 | 删除范围、盒型、晶粒原子数均不同 |
| 推荐保留 | 新（3600） |
| 理由 | 旧模型被回归测试明确判为错误场景（legacy 3600→3224 不得复现） |

## C7. 势函数：无冲突

| 项 | 值 |
|---|---|
| 检查项 | 旧 in.bicrystal_relax、v2_final、v2_obsolete 全部为 `SNAP_ZuoChenLi_2019_Li__MO_732106099012_000` |
| 结论 | 一致，无冲突，无需处理 |

## C8. 晶胞松弛策略冲突（v2 中间版 vs v2_final）

| 项 | 值 |
|---|---|
| 文件A | `v2_obsolete/lammps/in.1_minimize`：x 固定，y/z 松弛 |
| 文件B | `v2_final/lammps/in.1_minimize`、in.4_final：x 松弛（零法向压力），y/z 固定（保持面内 CSL 周期） |
| 区别 | 方向策略相反，均为文档化决定 |
| 推荐保留 | v2_final 策略 |
| 理由 | 实际跑完 Stage 0–4，Pxx≈3.6e-4 bar 收敛；v2_final README/FINAL_STATUS 明确记录 |

## C9. 重复文件（去重方案，不删除）

- 旧生成器 5 个核心 .py：3 份（legacy_results、Li_Project、Li_GB_Project）SHA256 完全一致 → 保留 1 份入 archive/legacy_generator/，其余移入 temp/duplicates/
- 旧输出 Σ5_310.cif/.lmp/.xyz/_grainID.xyz：2 份（legacy_results/output、Li_Project/output）→ 同上去重
- 早期 LAMMPS 环境验证：2 份（Li_Project 根、Li_GB_Project/lammps_test）→ 同上去重
- legacy_results 是 Li_GB_Project 的子集（缺 utils.py/__init__.py/__main__.py/examples/）→ 以 Li_GB_Project 版本为完整归档源
