# 项目整理报告 PROJECT_CLEANUP_REPORT.md

> 生成日期：2026-08-12
> 状态：**扫描与报告完成，等待确认；尚未移动/删除任何文件。**

## 1. 整理前文件数量

| 位置 | 文件数 |
|---|---|
| 工作区 v2_final | 81 |
| 工作区 legacy_results | 30 |
| 工作区 v2_obsolete_20260802_150439 | 41 |
| 工作区 work | 1 |
| 工作区 outputs（根，空） | 0 |
| Li_Project（不含 venv） | 55 |
| Li_Project/venv | 9215 |
| Li_GB_Project | 52 |
| **合计（不含 venv）** | **260** |

## 2. 目标结构（待确认根目录：C:\Users\1\Documents\Li_Project\）

```
Li_Project/
├── production/                 ← v2_final（81 个文件）
│   ├── generator/  lammps/  tests/  outputs/  logs/  validation/
│   └── README.md  FINAL_STATUS.md  ENVIRONMENT.txt  SHA256SUMS.txt  pytest.ini  requirements.txt
├── archive/
│   ├── legacy_generator/       ← 旧生成器 1 份（去重后）
│   ├── old_lammps/             ← legacy_results/lammps_workflow + 早期环境验证 1 份
│   ├── old_outputs/            ← 旧 3224 原子结构 1 份
│   └── debugging_history/      ← v2_obsolete(41) + work/cs5_310_ref_result.json
└── temp/
    ├── duplicates/             ← 重复拷贝（不删除）
    ├── unrelated/              ← markdown file/test 1.md
    └── venv/                   ← 建议原地保留，不移动
```

## 3. 保留的生产文件（FINAL → production/）

- `v2_final/generator/`（10 个）：config、geometry、builder、exporters、verify、lammps_io、generate、bulk、reimport_check、__init__
- `v2_final/lammps/`（11 个）：in.0_check、in.1_minimize、in.2_anneal、in.2_anneal_resume、in.3_quench、in.3_quench_resume、in.4_final、in.bulk_energy、in.dump2data、run_all.sh、resume.sh
- `v2_final/tests/`（14 个）：conftest + 13 个测试
- `v2_final/outputs/`（24 个）：初始结构 + Stage1–4 输出 + bulk 参考 + 报告
- `v2_final/logs/`（16 个）：各阶段日志与 .done
- `v2_final/validation/rebuilt_stage2.data`
- 根文件 6 个：README.md、FINAL_STATUS.md、ENVIRONMENT.txt、SHA256SUMS.txt、pytest.ini、requirements.txt

## 4. 归档文件（ARCHIVE → archive/）

| 归档目录 | 内容 | 文件数 |
|---|---|---|
| archive/legacy_generator/ | 旧生成器源码 1 份（建议 Li_GB_Project/Li_Bicrystal_Generator 完整版） | 约 17 |
| archive/old_lammps/ | legacy_results/lammps_workflow（17）+ 早期环境验证 1 份 | 约 23 |
| archive/old_outputs/ | 旧 Σ5_310 结构 1 份 | 4 |
| archive/debugging_history/ | v2_obsolete（41）+ cs5_310_ref_result.json（1） | 42 |

## 5. 建议删除 / 归入 temp 的文件（TEMP，均不直接删除）

| 内容 | 处理 |
|---|---|
| v2_obsolete（DO_NOT_USE） | 入 archive/debugging_history（不删） |
| 重复的旧生成器/旧输出/早期测试（2 份多余拷贝） | 移入 temp/duplicates/（不删） |
| markdown file/test 1.md | 移入 temp/unrelated/（不删） |
| Li_Bicrysmkdir/（空目录，拼写错误） | 空目录；是否删除目录待确认 |
| venv/（9215 文件） | 原地保留不移动（含绝对路径）；清理需手动 |
| __pycache__ / *.pyc | 自动缓存，可归 temp/cache/（不删） |
| 工作区空 outputs/、Li_GB_Project/notes/（空） | 空目录，无文件可移动 |

## 6. 关键文件检查结果（用户要求清单）

| 期望文件 | 是否存在于 production/lammps/（对应文件） | 说明 |
|---|---|---|
| stage0_test.in | 否（对应 in.0_check） | 命名不同，见 C4 |
| stage1_minimize.in | 否（对应 in.1_minimize） | 同上 |
| stage2_anneal.in | 否（对应 in.2_anneal） | 同上 |
| stage3_quench.in | 否（对应 in.3_quench） | 同上 |
| stage4_relax.in | 否（对应 in.4_final） | Stage 4 实际是 0 K 最小化，非 relax |

## 7. 发现的问题清单

1. 旧生成器在 3 个位置重复（legacy_results / Li_Project / Li_GB_Project），需去重。
2. 旧流程存在明确错误（3224 原子、10.44% 删除、三斜盒、全局删原子、`nvt 300.0 0.0`、`in.resume` 复制主脚本、dump2data 空盒）。
3. 用户期望的 stage*.in 文件名与 v2_final 实际文件名不一致（建议保留实际名并文档化映射）。
4. `C:\Users\1\Documents\Li_Project` 根目录混乱（venv、教程 md、空目录、早期测试），且与 Li_GB_Project 高度重复。
5. v2_obsolete 的松弛策略（x 固定/y,z 松弛）与 v2_final（x 松弛/y,z 固定）相反——以 v2_final 为准。
6. 势函数无冲突；晶格常数差异（3.51 vs 3.506186）为有意更新。
7. v2_final 的 README 明确声明最终结构是「单条轨迹的局部极小候选」，不是统计代表性系综——引用结果时需注意。

## 8. 待确认事项（确认后再执行移动）

- [ ] 目标根目录：`C:\Users\1\Documents\Li_Project\` 是否 OK？（或改为在工作区内创建 Li_Project/）
- [ ] 是否同意「保留 1 份归档 + 其余重复拷贝移入 temp/duplicates/」的去重方式（不删除）
- [ ] 是否同意 v2_obsolete 入 archive/debugging_history、venv 原地保留
- [ ] 是否同意保留 v2_final 的 in.* 命名并只文档化映射（不新建 stage*.in 别名）
