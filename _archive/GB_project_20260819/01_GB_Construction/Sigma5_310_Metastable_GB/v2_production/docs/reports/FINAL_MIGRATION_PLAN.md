# 最终迁移计划 FINAL_MIGRATION_PLAN.md（修订版 v2）

> 更新日期：2026-08-12（按用户确认的最终结构修订）
> 状态：已确认执行
> 目标根：`D:\Research\Lithium_GrainBoundary_Project`（不创建 Li_Project / Li_GB_Project 平级目录）

## 最终结构

```
D:\Research\Lithium_GrainBoundary_Project\
├── 00_Literature/
│   ├── grain_boundary/        ← GB 文献（3 篇）
│   ├── lithium_metal/         ← 锂金属文献（4 篇）
│   ├── ml_potential/          ← ML 势文献（空，占位）
│   ├── soap/                  ← SOAP 文献（2 篇）
│   └── md_methodology/        ← MD 方法文献（空，占位）
├── 01_GB_Construction/
│   └── Sigma5_310_Metastable_GB/
│       └── v2_production/     ← 当前 Sigma5(310)[001] 亚稳晶界项目（生产）
│           ├── src/generator/   （完整 Python 包，10 个文件）
│           ├── lammps/          （stage0_test…stage4_relax + 工具 + run_all.sh/resume.sh）
│           ├── data/            （initial/minimized/annealed/quenched/final/bulk/validation）
│           ├── tests/           （14 个测试 + pytest.ini）
│           ├── logs/            （16 个日志）
│           ├── docs/            （README.md + reports/ 含 5 份整理文档）
│           └── requirements.txt
├── 02_MD_Simulation/         （空，占位）
├── 03_Analysis/              （空，占位）
├── 04_Potential_Validation/  （空，占位）
├── 05_DFT_Validation/        （空，占位）
├── 06_Publication/           （空，占位）
├── 07_Documentation/
│   └── presentations/        ← Li_GB_Project/docs/（开题 PPTX + PDF，2 个）
├── archive/
│   ├── old_versions/
│   │   ├── v2_obsolete_20260802_150439/   （41 个文件）
│   │   └── legacy_generator/              （Li_GB_Project/Li_Bicrystal_Generator 完整版，31 个文件）
│   ├── duplicate_candidates/
│   │   ├── legacy_generator_partial/      （legacy_results 根文件，9 个）
│   │   ├── old_lammps_workflow/           （legacy_results/lammps_workflow，17 个）
│   │   ├── old_outputs/legacy_results_copy/（4 个）
│   │   ├── old_outputs/li_project_copy/   （4 个）
│   │   ├── Li_Bicrystal_Generator/        （Li_Project 重复副本，27 个，不含 __pycache__）
│   │   └── early_env_test/                （Li_Project 根早期测试，6 个）
│   ├── legacy/li_gb_project/              （Li_GB_Project 根文档 5 + lammps_test 5，共 10 个）
│   ├── debugging/cs5_310_ref_result.json  （1 个）
│   └── not_migrated/NOTES.md              （记录 venv 等未迁移项）
└── README.md                （项目总览）
```

## 移动数量（预计 247 个文件 + 新建少量说明文件）

| 来源 | 数量 | 去向 |
|---|---|---|
| 工作区 v2_final | 81 | v2_production |
| 工作区整理文档（5 份） | 5 | v2_production/docs/reports/ |
| 工作区 v2_obsolete | 41 | archive/old_versions/ |
| 工作区 legacy_results | 30 | archive（去重后进 duplicate_candidates） |
| 工作区 work | 1 | archive/debugging/ |
| Li_Project（旧目录） | 37 | archive/duplicate_candidates/ |
| Li_GB_Project | 52 | 拆分：archive/legacy_generator、archive/legacy/li_gb_project、00_Literature、07_Documentation |

## 不迁移列表（archive/not_migrated/NOTES.md 记录）

- Li_Project/venv（9215 文件，337 MB）：environment can be recreated from requirements.txt
- Li_Project/markdown file/test 1.md（无关教程）
- __pycache__ / *.pyc（Li_Project 生成器内，留原处，候选手动清理）
- 空目录：工作区 outputs/、v2_obsolete 的 logs/+outputs/、Li_Project/Li_Bicrysmkdir/、Li_GB_Project/notes/ 等

## 会修改的文件（仅 2 个，只改路径字符串）

- v2_production/lammps/run_all.sh：OUTDIR 默认 outputs→data；5 个 stage 脚本路径加 stage 子目录
- v2_production/lammps/resume.sh：OUTDIR 默认 outputs→data；4 个 SCRIPT 路径加 stage 子目录
- 禁止修改 LAMMPS 参数 / timestep / temperature / fix / pair_style / thermo / .in 内容 / Python 文件

## 执行方式

1. 创建 D:\Research\Lithium_GrainBoundary_Project 目录骨架。
2. 按上表逐项 Move-Item（源位置留空目录，不删除文件）。
3. 修改 run_all.sh / resume.sh 路径字符串。
4. 生成 archive/not_migrated/NOTES.md、根 README.md、占位 README。
5. 生成 PROJECT_MIGRATION_REPORT.md（放 07_Documentation/）。
6. 迁移后校验目标文件数与来源清空情况。
