# PROJECT_MIGRATION_REPORT.md

> 生成：2026-08-12（由 MIGRATION_RESUME.ps1 完成未迁移部分后生成）

## 1. 原项目位置
- C:\Users\1\Documents\Codex\2026-08-01\bcc-li-5-001-sigma5-310（v2_final / legacy_results / v2_obsolete / work）
- C:\Users\1\Documents\Li_Project（旧目录：生成器重复副本、早期环境测试、venv 等）
- C:\Users\1\Documents\Li_GB_Project（文献、演示、旧生成器、lammps_test）

## 2. 新项目位置
- D:\Research\Lithium_GrainBoundary_Project（根）
- 生产包：01_GB_Construction\Sigma5_310_Metastable_GB\v2_production
- 文献：00_Literature；文档：07_Documentation；归档：archive

## 3. 迁移文件数量
- 目标盘文件总数：257
- 说明：v2_final 81、legacy_results 30、v2_obsolete 41、work 1、整理文档 5、Li_Project 旧目录 37、Li_GB_Project 52

## 4. 空目录（保留在原处，未迁移）
- 工作区 outputs/；v2_obsolete 的 logs/ 与 outputs/；Li_Project\Li_Bicrysmkdir/；Li_GB_Project\notes/

## 5. 删除候选（只标记，未删除）
- Li_Project\venv（9215 文件，约 337 MB；可由 requirements.txt 重建）
- Li_Project\markdown file\test 1.md（无关教程）
- __pycache__ / *.pyc（Li_Project 生成器内，已留原处）
- 重复文件：已移入 archive\duplicate_candidates\（不删除）

## 6. 版本冲突
- 旧生成器 3 份相同：保留 Li_GB 完整版于 archive\old_versions\legacy_generator；其余进 duplicate_candidates
- 旧输出/旧工作流/早期测试：各保留 1 份，其余进 duplicate_candidates
- 势函数三处一致（SNAP_ZuoChenLi_2019_Li__MO_732106099012_000），无冲突
- 晶格常数 3.51 vs 3.506186 Å、原子数 3224 vs 3600：为几何修复后的有意差异
- run_all.sh / resume.sh 仅修改路径字符串；LAMMPS 参数 / .in 文件 / Python 文件未改动

## 7. 未迁移/缺失记录
- 见 archive\not_migrated\NOTES.md
- 缺失文件记录：migration_missing_files.txt（本次运行共有 0 条）

## 8. 后续推荐工作目录
- D:\Research\Lithium_GrainBoundary_Project\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production
- 复现命令参考 v2_production\docs\README.md