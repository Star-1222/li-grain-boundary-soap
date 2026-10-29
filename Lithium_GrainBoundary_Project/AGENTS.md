# 锂晶界科研项目 — 项目级上下文

本文件只在本目录（或其子目录）工作时由 agent 加载。全局环境与偏好见 `~/.dsh/AGENTS.md`。

## 研究主题

**锂金属亚稳晶界（Li metastable grain boundary）的原子结构与力学响应**。

- 核心研究链条：**SOAP 结构描述符 → 亚稳程度 M → 不同受力机制 → 临界剪切强度 τ_c → M–τ_c 映射**。方法组合为 LAMMPS + SOAP + 剪切加载。
- 相关文档：`Li_metal_metastable_atomic_structure_mechanical_analysis.pdf`（14 页，2026 年大创立项报告）。
- 模拟数据规模：`natoms=14400`（Li Σ5 晶界模型），6 个 MPI 进程，5 ps 时长，时间步 0.001。
- 控压（barostat）方案对照实验（2026-10-15，归档在仓库根 `Li_Sigma5_控压方案对照_2026-10-15/`）：`xonly` 与 `aniso` 两个分支，`fix npt temp 10 10 0.1 ${scheme} 0.0 0.0 ${damp}`；日志在 `03_原始日志/`（`baro_xonly_fast.lammps`、`baro_aniso_fast.lammps`）。结论：偏应力 627.6 → 15.9 bar，代价是 ε_zz = −0.99%。
- 结构分析用配位数（CN）判据：第一近邻截断 3.25 Å，取 g(r) 第一谷。

## 方法论争议点（讨论时主动指出）

「结构无序程度」不等于「热力学/动力学意义上的亚稳程度」。一个晶界可以在结构上很无序但处于稳定非晶态，也可能结构看起来有序却因制备动力学锁在高于平衡态的构型上。若 M 只由 SOAP 距离或无序度导出，则所有 M–τ_c 映射都建立在一个语义不牢的变量上。

## 已确认的技术坑（不要重复踩）

1. **`v_fBCC` 恒为 0**：`compute cna/atom 3.27` 配合 `c_cCNA == 3` 在 MD 快照上不工作，26 行采样全部输出 0.000，而同一结构用配位数判据是 80.6%。这个读数不能作为证据引用，要么修要么从 thermo 里删。
2. **两套第一近邻截断不一致**：CNA 用 3.27 Å，后处理统一用 3.25 Å。同一研究里应统一口径。
3. **采样点数误算**：`thermo 200` 在 5000 步内产生 26 行，后处理取 frac=0.5 即最后 13 行（窗口 2600–5000 步 = 2.4 ps）。旧报告把「全程 26 行」当成了「窗口内 26 个采样点」。按 13 点算，偏应力均值的标准误约 1.7 bar（未做自相关校正），不是 6.3 bar。
4. **dump 太稀**：`dump ... 2000` 在 5 ps 里只有 3 帧，无法做 εzz 的局域分布分析。下轮应改为 200–500 步，并在中途加 `restart` checkpoint。
5. **Pdamp 单位陷阱**：旧协议 `Pdamp = 1000 ps`，新协议 `0.1 ps`。日志的 `BRANCH` 行会回显实际值，务必核对。
6. **脚本注释与实际运行不符**：`in.baro_hold`、`run_barostat.sh` 的注释写着 20000 步 / 12 进程 / `Pdamp 1000`，实际跑的是 5000 步 / 6 进程 / `*_fast`。归档前要把注释改成一致。

## 目录

- `00_Literature/` — 文献，按 `grain_boundary/`、`lithium_metal/`、`soap/` 分主题
- `01_GB_Construction/Sigma5_310_Metastable_GB/` — 晶界构建
- `Li_Sigma5_310_metastable_2x2/` — 主模拟：`PARAMETERS.json`、`generator/`、`lammps/`、`scripts/`、`stage0_build_validate/`
- `Reports/` — `GB_Construction/`、`Li_2x2_Workflow/`、`Team_Reviews/`、`README_BI.md`
- `07_Documentation/presentations/` — 汇报材料
