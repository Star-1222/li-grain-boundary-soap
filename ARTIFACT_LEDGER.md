# 研究产物台账

> 本文件由 `artifact_ledger.py` 自动生成，请勿手工编辑。
> 每次重新生成都会刷新文件清单、校验和与内容分析。

- 文件总数：**521**
- 总体积：**54.40 MB**

## 层次分布

| 层次 | 文件数 | 体积 | 占总比 | 处置 |
|---|---:|---:|---:|---|
| **结论** | 192 | 35.05 MB | 64.4% | 必须保留 |
| **续算** | 11 | 6.60 MB | 12.1% | 保留至后续阶段完成 |
| **追溯** | 39 | 1.06 MB | 2.0% | 随代码进版本控制 |
| **代码** | 134 | 902.9 KB | 1.6% | 进版本控制 |
| **输入** | 51 | 8.10 MB | 14.9% | 进版本控制 |
| **文档** | 94 | 2.71 MB | 5.0% | 进版本控制 |

## 自动识别的结构特征

- 结论性产物只有 35.05 MB，占 64.4%。
- 决定“可复现”的追溯+代码+输入合计仅 10.04 MB（18.5%），是最该保护的一层。
- 双语拼接文档 16 份。中英文同处一文件，改动一侧另一侧会静默过期，建议拆分为 `README.md` + `README.en.md`。

## 交叉核对：报告声称 vs 数据实测

- **能量模型自洽**：`stage6_energy_models.csv` 实测 Model0_D: R²=0.0393 · Model1_D_phys: R²=0.9869 · Model2_D_phys_dir: R²=0.9886 · ModelD_D_dir: R²=0.7497。与最终报告声称的 D-only 0.0393、D+物理 0.9869、D+物理+方向 0.9886、D+方向 0.7497 一致。
- **验收口径风险已量化**：`p1_hold_pressure.stats` 显示 2 个采样点 · 标量时均压力 -0.87 bar · 分量 Pxx=-187, Pyy=+383, Pzz=-199 · ⚠ 分量极差 582 bar，各向异性显著。报告把“0 bar”定义为标量时均压力，若改为要求应力张量各分量趋零，需补做 aniso NPT 敏感性复核。
- **无序化进程可量化**：`cna_phase1_report.txt` 给出 晶界 BCC: initial 73.3% → p1_min0 60.0% → casecheck 43.7% → p1_final 38.1% · 晶内 BCC: 100.0% → 100.0% → 99.9% → 98.7% · 无序带宽: 2.5 Å → 4.5 Å → 4.5 Å → 6.5 Å。晶界 BCC 占比持续下降而晶内基本不变，这正是“晶界真实无序、晶粒未熔化”的证据链。
- **重大发现（需人工确认）**：`stage6_correlations.csv` 显示 无序度 D 与位移范数相关 0.99998 · D 与过剩能相关 0.277。若 D 与位移范数几乎完全重合，则 SOAP 派生坐标相对简单几何量并未提供额外信息——这会动摇整个方法线的前提，建议优先核实。

## 关键产物内容分析

以下是直接从文件内容读出的事实，不是人工描述。

| 文件 | 内容分析 |
|---|---|
| `ni_stable.cfg` | 47,120 原子 |
| `ni_metastable.cfg` | 47,120 原子 |
| `p1_held.data` | 14,400 原子 · 盒子 66.17 × 66.10 × 69.66 Å · timestep 5,000 |
| `p1_held_0K.data` | 14,400 原子 · 盒子 66.17 × 66.10 × 69.66 Å · timestep 187 |
| `p1_quenched.data` | 14,400 原子 · 盒子 66.17 × 66.11 × 69.66 Å · timestep 35,300 |
| `p1_annealed.data` | 14,400 原子 · 盒子 66.48 × 66.41 × 69.98 Å · timestep 25,000 |
| `p1_casecheck_0K.data` | 14,400 原子 · 盒子 66.48 × 66.41 × 69.98 Å · timestep 934 |
| `stage5a_matched_neighborhoods.xyz` | 63 原子 · 附加列: central=290 state=stable D10=9.768196778712806 dE=0.3807500000000008 Eres=0.003080000000000638 CN=10.0 RMS=0.052540184592867356 CSP=27.686566695884597 dV=-0.11575476290275333 |
| `stage5a_atomwise_physical_descriptors.csv` | 4,659 行 × 23 列 [state, atom_id, x, y, z, interface_id …] |
| `stage5b_atomwise_soap_direction.csv` | 4,659 行 × 25 列 [structure, atom_id, x, y, z, interface_id …] |
| `Sigma5_310_K6_J6_L20.data` | 14,400 原子 · 盒子 66.53 × 66.53 × 70.12 Å |
| `structure.data` | 14,400 原子 · 盒子 66.53 × 66.53 × 70.12 Å |
| `stage4b_atomwise_energy_disorder.csv` | 4,659 行 × 21 列 [structure, region, atom_id, x, y, z …] |
| `stage2_75_atom_metadata.csv` | 4,659 行 × 19 列 [structure, region, atom_id, x, y, z …] |
| `p1_min0.data` | 14,400 原子 · 盒子 66.24 × 66.18 × 69.73 Å · timestep 12,500 |
| `stage4a_atomwise_disorder.csv` | 4,659 行 × 17 列 [structure, region, atom_id, x, y, z …] |
| `Sigma5_310_K6_J6_L20.xyz` | 14,400 原子 · 附加列: Lattice="66.5252109796 0.0 0.0 0.0 66.5252109796 0.0 0.0 0.0 70.1237295062" Properties=species:S:1:pos:R:3:grain_id:I:1 pbc="T T T" |
| `Sigma5_310_K6_J6_L20.xyz` | 14,400 原子 · 附加列: Lattice="66.5252109796 0.0 0.0 0.0 66.5252109796 0.0 0.0 0.0 70.1237295062" Properties=species:S:1:pos:R:3:grain_id:I:1 pbc="T T T" |
| `structure.xyz` | 14,400 原子 · 附加列: Lattice="66.5252109796 0.0 0.0 0.0 66.5252109796 0.0 0.0 0.0 70.1237295062" Properties=species:S:1:pos:R:3:grain_id:I:1 pbc="T T T" |
| `Sigma5_310_K6_J6_L20_grain_id.xyz` | 14,400 原子 · 附加列: Li x y z grain_id |
| `Sigma5_310_K6_J6_L20_grain_id.xyz` | 14,400 原子 · 附加列: Li x y z grain_id |
| `stage3b_atomwise_mid_epsilon.csv` | 4,659 行 × 14 列 [row_index, atom_id, x, y, z, class_id …] |
| `stage3b_mid_epsilon_lae.xyz` | 4,659 原子 · 附加列: Ni LAE_LER representation for epsilon_mid (visual inspection only) |
| `stage2_annealed.data` | 3,600 原子 · 盒子 65.95 × 33.26 × 35.06 Å · timestep 20,000 |
| `final_gb.data` | 3,224 原子 · 盒子 85.98 × 61.20 × 15.30 Å · timestep 50,044 |
| `annealed.data` | 3,224 原子 · 盒子 85.98 × 61.20 × 15.30 Å · timestep 20,277 |
| `stage3a_nearest_neighbor_distances.csv` | 4,659 行 × 9 列 [row_index, structure, region, atom_id, nearest_neighbor_row, nearest_neighbor_structure …] |
| `Sigma5_310.data` | 3,600 原子 · 盒子 66.53 × 33.26 × 35.06 Å |
| `stage1_minimized.data` | 3,600 原子 · 盒子 65.95 × 33.26 × 35.06 Å · timestep 5,000 |
| `stage4_final.data` | 3,600 原子 · 盒子 66.01 × 33.26 × 35.06 Å · timestep 3,000 |
| `stage3_quenched.data` | 3,600 原子 · 盒子 65.95 × 33.26 × 35.06 Å · timestep 50,000 |
| `Sigma5_310.lmp` | 3,224 原子 · 盒子 85.98 × 60.79 × 14.89 Å |
| `Sigma5_310.lmp` | 3,224 原子 · 盒子 85.98 × 60.79 × 14.89 Å |
| `Sigma5_310_minimized.data` | 3,224 原子 · 盒子 85.98 × 61.20 × 15.30 Å · timestep 277 |
| `rebuilt_stage2.data` | 3,600 原子 · 盒子 65.95 × 33.26 × 35.06 Å · timestep 20,000 |
| `Sigma5_310.xyz` | 3,224 原子 |
| `Sigma5_310.xyz` | 3,600 原子 · 附加列: Lattice="66.5252109796 0.0 0.0 0.0 33.2626054898 0.0 0.0 0.0 35.0618647531" Properties=species:S:1:pos:R:3:grain_id:I:1 pbc="T T T" |
| `Sigma5_310_grain_id.xyz` | 3,600 原子 · 附加列: Li x y z grain_id |
| `Sigma5_310_grainID.xyz` | 3,224 原子 · 附加列: Lattice="-35.1000000000 35.1000000000 70.2000000000 42.1200000000 -42.1200000000 21.0600000000 10.5300000000 10.5300000000 -10.5300000000" Properties=species:S:1:pos:R:3:grain_id:R:1 |
| `stage4a_metastable_disorder.xyz` | 2,559 原子 · 附加列: Ni D10 continuous SOAP disorder (visualization only) |

## 项目二：阶段一验收包（14400 原子 2x2）

`Li_Sigma5_310_2x2_阶段一验收_2026-09-05_new` — 39 个文件，16.59 MB

### 结论（10 个，11.43 MB）

| 文件 | 大小 | 用途 | 内容分析 |
|---|---|---|---|
| `0_初始XYZ与检测模块\Sigma5_310_K6_J6_L20.xyz` | 668.3 KB | 14400 原子 2×2 初始模型（K=6, J=6, L=20）。本验收包所有计算的起点。 | 14,400 原子 · 附加列: Lattice="66.5252109796 0.0 0.0 0.0 66.5252109796 0.0 0.0 0.0 70.1237295062" Properties=species:S:1:pos:R:3:grain_id:I:1 pbc="T T T" |
| `0_初始XYZ与检测模块\Sigma5_310_K6_J6_L20_grain_id.xyz` | 668.2 KB | 初始模型带晶粒标注的版本，用于确认两条晶界等价、晶粒内部无缺陷。 | 14,400 原子 · 附加列: Li x y z grain_id |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\1_零应力弛豫0K\p1_min0.data` | 1.09 MB | 验收步骤 A：0 K 全维度零应力弛豫的末态。判据是三个方向压力约 0、最大力 8.5e-6 eV/Å。 | 14,400 原子 · 盒子 66.24 × 66.18 × 69.73 Å · timestep 12,500 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\2_退火363K零压\p1_annealed.data` | 1.79 MB | ★ 验收步骤 B：363 K 等压退火 25 ps 后的末态。晶界重构就是在这 25 ps 里发生的。 | 14,400 原子 · 盒子 66.48 × 66.41 × 69.98 Å · timestep 25,000 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\3_末帧0K检查\p1_casecheck_0K.data` | 1.79 MB | ★ 验收步骤 C：取退火末帧做 0 K 检查。这是判定“Case B”（晶界真实无序、晶粒内部仍为 BCC）的关键结构。 | 14,400 原子 · 盒子 66.48 × 66.41 × 69.98 Å · timestep 934 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\4_淬火10K\p1_quenched.data` | 1.82 MB | 验收步骤 D：以 10 K/ps 从 363 K 淬火到 10 K 后的末态，历时 35.3 ps。 | 14,400 原子 · 盒子 66.17 × 66.11 × 69.66 Å · timestep 35,300 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\5_保持10K零压\p1_held.data` | 1.82 MB | ★ 验收的核心产物：10 K 保持结束时的结构。此时晶界无序、晶粒内部仍为 BCC。 | 14,400 原子 · 盒子 66.17 × 66.10 × 69.66 Å · timestep 5,000 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\5_保持10K零压\p1_held_0K.data` | 1.82 MB | 终态在 0 K 下再做一次复核的结构，对应验收步骤 F（晶界 38.1%/38.0% BCC，内部 98.7% BCC）。 | 14,400 原子 · 盒子 66.17 × 66.10 × 69.66 Å · timestep 187 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\5_保持10K零压\p1_hold_pressure.stats` | 171 B | ★ 保持段的压力时间序列。“10 K、0 bar 下维持”这个结论的直接证据来源。 | 2 个采样点 · 标量时均压力 -0.87 bar · 分量 Pxx=-187, Pyy=+383, Pzz=-199 · ⚠ 分量极差 582 bar，各向异性显著 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\6_结构分析\cna_phase1_report.txt` | 3.3 KB | ★ CNA 结构分析结论：各阶段 BCC 分数（晶界约 38%、晶粒内部 98.7%）。验收判据“晶界真实无序、内部保持 BCC”靠它成立。 | 晶界 BCC: initial 73.3% → p1_min0 60.0% → casecheck 43.7% → p1_final 38.1% · 晶内 BCC: 100.0% → 100.0% → 99.9% → 98.7% · 无序带宽: 2.5 Å → 4.5 Å → 4.5 Å → 6.5 Å |

### 续算（4 个，4.84 MB）

| 文件 | 大小 | 用途 |
|---|---|---|
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\1_零应力弛豫0K\p1_min0.restart` | 1.21 MB | LAMMPS 续算检查点：保存速度场与热力学状态，只能在 LAMMPS 内读取。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\2_退火363K零压\p1_annealed.restart` | 1.21 MB | LAMMPS 续算检查点：保存速度场与热力学状态，只能在 LAMMPS 内读取。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\4_淬火10K\p1_quenched.restart` | 1.21 MB | LAMMPS 续算检查点：保存速度场与热力学状态，只能在 LAMMPS 内读取。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\5_保持10K零压\p1_held.restart` | 1.21 MB | ★ 阶段二的起点。验收报告明确指定用这个检查点继续后续计算。 |

### 追溯（10 个，288.9 KB）

| 文件 | 大小 | 用途 |
|---|---|---|
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\7_运行日志\p1_anneal.barofail.lammps` | 25.8 KB | 柱压器失败的运行日志。保留了故障现场——本次运行修正了 Pdamp 参数，是问题定位的痕迹。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\7_运行日志\p1_anneal.lammps` | 25.8 KB | 该步骤的 LAMMPS 运行日志，含温度/压力/能量时间序列与断言检查。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\7_运行日志\p1_casecheck.barofail.lammps` | 12.0 KB | 柱压器失败的运行日志。保留了故障现场——本次运行修正了 Pdamp 参数，是问题定位的痕迹。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\7_运行日志\p1_casecheck.lammps` | 12.0 KB | 该步骤的 LAMMPS 运行日志，含温度/压力/能量时间序列与断言检查。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\7_运行日志\p1_finalcheck.lammps` | 8.3 KB | 该步骤的 LAMMPS 运行日志，含温度/压力/能量时间序列与断言检查。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\7_运行日志\p1_hold.barofail.lammps` | 17.3 KB | 柱压器失败的运行日志。保留了故障现场——本次运行修正了 Pdamp 参数，是问题定位的痕迹。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\7_运行日志\p1_hold.lammps` | 17.3 KB | 该步骤的 LAMMPS 运行日志，含温度/压力/能量时间序列与断言检查。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\7_运行日志\p1_min0.lammps` | 129.1 KB | 该步骤的 LAMMPS 运行日志，含温度/压力/能量时间序列与断言检查。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\7_运行日志\p1_quench.barofail.lammps` | 20.7 KB | 柱压器失败的运行日志。保留了故障现场——本次运行修正了 Pdamp 参数，是问题定位的痕迹。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\7_运行日志\p1_quench.lammps` | 20.7 KB | 该步骤的 LAMMPS 运行日志，含温度/压力/能量时间序列与断言检查。 |

### 代码（3 个，13.3 KB）

| 文件 | 大小 | 用途 |
|---|---|---|
| `0_初始XYZ与检测模块\cna_phase1.py` | 5.5 KB | CNA（共近邻分析）结构分析脚本，用于统计 BCC 分数。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\6_结构分析\cna_phase1.py` | 5.5 KB | CNA 结构分析脚本，统计各阶段晶界与晶粒内部的 BCC 比例。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\6_结构分析\run_phase1.sh` | 2.3 KB | 阶段一验收的完整编排脚本（依次跑五个步骤并收集结果）。 |

### 输入（6 个，7.3 KB）

| 文件 | 大小 | 用途 |
|---|---|---|
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\6_结构分析\in.p1_anneal` | 1.2 KB | 该步骤的 LAMMPS 输入脚本。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\6_结构分析\in.p1_casecheck` | 1.1 KB | 该步骤的 LAMMPS 输入脚本。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\6_结构分析\in.p1_finalcheck` | 1.1 KB | 该步骤的 LAMMPS 输入脚本。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\6_结构分析\in.p1_hold` | 1.3 KB | 该步骤的 LAMMPS 输入脚本。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\6_结构分析\in.p1_min0` | 1.6 KB | 该步骤的 LAMMPS 输入脚本。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\6_结构分析\in.p1_quench` | 1.0 KB | 该步骤的 LAMMPS 输入脚本。 |

### 文档（6 个，16.9 KB）

| 文件 | 大小 | 用途 |
|---|---|---|
| `0_初始XYZ与检测模块\README.md` | 448 B | 初始结构与检测模块的说明文档。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\README.md` | 1.2 KB | 阶段一结果包的入口说明。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\RUN_META.md` | 3.3 KB | 本次运行的元信息：时间、机器、参数、势函数标识。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\RUN_STATUS.md` | 4.6 KB | 本次运行的实时执行记录与各步骤状态。 |
| `Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\总结报告.md` | 5.4 KB | 阶段一零压运行的中文总结报告。 |
| `阶段一验收汇报_2026-09-04.md` | 2.0 KB | ★ 阶段一验收结论：原目标通过，可进入阶段二；并记录了 0 bar 口径的保留意见。 |

## 项目一：锂晶界主项目

`Lithium_GrainBoundary_Project` — 123 个文件，12.09 MB

### 结论（21 个，7.78 MB）

| 文件 | 大小 | 用途 | 内容分析 |
|---|---|---|---|
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\annealed\stage2_annealed.data` | 456.2 KB | ★ Stage 2 退火末态：363 K、0 bar 下 25 ps 后的结构。晶界在此阶段发生热重构。 | 3,600 原子 · 盒子 65.95 × 33.26 × 35.06 Å · timestep 20,000 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\bulk\bulk_Li.cif` | 6.8 KB | 体相参考结构的晶体学格式。 |  |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\bulk\bulk_Li.data` | 11.2 KB | 体相参考结构的 LAMMPS 格式（128 原子 BCC）。 | 128 原子 · 盒子 14.04 × 14.04 × 14.04 Å |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\bulk\bulk_Li.xyz` | 5.7 KB | 体相参考结构的通用格式。 | 128 原子 · 附加列: Lattice="14.0400000000 0.0 0.0 0.0 14.0400000000 0.0 0.0 0.0 14.0400000000" Properties=species:S:1:pos:R:3:grain_id:I:1 pbc="T T T" |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\bulk\bulk_Li_minimized.data` | 9.9 KB | ★ 体相参考的零压弛豫结果。给出该势函数下的平衡晶格常数 3.506186 Å —— 全部双晶都用它构建，避免把体相应变强加给晶界周期。 | 128 原子 · 盒子 14.02 × 14.02 × 14.02 Å · timestep 12 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\final\stage4_final.data` | 277.7 KB | ★ 全项目主产物：最终弛豫结构。单条退火-淬火-极小化轨迹收敛到的一个局域极小候选。 | 3,600 原子 · 盒子 66.01 × 33.26 × 35.06 Å · timestep 3,000 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\initial\Sigma5_310.cif` | 185.6 KB | 晶体学标准格式。给外部软件或合作者读，通用性最好。 |  |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\initial\Sigma5_310.data` | 322.3 KB | 初始双晶模型的 LAMMPS 格式（3600 原子）。整个项目的结构基准，后续每一级计算都从这里出发。 | 3,600 原子 · 盒子 66.53 × 33.26 × 35.06 Å |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\initial\Sigma5_310.grain_map` | 23.5 KB | 晶粒归属图，用于快速核对两条晶界的位置和对称性。 |  |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\initial\Sigma5_310.xyz` | 162.6 KB | 同一模型的扩展 XYZ 格式，比 .data 更通用，便于其他软件读取。 | 3,600 原子 · 附加列: Lattice="66.5252109796 0.0 0.0 0.0 33.2626054898 0.0 0.0 0.0 35.0618647531" Properties=species:S:1:pos:R:3:grain_id:I:1 pbc="T T T" |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\initial\Sigma5_310_grain_id.xyz` | 162.4 KB | 每原子标注属于哪个晶粒。用来验证两条晶界在晶体学上等价——这是“对称倾转晶界”成立的前提。 | 3,600 原子 · 附加列: Li x y z grain_id |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\minimized\stage1_minimized.data` | 277.9 KB | Stage 1 极小化后的结构：初始双晶在 0 K 下消除残余力的结果。 | 3,600 原子 · 盒子 65.95 × 33.26 × 35.06 Å · timestep 5,000 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\quenched\stage3_quenched.data` | 277.7 KB | Stage 3 淬火末态：以 10 K/ps 降到 10 K 后的结构。无序被快速“冻”住。 | 3,600 原子 · 盒子 65.95 × 33.26 × 35.06 Å · timestep 50,000 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\validation\rebuilt_stage2.data` | 234.8 KB | 从 dump 轨迹反推重建的结构。用来证明 dump 没有丢信息——是“模拟结果可复现”的直接证据。 | 3,600 原子 · 盒子 65.95 × 33.26 × 35.06 Å · timestep 20,000 |
| `Li_Sigma5_310_metastable_2x2\stage0_build_validate\Sigma5_310_K6_J6_L20.cif` | 763.0 KB | 2×2 模型的晶体学标准格式。 |  |
| `Li_Sigma5_310_metastable_2x2\stage0_build_validate\Sigma5_310_K6_J6_L20.data` | 1.30 MB | 2×2 模型的 Stage 0 构建产物（14400 原子），是验收包的前身。 | 14,400 原子 · 盒子 66.53 × 66.53 × 70.12 Å |
| `Li_Sigma5_310_metastable_2x2\stage0_build_validate\Sigma5_310_K6_J6_L20.grain_map` | 115.7 KB | 晶粒归属图。 |  |
| `Li_Sigma5_310_metastable_2x2\stage0_build_validate\Sigma5_310_K6_J6_L20.xyz` | 668.3 KB | 2×2 模型的通用结构格式。 | 14,400 原子 · 附加列: Lattice="66.5252109796 0.0 0.0 0.0 66.5252109796 0.0 0.0 0.0 70.1237295062" Properties=species:S:1:pos:R:3:grain_id:I:1 pbc="T T T" |
| `Li_Sigma5_310_metastable_2x2\stage0_build_validate\Sigma5_310_K6_J6_L20_grain_id.xyz` | 668.2 KB | 扩展 XYZ 结构文件，可携带每原子附加列。 | 14,400 原子 · 附加列: Li x y z grain_id |
| `Li_Sigma5_310_metastable_2x2\stage0_build_validate\structure.data` | 1.30 MB | 初始结构的 LAMMPS 格式。 | 14,400 原子 · 盒子 66.53 × 66.53 × 70.12 Å |
| `Li_Sigma5_310_metastable_2x2\stage0_build_validate\structure.xyz` | 668.3 KB | 初始结构的通用格式。 | 14,400 原子 · 附加列: Lattice="66.5252109796 0.0 0.0 0.0 66.5252109796 0.0 0.0 0.0 70.1237295062" Properties=species:S:1:pos:R:3:grain_id:I:1 pbc="T T T" |

### 续算（5 个，1.22 MB）

| 文件 | 大小 | 用途 |
|---|---|---|
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\annealed\stage2_annealed.restart` | 310.3 KB | LAMMPS 续算检查点：保存速度场与热力学状态，只能在 LAMMPS 内读取。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\bulk\bulk_Li_minimized.restart` | 11.9 KB | 体相弛豫的续算检查点。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\final\stage4_final.restart` | 310.3 KB | 最终结构的续算检查点，含速度与热力学状态。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\minimized\stage1_minimized.restart` | 310.3 KB | LAMMPS 续算检查点：保存速度场与热力学状态，只能在 LAMMPS 内读取。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\quenched\stage3_quenched.restart` | 310.3 KB | LAMMPS 续算检查点：保存速度场与热力学状态，只能在 LAMMPS 内读取。 |

### 追溯（24 个，525.2 KB）

| 文件 | 大小 | 用途 | 内容分析 |
|---|---|---|---|
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\initial\params.json` | 557 B | 构建参数：超胞尺寸 K/J/L、晶格常数、势函数名。复现这个模型所需的最小信息集。 | K=6 · J=3 · L=10 · a_A=3.506186475312 · theta_deg=36.86989764584402 · final_total=3600 · total_before_removal=3720 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\data\initial\removal_log.csv` | 8.2 KB | 重合位点删除记录：每个晶界删了哪些原子、删了多少。保证删除对称、可审计，是可追溯性的核心。 | 120 行 × 8 列 [atom_id, grain_id, x, y, z, nearest_distance_A …] |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\docs\reports\ENVIRONMENT.txt` | 180 B | 运行环境快照（LAMMPS 版本、MPI、势函数等），复现的必要条件。 |  |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\docs\reports\SHA256SUMS.txt` | 7.3 KB | 本包内文件的校验清单。日后任何人可据此确认文件未被改动。 | 收录 80 个文件的校验和 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\logs\bulk_energy.lammps` | 8.9 KB | 体相能量计算的运行日志，给出每原子势能参考值。 |  |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\logs\dump2data_validation.lammps` | 1.7 KB | dump→data 重建验证的运行日志，用于确认轨迹文件未丢信息。 |  |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\logs\kim.log` | 277.9 KB | OpenKIM 势函数调用日志，记录实际加载的模型标识。 |  |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\logs\stage0_check.done` | 39 B | 阶段完成标记，供编排脚本判断是否跳过。 |  |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\logs\stage0_check_a0_mpi12.lammps` | 6.8 KB | 该阶段的 LAMMPS 运行日志，含数值历史与断言检查。 |  |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\logs\stage1_minimize.done` | 39 B | 阶段完成标记，供编排脚本判断是否跳过。 |  |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\logs\stage1_minimize_checkpoint_mpi12.lammps` | 55.9 KB | 该阶段的 LAMMPS 运行日志，含数值历史与断言检查。 |  |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\logs\stage2_anneal.done` | 39 B | 阶段完成标记，供编排脚本判断是否跳过。 |  |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\logs\stage2_anneal_mpi12.lammps` | 15.4 KB | 该阶段的 LAMMPS 运行日志，含数值历史与断言检查。 |  |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\logs\stage2_anneal_mpi12.screen` | 10.9 KB | 进程标准输出记录。 |  |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\logs\stage3_quench.done` | 39 B | 阶段完成标记，供编排脚本判断是否跳过。 |  |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\logs\stage3_quench_mpi12.lammps` | 17.1 KB | 该阶段的 LAMMPS 运行日志，含数值历史与断言检查。 |  |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\logs\stage3_quench_mpi12.screen` | 12.6 KB | 进程标准输出记录。 |  |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\logs\stage4_final.done` | 39 B | 阶段完成标记，供编排脚本判断是否跳过。 |  |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\logs\stage4_final_mpi12.lammps` | 37.1 KB | 该阶段的 LAMMPS 运行日志，含数值历史与断言检查。 |  |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\logs\stage4_final_mpi12.screen` | 28.6 KB | 进程标准输出记录。 |  |
| `Li_Sigma5_310_metastable_2x2\PARAMETERS.json` | 1.9 KB | 2×2 扩展的参数记录：K 保持 6，J 由 3 加到 6，L 由 10 加到 20。 | 字段: experiment, created, material, gb, gb_plane_label, misorientation_deg, potential, lattice_constant_A |
| `Li_Sigma5_310_metastable_2x2\stage0_build_validate\parameters.json` | 563 B | 构建参数记录（副本）。 | K=6 · J=6 · L=20 · a_A=3.506186475312 · theta_deg=36.86989764584402 · final_total=14400 · total_before_removal=14880 |
| `Li_Sigma5_310_metastable_2x2\stage0_build_validate\params.json` | 563 B | 构建参数记录。 | K=6 · J=6 · L=20 · a_A=3.506186475312 · theta_deg=36.86989764584402 · final_total=14400 · total_before_removal=14880 |
| `Li_Sigma5_310_metastable_2x2\stage0_build_validate\removal_log.csv` | 33.0 KB | 2×2 模型的重合位点删除记录。 | 480 行 × 8 列 [atom_id, grain_id, x, y, z, nearest_distance_A …] |

### 代码（38 个，152.2 KB）

| 文件 | 大小 | 用途 |
|---|---|---|
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\lammps\resume.sh` | 2.0 KB | 从中断处续算的脚本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\lammps\run_all.sh` | 1.3 KB | 完整热历史编排脚本；带 .done 标记实现断点续算。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\src\generator\__init__.py` | 86 B | Python 脚本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\src\generator\builder.py` | 8.8 KB | Python 脚本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\src\generator\bulk.py` | 2.3 KB | Python 脚本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\src\generator\config.py` | 2.1 KB | Python 脚本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\src\generator\exporters.py` | 5.6 KB | Python 脚本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\src\generator\generate.py` | 10.1 KB | Python 脚本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\src\generator\geometry.py` | 9.6 KB | Python 脚本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\src\generator\lammps_io.py` | 3.4 KB | Python 脚本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\src\generator\reimport_check.py` | 4.2 KB | Python 脚本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\src\generator\verify.py` | 16.7 KB | Python 脚本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\tests\conftest.py` | 223 B | 测试公共夹具。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\tests\test_builder.py` | 2.5 KB | Python 脚本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\tests\test_conventional_bcc.py` | 1.1 KB | Python 脚本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\tests\test_grain_id.py` | 1.3 KB | Python 脚本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\tests\test_legacy_regression.py` | 1.0 KB | Python 脚本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\tests\test_no_nvt_zero.py` | 922 B | Python 脚本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\tests\test_outputs_consistency.py` | 1.5 KB | Python 脚本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\tests\test_overlap_removal.py` | 1.2 KB | Python 脚本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\tests\test_periodicity.py` | 1.4 KB | Python 脚本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\tests\test_sigma5_geometry.py` | 1.7 KB | Python 脚本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\tests\test_velocities_zero.py` | 1.2 KB | Python 脚本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\tests\test_verify.py` | 2.5 KB | Python 脚本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\tests\test_workflow_scripts.py` | 3.7 KB | Python 脚本。 |
| `Li_Sigma5_310_metastable_2x2\generator\__init__.py` | 86 B | Python 脚本。 |
| `Li_Sigma5_310_metastable_2x2\generator\builder.py` | 8.8 KB | Python 脚本。 |
| `Li_Sigma5_310_metastable_2x2\generator\bulk.py` | 2.3 KB | Python 脚本。 |
| `Li_Sigma5_310_metastable_2x2\generator\config.py` | 2.1 KB | Python 脚本。 |
| `Li_Sigma5_310_metastable_2x2\generator\exporters.py` | 5.6 KB | Python 脚本。 |
| `Li_Sigma5_310_metastable_2x2\generator\generate.py` | 10.1 KB | Python 脚本。 |
| `Li_Sigma5_310_metastable_2x2\generator\geometry.py` | 9.6 KB | Python 脚本。 |
| `Li_Sigma5_310_metastable_2x2\generator\lammps_io.py` | 3.4 KB | Python 脚本。 |
| `Li_Sigma5_310_metastable_2x2\generator\reimport_check.py` | 4.2 KB | Python 脚本。 |
| `Li_Sigma5_310_metastable_2x2\generator\verify.py` | 16.7 KB | Python 脚本。 |
| `Li_Sigma5_310_metastable_2x2\lammps\resume.sh` | 1.2 KB | 从中断处续算的脚本。 |
| `Li_Sigma5_310_metastable_2x2\lammps\run_all.sh` | 1.3 KB | 完整热历史编排脚本；带 .done 标记实现断点续算。 |
| `Li_Sigma5_310_metastable_2x2\scripts\run_windows.ps1` | 788 B | Windows 侧的辅助运行脚本。 |

### 输入（15 个，13.9 KB）

| 文件 | 大小 | 用途 |
|---|---|---|
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\lammps\in.bulk_energy` | 990 B | 体相能量参考计算：给出每原子势能基准，用于扣除体相贡献、得到晶界过剩量。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\lammps\in.dump2data` | 546 B | dump→data 重建脚本：把轨迹转回结构文件，用来验证轨迹未丢失信息。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\lammps\stage0_test\in.0_check` | 662 B | Stage 0 输入脚本：核对模型几何、原子数与势函数是否被正确加载。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\lammps\stage1_minimize\in.1_minimize` | 1.5 KB | Stage 1 输入脚本：0 K 极小化，仅沿晶界法向做零压弛豫（y/z 固定在零压 CSL 周期），消除初始残余力。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\lammps\stage2_anneal\in.2_anneal` | 953 B | Stage 2 输入脚本：363 K 等压（NPT）退火，是诱发晶界热重构的关键一步。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\lammps\stage2_anneal\in.2_anneal_resume` | 1.1 KB | Stage 2 的续算版本：从中断的检查点继续退火，不重跑已完成的时段。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\lammps\stage3_quench\in.3_quench` | 964 B | Stage 3 输入脚本：以 10 K/ps 从 363 K 淬火到 10 K，把高温无序结构冻结下来。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\lammps\stage3_quench\in.3_quench_resume` | 1.0 KB | Stage 3 的续算版本。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\lammps\stage4_relax\in.4_final` | 1.6 KB | Stage 4 输入脚本：最终弛豫与保持，产出用于结构分析的末态。 |
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\requirements.txt` | 34 B | Python 依赖清单，用于重建运行环境。 |
| `Li_Sigma5_310_metastable_2x2\lammps\in.0_check` | 514 B | Stage 0 输入脚本：核对模型几何、原子数与势函数是否被正确加载。 |
| `Li_Sigma5_310_metastable_2x2\lammps\in.1_minimize` | 1.4 KB | Stage 1 输入脚本：0 K 极小化，仅沿晶界法向做零压弛豫（y/z 固定在零压 CSL 周期），消除初始残余力。 |
| `Li_Sigma5_310_metastable_2x2\lammps\in.2_highT_npt` | 995 B | 2×2 实验的高温 NPT 输入脚本（对应 in.2_anneal 在旧版命名下的同一角色）。 |
| `Li_Sigma5_310_metastable_2x2\lammps\in.3_quench` | 940 B | Stage 3 输入脚本：以 10 K/ps 从 363 K 淬火到 10 K，把高温无序结构冻结下来。 |
| `Li_Sigma5_310_metastable_2x2\lammps\in.4_10K_hold` | 863 B | 2×2 实验的 10 K 保持输入脚本（对应 in.4_final 在旧版命名下的同一角色）。 |

### 文档（20 个，2.40 MB）

| 文件 | 大小 | 用途 |
|---|---|---|
| `01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\tests\pytest.ini` | 40 B | 测试框架配置。 |
| `07_Documentation\presentations\Li_metal_metastable_GB_atomic_mechanical_analysis.pptx` | 2.20 MB | 汇报幻灯片。 |
| `AGENTS.md` | 3.1 KB | Markdown 文档。 |
| `Li_Sigma5_310_metastable_2x2\stage0_build_validate\grain_map` | 101.7 KB | 2×2 模型的晶粒归属图。 |
| `Reports\GB_Construction\FINAL_STATUS_BI.md` | 1.0 KB | 双语版报告（中文版与英文版首尾拼接在同一文件内）。 |
| `Reports\GB_Construction\README_BI.md` | 20.7 KB | 双语版报告（中文版与英文版首尾拼接在同一文件内）。 |
| `Reports\GB_Construction\build_report_BI.md` | 4.4 KB | 双语版报告（中文版与英文版首尾拼接在同一文件内）。 |
| `Reports\GB_Construction\verification_report_BI.md` | 3.8 KB | 双语版报告（中文版与英文版首尾拼接在同一文件内）。 |
| `Reports\Li_2x2_Workflow\PROTOCOL_BI.md` | 10.5 KB | 双语版报告（中文版与英文版首尾拼接在同一文件内）。 |
| `Reports\Li_2x2_Workflow\PROVENANCE_BI.md` | 4.2 KB | 双语版报告（中文版与英文版首尾拼接在同一文件内）。 |
| `Reports\Li_2x2_Workflow\README_BI.md` | 3.2 KB | 双语版报告（中文版与英文版首尾拼接在同一文件内）。 |
| `Reports\Li_2x2_Workflow\RUN_STATUS_BI.md` | 4.6 KB | 双语版报告（中文版与英文版首尾拼接在同一文件内）。 |
| `Reports\Li_2x2_Workflow\stage0_build_report_BI.md` | 4.8 KB | 双语版报告（中文版与英文版首尾拼接在同一文件内）。 |
| `Reports\Li_2x2_Workflow\stage0_validation_report_BI.md` | 3.8 KB | 双语版报告（中文版与英文版首尾拼接在同一文件内）。 |
| `Reports\Li_2x2_Workflow\stage0_verification_report_BI.md` | 3.8 KB | 双语版报告（中文版与英文版首尾拼接在同一文件内）。 |
| `Reports\Li_2x2_Workflow\阶段一_原目标_执行指南_2026-09-02_BI.md` | 7.4 KB | 双语版报告（中文版与英文版首尾拼接在同一文件内）。 |
| `Reports\README_BI.md` | 796 B | 双语版报告（中文版与英文版首尾拼接在同一文件内）。 |
| `Reports\Team_Reviews\TEAM_REVIEW_Li_GB_Deliverable_20260823_BI.md` | 14.5 KB | 双语版报告（中文版与英文版首尾拼接在同一文件内）。 |
| `Reports\Team_Reviews\TEAM_REVIEW_Li_Sigma5_2x2_run_20260830_BI.md` | 11.0 KB | 双语版报告（中文版与英文版首尾拼接在同一文件内）。 |
| `Reports\Team_Reviews\TEAM_REVIEW_Li_Sigma5_310_metastable_2x2_20260901_BI.md` | 8.1 KB | 双语版报告（中文版与英文版首尾拼接在同一文件内）。 |

## 项目三：Ni SOAP 方法探索

`Ni_SOAP_Method_Exploration_by_Stage` — 190 个文件，22.15 MB

### 结论（150 个，13.62 MB）

| 文件 | 大小 | 用途 | 内容分析 |
|---|---|---|---|
| `00_汇总报告\项目文件索引_按阶段.csv` | 60.2 KB | 数据表。 | 209 行 × 11 列 [Stage, Category, Subcategory, RelativePath, FileName, Extension …] |
| `Stage01_数据准备与SOAP初测\03_输出\inspection_summary.txt` | 2.8 KB | 纯文本结果或检查报告。 |  |
| `Stage01_数据准备与SOAP初测\03_输出\soap_test_summary.txt` | 1.4 KB | 纯文本结果或检查报告。 |  |
| `Stage01_数据准备与SOAP初测\03_输出\stage1_qc_report.txt` | 7.3 KB | Stage01 数据准备与初测产物。 |  |
| `Stage02_SOAP表示与审计\03_输出\representative_atoms.csv` | 4.8 KB | 数据表。 | 33 行 × 12 列 [structure, atom_index, id, ptm, cna, peatom …] |
| `Stage02_SOAP表示与审计\03_输出\soap_similarity_heatmap.png` | 181.7 KB | 图表，可直接用于报告或论文。 |  |
| `Stage02_SOAP表示与审计\03_输出\soap_similarity_matrix.csv` | 11.4 KB | 数据表。 | 33 行 × 34 列 [label, stable|15487|bulk, stable|3157|bulk, stable|2287|bulk, stable|8226|bulk, stable|8836|non_dominant …] |
| `Stage02_SOAP表示与审计\03_输出\spatial_profile_metastable.png` | 121.9 KB | 图表，可直接用于报告或论文。 |  |
| `Stage02_SOAP表示与审计\03_输出\spatial_profile_stable.png` | 109.9 KB | 图表，可直接用于报告或论文。 |  |
| `Stage02_SOAP表示与审计\03_输出\spatial_profile_summary.txt` | 1.3 KB | 纯文本结果或检查报告。 |  |
| `Stage02_SOAP表示与审计\03_输出\stage2_5_atom_metadata.csv` | 4.4 KB | Stage02.5 表示审计产物。 | 33 行 × 12 列 [row_index, structure, atom_index, atom_id, x, y …] |
| `Stage02_SOAP表示与审计\03_输出\stage2_5_metric_heatmaps.png` | 218.3 KB | Stage02.5 表示审计产物。 |  |
| `Stage02_SOAP表示与审计\03_输出\stage2_5_pair_category_summary.csv` | 4.0 KB | Stage02.5 表示审计产物。 | 28 行 × 8 列 [category, metric, n_pairs, min, median, mean …] |
| `Stage02_SOAP表示与审计\03_输出\stage2_5_pairwise_metrics.csv` | 104.2 KB | Stage02.5 表示审计产物。 | 528 行 × 20 列 [row_i, structure_i, atom_id_i, ptm_i, cna_i, peatom_i …] |
| `Stage02_SOAP表示与审计\03_输出\stage2_5_parameter_sensitivity.csv` | 254 B | Stage02.5 表示审计产物。 | 4 行 × 7 列 [n_max, l_max, feature_dimension, linalg_warning, warning_count, vectors_finite …] |
| `Stage02_SOAP表示与审计\03_输出\stage2_75_atom_metadata.csv` | 1.17 MB | Stage02.75 全人口 SOAP 与距离审计产物。 | 4,659 行 × 19 列 [structure, region, atom_id, x, y, z …] |
| `Stage02_SOAP表示与审计\03_输出\stage2_75_bulk_distance_distributions.png` | 35.4 KB | Stage02.75 全人口 SOAP 与距离审计产物。 |  |
| `Stage02_SOAP表示与审计\03_输出\stage2_75_distance_decomposition.png` | 161.5 KB | ★ 距离分解图：把 SOAP 距离拆成范数贡献与形状贡献。这是“原始距离含义不纯”的关键证据图。 |  |
| `Stage02_SOAP表示与审计\03_输出\stage2_75_distance_vs_PEATOM.png` | 149.4 KB | Stage02.75 全人口 SOAP 与距离审计产物。 |  |
| `Stage02_SOAP表示与审计\03_输出\stage2_75_distance_vs_z.png` | 111.7 KB | Stage02.75 全人口 SOAP 与距离审计产物。 |  |
| `Stage02_SOAP表示与审计\03_输出\stage2_75_norm_distributions.png` | 29.2 KB | Stage02.75 全人口 SOAP 与距离审计产物。 |  |
| `Stage02_SOAP表示与审计\03_输出\stage2_75_pairwise_summary.csv` | 4.5 KB | Stage02.75 全人口 SOAP 与距离审计产物。 | 21 行 × 12 列 [category, metric, N, min, p5, p25 …] |
| `Stage02_SOAP表示与审计\03_输出\stage2_75_population_summary.csv` | 7.7 KB | Stage02.75 全人口 SOAP 与距离审计产物。 | 40 行 × 12 列 [category, metric, N, min, p5, p25 …] |
| `Stage02_SOAP表示与审计\03_输出\stage2_75_ptm_cna_summary.csv` | 2.8 KB | Stage02.75 全人口 SOAP 与距离审计产物。 | 16 行 × 12 列 [structure, group_type, group_value, N, soap_norm_median, soap_norm_iqr …] |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3a_bulk_behavior_vs_epsilon.png` | 32.0 KB | Stage03A 阈值扫描产物：衡量 epsilon 取值对结果的影响。 |  |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3a_class_count_vs_epsilon.png` | 21.0 KB | Stage03A 阈值扫描产物：衡量 epsilon 取值对结果的影响。 |  |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3a_epsilon_grid.csv` | 5.0 KB | Stage03A 阈值扫描产物：衡量 epsilon 取值对结果的影响。 | 115 行 × 3 列 [index, epsilon, source] |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3a_epsilon_sweep.csv` | 89.4 KB | epsilon 阈值扫描结果。证明了阈值轻微变动就会改变分类数——离散 LAE 不稳。 | 115 行 × 55 列 [epsilon, n_classes_mean, n_classes_std, singleton_class_fraction_mean, singleton_class_fraction_std, atom_singleton_fraction_mean …] |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3a_nearest_neighbor_distances.csv` | 331.0 KB | Stage03A 阈值扫描产物：衡量 epsilon 取值对结果的影响。 | 4,659 行 × 9 列 [row_index, structure, region, atom_id, nearest_neighbor_row, nearest_neighbor_structure …] |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3a_nearest_neighbor_summary.csv` | 2.0 KB | Stage03A 阈值扫描产物：衡量 epsilon 取值对结果的影响。 | 7 行 × 15 列 [group, N, min, p1, p5, p10 …] |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3a_nn_distributions.png` | 25.1 KB | Stage03A 阈值扫描产物：衡量 epsilon 取值对结果的影响。 |  |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3a_order_dependence_vs_epsilon.png` | 22.4 KB | Stage03A 阈值扫描产物：衡量 epsilon 取值对结果的影响。 |  |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3a_partition_stability.csv` | 2.4 KB | Stage03A 阈值扫描产物：衡量 epsilon 取值对结果的影响。 | 15 行 × 11 列 [epsilon, epsilon_index, n_pairs, n_classes_mean, ari_mean, ari_median …] |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3a_partition_stability_vs_epsilon.png` | 24.8 KB | Stage03A 阈值扫描产物：衡量 epsilon 取值对结果的影响。 |  |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3a_refined_epsilon_sweep.csv` | 28.1 KB | Stage03A 阈值扫描产物：衡量 epsilon 取值对结果的影响。 | 35 行 × 55 列 [epsilon, n_classes_mean, n_classes_std, singleton_class_fraction_mean, singleton_class_fraction_std, atom_singleton_fraction_mean …] |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3a_refined_threshold_plot.png` | 25.9 KB | Stage03A 阈值扫描产物：衡量 epsilon 取值对结果的影响。 |  |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3a_singleton_fraction_vs_epsilon.png` | 19.2 KB | Stage03A 阈值扫描产物：衡量 epsilon 取值对结果的影响。 |  |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3a_threshold_candidates.csv` | 281 B | Stage03A 阈值扫描产物：衡量 epsilon 取值对结果的影响。 | 1 行 × 10 列 [interval_id, epsilon_min, epsilon_max, n_grid_points, n_classes_median_min, n_classes_median_max …] |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3b_atomwise_mid_epsilon.csv` | 611.5 KB | Stage03B 类别构建产物：LAE 分类与 LER 人口指纹。 | 4,659 行 × 14 列 [row_index, atom_id, x, y, z, class_id …] |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3b_canonical_classes.csv` | 2.6 KB | Stage03B 类别构建产物：LAE 分类与 LER 人口指纹。 | 25 行 × 12 列 [epsilon, class_id, representative_atom_id, total, stable_interface, metastable_interface …] |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3b_class_compactness.csv` | 3.0 KB | Stage03B 类别构建产物：LAE 分类与 LER 人口指纹。 | 25 行 × 9 列 [epsilon, class_id, n_atoms, n_pairs, diameter, median_within_distance …] |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3b_class_compactness.png` | 18.2 KB | Stage03B 类别构建产物：LAE 分类与 LER 人口指纹。 |  |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3b_class_population_vs_epsilon.png` | 18.3 KB | Stage03B 类别构建产物：LAE 分类与 LER 人口指纹。 |  |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3b_class_statistics.csv` | 33.8 KB | Stage03B 类别构建产物：LAE 分类与 LER 人口指纹。 | 821 行 × 5 列 [epsilon, class_id, stat, group, value] |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3b_epsilon_set.csv` | 544 B | Stage03B 类别构建产物：LAE 分类与 LER 人口指纹。 | 5 行 × 4 列 [epsilon, role, source, rationale] |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3b_lae_vs_z_metastable.png` | 13.1 KB | Stage03B 类别构建产物：LAE 分类与 LER 人口指纹。 |  |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3b_lae_vs_z_stable.png` | 13.2 KB | Stage03B 类别构建产物：LAE 分类与 LER 人口指纹。 |  |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3b_ler_by_run.csv` | 42.1 KB | Stage03B 类别构建产物：LAE 分类与 LER 人口指纹。 | 772 行 × 8 列 [epsilon, run, seed, class_id_in_run, class_ordinal_in_run, occupancy …] |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3b_ler_comparison.csv` | 13.1 KB | Stage03B 类别构建产物：LAE 分类与 LER 人口指纹。 | 155 行 × 10 列 [epsilon, run, seed, level, l1, jsd …] |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3b_ler_stable_vs_metastable.png` | 21.6 KB | Stage03B 类别构建产物：LAE 分类与 LER 人口指纹。 |  |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3b_ler_summary.csv` | 3.5 KB | Stage03B 类别构建产物：LAE 分类与 LER 人口指纹。 | 25 行 × 11 列 [epsilon, class_id, class_rank, n_runs_matched, occupancy_median_medoid, occupancy_mean …] |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3b_ler_vs_epsilon.png` | 27.1 KB | Stage03B 类别构建产物：LAE 分类与 LER 人口指纹。 |  |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3b_mid_epsilon_lae.xyz` | 543.2 KB | Stage03B 类别构建产物：LAE 分类与 LER 人口指纹。 | 4,659 原子 · 附加列: Ni LAE_LER representation for epsilon_mid (visual inspection only) |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3b_partition_ari.csv` | 67.3 KB | ★ 不同随机顺序下分区的一致性（ARI）。这是判定贪心 LAE “有顺序依赖、不可作为物理分类”的量化依据。 | 2,175 行 × 4 列 [epsilon, run_i, run_j, ari] |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3b_representative_atoms.csv` | 6.5 KB | Stage03B 类别构建产物：LAE 分类与 LER 人口指纹。 | 50 行 × 13 列 [epsilon, class_id, role, atom_id, structure, region …] |
| `Stage03_LAE阈值与LAE_LER构建\03_输出\stage3b_run_summary.csv` | 7.4 KB | Stage03B 类别构建产物：LAE 分类与 LER 人口指纹。 | 150 行 × 7 列 [epsilon, run, seed, n_classes, mean_ari_to_others, min_ari_to_others …] |
| `Stage04_无序度与能量关联\03_输出\stage4a_atomwise_disorder.csv` | 997.6 KB | ★ 逐原子无序度 D。本项目提出的连续结构坐标，是后续所有关联分析的自变量。 | 4,659 行 × 17 列 [structure, region, atom_id, x, y, z …] |
| `Stage04_无序度与能量关联\03_输出\stage4a_bootstrap_summary.csv` | 942 B | Stage04A 连续无序度产物。 | 6 行 × 8 列 [parameter, observed, bootstrap_mean, bootstrap_std, ci_low, ci_high …] |
| `Stage04_无序度与能量关联\03_输出\stage4a_bulk_disorder_distribution.png` | 19.7 KB | Stage04A 连续无序度产物。 |  |
| `Stage04_无序度与能量关联\03_输出\stage4a_four_group_disorder.png` | 23.8 KB | Stage04A 连续无序度产物。 |  |
| `Stage04_无序度与能量关联\03_输出\stage4a_group_statistics.csv` | 3.4 KB | Stage04A 连续无序度产物。 | 16 行 × 13 列 [group, k, n, mean, median, std …] |
| `Stage04_无序度与能量关联\03_输出\stage4a_interface_disorder_distribution.png` | 22.8 KB | Stage04A 连续无序度产物。 |  |
| `Stage04_无序度与能量关联\03_输出\stage4a_k_robustness.csv` | 1.5 KB | Stage04A 连续无序度产物。 | 15 行 × 7 列 [scope, k_a, k_b, pearson, spearman, median_abs_diff …] |
| `Stage04_无序度与能量关联\03_输出\stage4a_k_sensitivity.png` | 17.2 KB | Stage04A 连续无序度产物。 |  |
| `Stage04_无序度与能量关联\03_输出\stage4a_metastable_disorder.xyz` | 136.2 KB | Stage04A 连续无序度产物。 | 2,559 原子 · 附加列: Ni D10 continuous SOAP disorder (visualization only) |
| `Stage04_无序度与能量关联\03_输出\stage4a_pca_disorder.png` | 61.5 KB | Stage04A 连续无序度产物。 |  |
| `Stage04_无序度与能量关联\03_输出\stage4a_pca_groups.png` | 31.8 KB | Stage04A 连续无序度产物。 |  |
| `Stage04_无序度与能量关联\03_输出\stage4a_pca_peatom.png` | 60.9 KB | Stage04A 连续无序度产物。 |  |
| `Stage04_无序度与能量关联\03_输出\stage4a_pca_variance.csv` | 459 B | Stage04A 连续无序度产物。 | 10 行 × 3 列 [pc, explained_variance_ratio, cumulative] |
| `Stage04_无序度与能量关联\03_输出\stage4a_reference_sensitivity.csv` | 1.2 KB | Stage04A 连续无序度产物。 | 10 行 × 7 列 [comparison, k, scope, pearson, spearman, median_abs_diff …] |
| `Stage04_无序度与能量关联\03_输出\stage4a_stable_disorder.xyz` | 111.7 KB | Stage04A 连续无序度产物。 | 2,100 原子 · 附加列: Ni D10 continuous SOAP disorder (visualization only) |
| `Stage04_无序度与能量关联\03_输出\stage4b_all_atoms_scatter.png` | 30.6 KB | Stage04B 无序度-能量关联产物。 |  |
| `Stage04_无序度与能量关联\03_输出\stage4b_atomwise_energy_disorder.csv` | 1.26 MB | 逐原子的无序度与能量配对表，用于考察 D 与能量的关系。 | 4,659 行 × 21 列 [structure, region, atom_id, x, y, z …] |
| `Stage04_无序度与能量关联\03_输出\stage4b_binned_interface_trend.csv` | 2.6 KB | Stage04B 无序度-能量关联产物。 | 18 行 × 9 列 [bin_id, n, disorder_low, disorder_high, median_disorder, median_dE …] |
| `Stage04_无序度与能量关联\03_输出\stage4b_binned_metastable_trend.csv` | 1.7 KB | Stage04B 无序度-能量关联产物。 | 12 行 × 9 列 [bin_id, n, disorder_low, disorder_high, median_disorder, median_dE …] |
| `Stage04_无序度与能量关联\03_输出\stage4b_binned_stable_trend.csv` | 1.7 KB | Stage04B 无序度-能量关联产物。 | 12 行 × 9 列 [bin_id, n, disorder_low, disorder_high, median_disorder, median_dE …] |
| `Stage04_无序度与能量关联\03_输出\stage4b_bootstrap_summary.csv` | 3.2 KB | Stage04B 无序度-能量关联产物。 | 16 行 × 9 列 [parameter, observed, bootstrap_mean, bootstrap_std, ci_low, ci_high …] |
| `Stage04_无序度与能量关联\03_输出\stage4b_common_support_comparison.csv` | 1.1 KB | Stage04B 无序度-能量关联产物。 | 8 行 × 9 列 [bin_id, disorder_low, disorder_high, median_disorder, n_stable, n_metastable …] |
| `Stage04_无序度与能量关联\03_输出\stage4b_conditional_state_difference.png` | 21.3 KB | Stage04B 无序度-能量关联产物。 |  |
| `Stage04_无序度与能量关联\03_输出\stage4b_correlations.csv` | 1.2 KB | Stage04B 无序度-能量关联产物。 | 18 行 × 5 列 [baseline, scope, n, pearson, spearman] |
| `Stage04_无序度与能量关联\03_输出\stage4b_energy_baselines.csv` | 90 B | Stage04B 无序度-能量关联产物。 | 3 行 × 3 列 [baseline, n, median_peatom] |
| `Stage04_无序度与能量关联\03_输出\stage4b_fixed_disorder_energy_spread.csv` | 761 B | Stage04B 无序度-能量关联产物。 | 5 行 × 8 列 [disorder_low, disorder_high, n, mean_dE, median_dE, std_dE …] |
| `Stage04_无序度与能量关联\03_输出\stage4b_group_energy_statistics.csv` | 877 B | Stage04B 无序度-能量关联产物。 | 4 行 × 11 列 [group, n, mean, median, std, q1 …] |
| `Stage04_无序度与能量关联\03_输出\stage4b_interface_binned_trend.png` | 22.5 KB | Stage04B 无序度-能量关联产物。 |  |
| `Stage04_无序度与能量关联\03_输出\stage4b_interface_scatter.png` | 28.8 KB | Stage04B 无序度-能量关联产物。 |  |
| `Stage04_无序度与能量关联\03_输出\stage4b_k_robustness.csv` | 601 B | Stage04B 无序度-能量关联产物。 | 9 行 × 5 列 [k, scope, n, pearson, spearman] |
| `Stage04_无序度与能量关联\03_输出\stage4b_metastable_interface_scatter.png` | 27.3 KB | Stage04B 无序度-能量关联产物。 |  |
| `Stage04_无序度与能量关联\03_输出\stage4b_pca_disorder.png` | 61.5 KB | Stage04B 无序度-能量关联产物。 |  |
| `Stage04_无序度与能量关联\03_输出\stage4b_pca_energy.png` | 62.3 KB | Stage04B 无序度-能量关联产物。 |  |
| `Stage04_无序度与能量关联\03_输出\stage4b_pca_interface_energy.png` | 60.4 KB | Stage04B 无序度-能量关联产物。 |  |
| `Stage04_无序度与能量关联\03_输出\stage4b_stable_interface_scatter.png` | 11.8 KB | Stage04B 无序度-能量关联产物。 |  |
| `Stage04_无序度与能量关联\03_输出\stage4b_state_binned_trends.png` | 23.8 KB | Stage04B 无序度-能量关联产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5a_atomwise_physical_descriptors.csv` | 1.44 MB | ★ 常规物理描述符（配位数、键应变、体积、CSP）。最终证明“解释力主要来自这些量而非 SOAP 方向”的就是它。 | 4,659 行 × 23 列 [state, atom_id, x, y, z, interface_id …] |
| `Stage05_物理空间与方向\03_输出\stage5a_descriptor_correlations.csv` | 2.3 KB | Stage05A 物理描述符与空间分布产物。 | 30 行 × 5 列 [scope, descriptor, n, pearson, spearman] |
| `Stage05_物理空间与方向\03_输出\stage5a_descriptor_group_statistics.csv` | 5.2 KB | Stage05A 物理描述符与空间分布产物。 | 40 行 × 8 列 [group, descriptor, n, mean, median, std …] |
| `Stage05_物理空间与方向\03_输出\stage5a_disorder_bin_conditional_stats.csv` | 1.5 KB | Stage05A 物理描述符与空间分布产物。 | 10 行 × 9 列 [bin_id, disorder_low, disorder_high, n, spearman_cn, spearman_rms_bond_strain …] |
| `Stage05_物理空间与方向\03_输出\stage5a_matched_neighborhoods.xyz` | 1.64 MB | Stage05A 物理描述符与空间分布产物。 | 63 原子 · 附加列: central=290 state=stable D10=9.768196778712806 dE=0.3807500000000008 Eres=0.003080000000000638 CN=10.0 RMS=0.052540184592867356 CSP=27.686566695884597 dV=-0.11575476290275333 |
| `Stage05_物理空间与方向\03_输出\stage5a_matched_pair_energy_vs_disorder.png` | 27.2 KB | Stage05A 物理描述符与空间分布产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5a_matched_pairs.csv` | 32.7 KB | Stage05A 物理描述符与空间分布产物。 | 72 行 × 32 列 [category, tolerance, atom_i_state, atom_i_id, atom_j_state, atom_j_id …] |
| `Stage05_物理空间与方向\03_输出\stage5a_metastable_bond_strain_map.png` | 62.3 KB | Stage05A 物理描述符与空间分布产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5a_metastable_disorder_map.png` | 62.5 KB | Stage05A 物理描述符与空间分布产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5a_metastable_physical_descriptors.png` | 56.3 KB | Stage05A 物理描述符与空间分布产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5a_metastable_residual_energy_map.png` | 60.0 KB | Stage05A 物理描述符与空间分布产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5a_residual_definition.csv` | 251 B | Stage05A 物理描述符与空间分布产物。 | 4 行 × 5 列 [primary_bins, sensitivity_bins, spearman_with_primary, pearson_with_primary, n_interface] |
| `Stage05_物理空间与方向\03_输出\stage5a_residual_vs_bond_strain.png` | 27.5 KB | Stage05A 物理描述符与空间分布产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5a_residual_vs_cn.png` | 14.9 KB | Stage05A 物理描述符与空间分布产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5a_residual_vs_csp.png` | 23.9 KB | Stage05A 物理描述符与空间分布产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5a_residual_vs_disorder.png` | 26.7 KB | Stage05A 物理描述符与空间分布产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5a_residual_vs_excess_volume.png` | 25.8 KB | Stage05A 物理描述符与空间分布产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5a_spatial_autocorrelation.csv` | 16.3 KB | Stage05A 物理描述符与空间分布产物。 | 240 行 × 7 列 [state, interface_id, field, r_mid, n_pairs, autocov …] |
| `Stage05_物理空间与方向\03_输出\stage5a_spatial_autocorrelation.png` | 59.2 KB | Stage05A 物理描述符与空间分布产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5a_spatial_patchiness.csv` | 982 B | Stage05A 物理描述符与空间分布产物。 | 8 行 × 8 列 [state, interface_id, set, n, nearest_neighbor_distance, fraction_within_3.0A …] |
| `Stage05_物理空间与方向\03_输出\stage5a_stable_bond_strain_map.png` | 49.5 KB | Stage05A 物理描述符与空间分布产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5a_stable_disorder_map.png` | 46.7 KB | Stage05A 物理描述符与空间分布产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5a_stable_physical_descriptors.png` | 55.1 KB | Stage05A 物理描述符与空间分布产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5a_stable_residual_energy_map.png` | 44.9 KB | Stage05A 物理描述符与空间分布产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5b_atomwise_soap_direction.csv` | 1.41 MB | Stage05B SOAP 方向相关产物。 | 4,659 行 × 25 列 [structure, atom_id, x, y, z, interface_id …] |
| `Stage05_物理空间与方向\03_输出\stage5b_direction_autocorrelation.csv` | 2.4 KB | Stage05B SOAP 方向相关产物。 | 60 行 × 5 列 [state, interface_id, r_mid, n_pairs, mean_cosine] |
| `Stage05_物理空间与方向\03_输出\stage5b_direction_similarity_distributions.png` | 22.0 KB | Stage05B SOAP 方向相关产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5b_direction_vs_disorder.png` | 123.3 KB | Stage05B SOAP 方向相关产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5b_direction_vs_residual_energy.png` | 60.5 KB | Stage05B SOAP 方向相关产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5b_direction_vs_stage5a_descriptors.png` | 49.4 KB | Stage05B SOAP 方向相关产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5b_directional_energy_statistics.csv` | 3.4 KB | Stage05B SOAP 方向相关产物。 | 12 行 × 7 列 [tolerance, category, n_pairs, spearman, pearson, median_dE_by_direction_quintile …] |
| `Stage05_物理空间与方向\03_输出\stage5b_directional_statistics.csv` | 652 B | Stage05B SOAP 方向相关产物。 | 3 行 × 10 列 [category, n_pairs, cosine_mean, cosine_median, cosine_std, cosine_q1 …] |
| `Stage05_物理空间与方向\03_输出\stage5b_matched_energy_difference_vs_direction.png` | 45.6 KB | Stage05B SOAP 方向相关产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5b_metastable_direction_autocorrelation.png` | 25.0 KB | Stage05B SOAP 方向相关产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5b_pca_direction.png` | 33.0 KB | Stage05B SOAP 方向相关产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5b_pca_displacement.png` | 30.9 KB | Stage05B SOAP 方向相关产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5b_pca_displacement_3d.png` | 37.4 KB | Stage05B SOAP 方向相关产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5b_pca_summary.csv` | 1.1 KB | Stage05B SOAP 方向相关产物。 | 20 行 × 4 列 [pca_type, pc, explained_variance_ratio, cumulative] |
| `Stage05_物理空间与方向\03_输出\stage5b_permutation_null.csv` | 227 B | Stage05B SOAP 方向相关产物。 | 1 行 × 8 列 [observed_spearman, null_mean, null_std, null_p5, null_p95, n_permutations …] |
| `Stage05_物理空间与方向\03_输出\stage5b_representative_cases.csv` | 3.9 KB | Stage05B SOAP 方向相关产物。 | 8 行 × 30 列 [category, tolerance, state_i, atom_i, state_j, atom_j …] |
| `Stage05_物理空间与方向\03_输出\stage5b_representative_cases.png` | 35.8 KB | Stage05B SOAP 方向相关产物。 |  |
| `Stage05_物理空间与方向\03_输出\stage5b_representative_direction_cases.xyz` | 117.4 KB | Stage05B SOAP 方向相关产物。 | 59 原子 · 附加列: central=1180 state=stable D=12.14652217671282 dE=-0.0661999999999993 Eres=-0.0015399999999998 CN=13.0 RMS=0.0613492856268347 CSP=9.284815819772238 dV=1.9134011445628063 |
| `Stage05_物理空间与方向\03_输出\stage5b_stable_direction_autocorrelation.png` | 25.0 KB | Stage05B SOAP 方向相关产物。 |  |
| `Stage06_最终评估\03_输出\stage6_correlations.csv` | 24.7 KB | ★ 全项目最终相关系数：原始 0.8317，控制物理描述符后降到 0.2058。这是方法线最重要的产出。 | 无序度 D 与位移范数相关 0.99998 · D 与过剩能相关 0.277 |
| `Stage06_最终评估\03_输出\stage6_direction_by_disorder_bin.csv` | 1.4 KB | Stage06 最终评估产物：相关系数、模型对比与稳健性检验。 | 20 行 × 6 列 [state, bin, n, median_D, mean_within_cosine, median_within_cosine] |
| `Stage06_最终评估\03_输出\stage6_directional_heterogeneity.png` | 21.7 KB | Stage06 最终评估产物：相关系数、模型对比与稳健性检验。 |  |
| `Stage06_最终评估\03_输出\stage6_energy_model_improvement.png` | 20.5 KB | Stage06 最终评估产物：相关系数、模型对比与稳健性检验。 |  |
| `Stage06_最终评估\03_输出\stage6_energy_models.csv` | 732 B | ★ 能量模型 R² 对比：仅 D 为 0.0393，加入物理描述符后 0.9869。量化了 SOAP 坐标的增量信息有限。 | Model0_D: R²=0.0393 · Model1_D_phys: R²=0.9869 · Model2_D_phys_dir: R²=0.9886 · ModelD_D_dir: R²=0.7497 |
| `Stage06_最终评估\03_输出\stage6_figure_summary.csv` | 536 B | Stage06 最终评估产物：相关系数、模型对比与稳健性检验。 | 6 行 × 2 列 [figure, content] |
| `Stage06_最终评估\03_输出\stage6_magnitude_matched_direction_energy.png` | 45.2 KB | Stage06 最终评估产物：相关系数、模型对比与稳健性检验。 |  |
| `Stage06_最终评估\03_输出\stage6_matched_analysis.csv` | 3.4 KB | Stage06 最终评估产物：相关系数、模型对比与稳健性检验。 | 16 行 × 12 列 [tolerance, category, n_pairs, spearman, pearson, bootstrap_ci_low …] |
| `Stage06_最终评估\03_输出\stage6_representative_cases.png` | 35.7 KB | Stage06 最终评估产物：相关系数、模型对比与稳健性检验。 |  |
| `Stage06_最终评估\03_输出\stage6_representative_directional_cases.csv` | 7.8 KB | Stage06 最终评估产物：相关系数、模型对比与稳健性检验。 | 12 行 × 42 列 [category, tolerance, state_i, atom_i, state_j, atom_j …] |
| `Stage06_最终评估\03_输出\stage6_robustness.csv` | 211 B | Stage06 最终评估产物：相关系数、模型对比与稳健性检验。 | 4 行 × 3 列 [condition, n_pairs, spearman] |
| `Stage06_最终评估\03_输出\stage6_scalar_disorder_vs_energy.png` | 28.5 KB | Stage06 最终评估产物：相关系数、模型对比与稳健性检验。 |  |
| `Stage06_最终评估\03_输出\stage6_spatial_direction_residual.png` | 80.3 KB | Stage06 最终评估产物：相关系数、模型对比与稳健性检验。 |  |
| `Stage06_最终评估\03_输出\stage6_spatial_summary.csv` | 955 B | Stage06 最终评估产物：相关系数、模型对比与稳健性检验。 | 8 行 × 9 列 [state, interface_id, field, n, mean, std …] |
| `Stage06_最终评估\03_输出\stage6_state_classification.csv` | 1.1 KB | Stage06 最终评估产物：相关系数、模型对比与稳健性检验。 | 10 行 × 8 列 [test_interface, n_train, n_test, roc_auc, balanced_accuracy, f1 …] |

### 代码（27 个，468.2 KB）

| 文件 | 大小 | 用途 |
|---|---|---|
| `06_环境与实现\setup_soap_env.sh` | 2.5 KB | SOAP 分析环境的一键配置脚本。 |
| `06_环境与实现\soap_analysis\__init__.py` | 1.8 KB | Python 脚本。 |
| `06_环境与实现\soap_analysis\bulk_reference.py` | 1.8 KB | Python 脚本。 |
| `06_环境与实现\soap_analysis\compare_structures.py` | 3.1 KB | Python 脚本。 |
| `06_环境与实现\soap_analysis\dedup.py` | 2.8 KB | Python 脚本。 |
| `06_环境与实现\soap_analysis\disorder_profile.py` | 5.1 KB | Python 脚本。 |
| `06_环境与实现\soap_analysis\io_lammps.py` | 5.7 KB | Python 脚本。 |
| `06_环境与实现\soap_analysis\run_analysis.py` | 11.1 KB | Python 脚本。 |
| `06_环境与实现\soap_analysis\similarity.py` | 2.1 KB | Python 脚本。 |
| `06_环境与实现\soap_analysis\soap_backend.py` | 12.3 KB | Python 脚本。 |
| `06_环境与实现\verify_soap_env.py` | 5.3 KB | 环境自检脚本，确认 SOAP 依赖可用。 |
| `90_历史归档\旧版材料\inspect_structures_v_old.py` | 3.9 KB | Python 脚本。 |
| `Stage01_数据准备与SOAP初测\02_脚本\cfg_reader.py` | 8.7 KB | CFG 格式解析脚本。 |
| `Stage01_数据准备与SOAP初测\02_脚本\inspect_structures.py` | 4.2 KB | 结构概览脚本，检查原子数与字段完整性。 |
| `Stage01_数据准备与SOAP初测\02_脚本\soap_test.py` | 5.0 KB | SOAP 初次可计算性测试脚本。 |
| `Stage01_数据准备与SOAP初测\02_脚本\stage1_qc.py` | 21.2 KB | Stage01 数据质量检查脚本。 |
| `Stage02_SOAP表示与审计\02_脚本\stage2_5_audit.py` | 26.4 KB | Stage02.5 表示审计产物。 |
| `Stage02_SOAP表示与审计\02_脚本\stage2_75_full.py` | 31.4 KB | Stage02.75 全人口 SOAP 与距离审计产物。 |
| `Stage02_SOAP表示与审计\02_脚本\stage2_profiles.py` | 8.9 KB | Stage02 SOAP 表示产物。 |
| `Stage02_SOAP表示与审计\02_脚本\stage2_soap.py` | 14.1 KB | Stage02 SOAP 表示产物。 |
| `Stage03_LAE阈值与LAE_LER构建\02_脚本\stage3a_lae_threshold.py` | 56.6 KB | Stage03A 阈值扫描产物：衡量 epsilon 取值对结果的影响。 |
| `Stage03_LAE阈值与LAE_LER构建\02_脚本\stage3b_lae_ler.py` | 39.4 KB | Stage03B 类别构建产物：LAE 分类与 LER 人口指纹。 |
| `Stage04_无序度与能量关联\02_脚本\stage4a_continuous_disorder.py` | 27.6 KB | Stage04A 连续无序度产物。 |
| `Stage04_无序度与能量关联\02_脚本\stage4b_disorder_energy.py` | 33.2 KB | Stage04B 无序度-能量关联产物。 |
| `Stage05_物理空间与方向\02_脚本\stage5a_physical_spatial_residual.py` | 49.9 KB | Stage05A 物理描述符与空间分布产物。 |
| `Stage05_物理空间与方向\02_脚本\stage5b_soap_direction.py` | 46.7 KB | Stage05B SOAP 方向相关产物。 |
| `Stage06_最终评估\02_脚本\stage6_final_analysis.py` | 37.4 KB | Stage06 最终评估产物：相关系数、模型对比与稳健性检验。 |

### 输入（3 个，8.03 MB）

| 文件 | 大小 | 用途 | 内容分析 |
|---|---|---|---|
| `06_环境与实现\requirements_soap.txt` | 1.5 KB | SOAP 分析的 Python 依赖清单。 |  |
| `Stage01_数据准备与SOAP初测\02_脚本\ni_metastable.cfg` | 4.01 MB | 输入构型之一：亚稳的 Ni 晶界构型（源数据）。 | 47,120 原子 |
| `Stage01_数据准备与SOAP初测\02_脚本\ni_stable.cfg` | 4.02 MB | 输入构型之一：稳定的 Ni 晶界构型（源数据）。 | 47,120 原子 |

### 文档（10 个，35.3 KB）

| 文件 | 大小 | 用途 |
|---|---|---|
| `00_汇总报告\Ni_SOAP_最终报告.md` | 6.3 KB | ★ 方法探索总报告：项目定位、六个 Stage 的结论、十条下一版建议、以及“不能确定”清单。 |
| `00_汇总报告\项目文件索引_脚本输出报告_按阶段.md` | 3.9 KB | Markdown 文档。 |
| `90_历史归档\归档说明.md` | 1.8 KB | Markdown 文档。 |
| `README.md` | 2.0 KB | 仓库总入口：项目组成、核心结论、文档索引、环境与日常操作。 |
| `Stage01_数据准备与SOAP初测\01_报告\Stage01_确定性结论与待讨论.md` | 3.2 KB | Markdown 文档。 |
| `Stage02_SOAP表示与审计\01_报告\Stage02_确定性结论与待讨论.md` | 4.0 KB | Markdown 文档。 |
| `Stage03_LAE阈值与LAE_LER构建\01_报告\Stage03_确定性结论与待讨论.md` | 3.1 KB | Markdown 文档。 |
| `Stage04_无序度与能量关联\01_报告\Stage04_确定性结论与待讨论.md` | 3.7 KB | Markdown 文档。 |
| `Stage05_物理空间与方向\01_报告\Stage05_确定性结论与待讨论.md` | 3.8 KB | Markdown 文档。 |
| `Stage06_最终评估\01_报告\Stage06_确定性结论与待讨论.md` | 3.5 KB | Markdown 文档。 |

## 归档区：为溯源保留的历史副本

`_archive` — 156 个文件，3.45 MB

### 结论（10 个，2.20 MB）

| 文件 | 大小 | 用途 | 内容分析 |
|---|---|---|---|
| `GB_project_20260819\archive\legacy\li_gb_project\lammps_test\final.data` | 728 B | LAMMPS 原子构型：某个阶段结束时的静态结构快照。 | 16 原子 · 盒子 7.02 × 7.02 × 7.02 Å · timestep 1 |
| `GB_project_20260819\archive\old_versions\legacy_generator\lammps_workflow\Sigma5_310.lmp` | 261.6 KB | LAMMPS 结构文件。 | 3,224 原子 · 盒子 85.98 × 60.79 × 14.89 Å |
| `GB_project_20260819\archive\old_versions\legacy_generator\lammps_workflow\Sigma5_310_minimized.data` | 248.1 KB | LAMMPS 原子构型：某个阶段结束时的静态结构快照。 | 3,224 原子 · 盒子 85.98 × 61.20 × 15.30 Å · timestep 277 |
| `GB_project_20260819\archive\old_versions\legacy_generator\lammps_workflow\annealed.data` | 408.6 KB | LAMMPS 原子构型：某个阶段结束时的静态结构快照。 | 3,224 原子 · 盒子 85.98 × 61.20 × 15.30 Å · timestep 20,277 |
| `GB_project_20260819\archive\old_versions\legacy_generator\lammps_workflow\final_gb.data` | 414.4 KB | LAMMPS 原子构型：某个阶段结束时的静态结构快照。 | 3,224 原子 · 盒子 85.98 × 61.20 × 15.30 Å · timestep 50,044 |
| `GB_project_20260819\archive\old_versions\legacy_generator\output\Sigma5_310.cif` | 280.1 KB | 晶体学标准结构文件。 |  |
| `GB_project_20260819\archive\old_versions\legacy_generator\output\Sigma5_310.lmp` | 261.6 KB | LAMMPS 结构文件。 | 3,224 原子 · 盒子 85.98 × 60.79 × 14.89 Å |
| `GB_project_20260819\archive\old_versions\legacy_generator\output\Sigma5_310.xyz` | 226.7 KB | 扩展 XYZ 结构文件，可携带每原子附加列。 | 3,224 原子 |
| `GB_project_20260819\archive\old_versions\legacy_generator\output\Sigma5_310_grainID.xyz` | 147.3 KB | 扩展 XYZ 结构文件，可携带每原子附加列。 | 3,224 原子 · 附加列: Lattice="-35.1000000000 35.1000000000 70.2000000000 42.1200000000 -42.1200000000 21.0600000000 10.5300000000 10.5300000000 -10.5300000000" Properties=species:S:1:pos:R:3:grain_id:R:1 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\DO_NOT_USE.txt` | 179 B | 纯文本结果或检查报告。 |  |

### 续算（2 个，556.0 KB）

| 文件 | 大小 | 用途 |
|---|---|---|
| `GB_project_20260819\archive\old_versions\legacy_generator\lammps_workflow\Sigma5_310_minimized.restart` | 278.0 KB | LAMMPS 续算检查点：保存速度场与热力学状态，只能在 LAMMPS 内读取。 |
| `GB_project_20260819\archive\old_versions\legacy_generator\lammps_workflow\final_gb.restart` | 278.0 KB | LAMMPS 续算检查点：保存速度场与热力学状态，只能在 LAMMPS 内读取。 |

### 追溯（5 个，272.2 KB）

| 文件 | 大小 | 用途 | 内容分析 |
|---|---|---|---|
| `GB_project_20260819\_unique_preserved\early_env_test\kim.log` | 6.2 KB | OpenKIM 势函数调用日志，记录实际加载的模型标识。 |  |
| `GB_project_20260819\archive\old_versions\legacy_generator\lammps_workflow\kim.log` | 240.6 KB | OpenKIM 势函数调用日志，记录实际加载的模型标识。 |  |
| `GB_project_20260819\archive\old_versions\legacy_generator\lammps_workflow\log.lammps` | 17.9 KB | LAMMPS 运行日志，含各阶段数值历史与断言结果。 |  |
| `Li_metastable_2x2_legacy_20260819\legacy\reports_v2\ENVIRONMENT.txt` | 180 B | 纯文本结果或检查报告。 |  |
| `Li_metastable_2x2_legacy_20260819\legacy\reports_v2\SHA256SUMS.txt` | 7.3 KB | 纯文本结果或检查报告。 | 收录 80 个文件的校验和 |

### 代码（62 个，201.9 KB）

| 文件 | 大小 | 用途 |
|---|---|---|
| `GB_project_20260819\archive\legacy\li_gb_project\lammps_test\run.sh` | 563 B | Shell 脚本（在 WSL/Linux 下运行）。 |
| `GB_project_20260819\archive\legacy\li_gb_project\setup_env.sh` | 1.4 KB | Shell 脚本（在 WSL/Linux 下运行）。 |
| `GB_project_20260819\archive\old_versions\legacy_generator\__init__.py` | 424 B | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\legacy_generator\__main__.py` | 97 B | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\legacy_generator\builder.py` | 9.3 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\legacy_generator\config.py` | 4.7 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\legacy_generator\examples\example.py` | 1.2 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\legacy_generator\exporter.py` | 5.0 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\legacy_generator\geometry.py` | 4.0 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\legacy_generator\lammps_workflow\run_relax.sh` | 1.7 KB | Shell 脚本（在 WSL/Linux 下运行）。 |
| `GB_project_20260819\archive\old_versions\legacy_generator\main.py` | 4.4 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\legacy_generator\utils.py` | 2.9 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\generator\__init__.py` | 86 B | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\generator\builder.py` | 8.6 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\generator\bulk.py` | 2.3 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\generator\config.py` | 2.1 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\generator\exporters.py` | 5.6 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\generator\generate.py` | 9.9 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\generator\geometry.py` | 9.6 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\generator\lammps_io.py` | 3.4 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\generator\reimport_check.py` | 4.2 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\generator\verify.py` | 15.3 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\lammps\resume.sh` | 787 B | 从中断处续算的脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\lammps\run_all.sh` | 1.4 KB | 完整热历史编排脚本；带 .done 标记实现断点续算。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\tests\conftest.py` | 223 B | 测试公共夹具。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\tests\test_builder.py` | 2.4 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\tests\test_conventional_bcc.py` | 1.1 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\tests\test_grain_id.py` | 1.3 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\tests\test_legacy_regression.py` | 683 B | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\tests\test_no_nvt_zero.py` | 922 B | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\tests\test_outputs_consistency.py` | 1.5 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\tests\test_overlap_removal.py` | 1.2 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\tests\test_periodicity.py` | 1.4 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\tests\test_sigma5_geometry.py` | 1.7 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\tests\test_velocities_zero.py` | 1.2 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\tests\test_verify.py` | 1.8 KB | Python 脚本。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\tests\test_workflow_scripts.py` | 1.7 KB | Python 脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\generator_old\__init__.py` | 86 B | Python 脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\generator_old\builder.py` | 8.8 KB | Python 脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\generator_old\bulk.py` | 2.3 KB | Python 脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\generator_old\config.py` | 2.1 KB | Python 脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\generator_old\exporters.py` | 5.6 KB | Python 脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\generator_old\generate.py` | 10.1 KB | Python 脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\generator_old\geometry.py` | 9.6 KB | Python 脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\generator_old\lammps_io.py` | 3.4 KB | Python 脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\generator_old\reimport_check.py` | 4.2 KB | Python 脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\generator_old\verify.py` | 16.7 KB | Python 脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\lammps_v2\resume.sh` | 2.0 KB | 从中断处续算的脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\lammps_v2\run_all.sh` | 1.3 KB | 完整热历史编排脚本；带 .done 标记实现断点续算。 |
| `Li_metastable_2x2_legacy_20260819\legacy\tests_v2\conftest.py` | 223 B | 测试公共夹具。 |
| `Li_metastable_2x2_legacy_20260819\legacy\tests_v2\test_builder.py` | 2.5 KB | Python 脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\tests_v2\test_conventional_bcc.py` | 1.1 KB | Python 脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\tests_v2\test_grain_id.py` | 1.3 KB | Python 脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\tests_v2\test_legacy_regression.py` | 1.0 KB | Python 脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\tests_v2\test_no_nvt_zero.py` | 922 B | Python 脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\tests_v2\test_outputs_consistency.py` | 1.5 KB | Python 脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\tests_v2\test_overlap_removal.py` | 1.2 KB | Python 脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\tests_v2\test_periodicity.py` | 1.4 KB | Python 脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\tests_v2\test_sigma5_geometry.py` | 1.7 KB | Python 脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\tests_v2\test_velocities_zero.py` | 1.2 KB | Python 脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\tests_v2\test_verify.py` | 2.5 KB | Python 脚本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\tests_v2\test_workflow_scripts.py` | 3.7 KB | Python 脚本。 |

### 输入（25 个，39.4 KB）

| 文件 | 大小 | 用途 |
|---|---|---|
| `GB_project_20260819\archive\legacy\li_gb_project\.gitignore` | 210 B | Git 忽略规则，指定哪些文件不纳入版本控制。 |
| `GB_project_20260819\archive\legacy\li_gb_project\lammps_test\in.test_li` | 1.4 KB | 早期环境连通性测试：确认 LAMMPS + OpenKIM 能正确调用 Li 势函数。 |
| `GB_project_20260819\archive\legacy\li_gb_project\requirements.txt` | 1.1 KB | Python 依赖清单，用于重建运行环境。 |
| `GB_project_20260819\archive\old_versions\legacy_generator\lammps_workflow\in.bicrystal_relax` | 9.2 KB | 早期版本的双晶弛豫脚本（已被 v2 工作流取代）。 |
| `GB_project_20260819\archive\old_versions\legacy_generator\lammps_workflow\in.quench` | 804 B | 早期版本的淬火脚本。 |
| `GB_project_20260819\archive\old_versions\legacy_generator\lammps_workflow\in.resume` | 9.2 KB | 早期版本的续算脚本。注意：旧版此文件曾直接复制主脚本，是已废弃流程中的一个已知缺陷。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\lammps\in.0_check` | 659 B | Stage 0 输入脚本：核对模型几何、原子数与势函数是否被正确加载。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\lammps\in.1_minimize` | 1.3 KB | Stage 1 输入脚本：0 K 极小化，仅沿晶界法向做零压弛豫（y/z 固定在零压 CSL 周期），消除初始残余力。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\lammps\in.2_anneal` | 939 B | Stage 2 输入脚本：363 K 等压（NPT）退火，是诱发晶界热重构的关键一步。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\lammps\in.3_quench` | 893 B | Stage 3 输入脚本：以 10 K/ps 从 363 K 淬火到 10 K，把高温无序结构冻结下来。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\lammps\in.4_final` | 1.1 KB | Stage 4 输入脚本：最终弛豫与保持，产出用于结构分析的末态。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\lammps\in.bulk_energy` | 997 B | 体相能量参考计算：给出每原子势能基准，用于扣除体相贡献、得到晶界过剩量。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\lammps\in.dump2data` | 725 B | dump→data 重建脚本：把轨迹转回结构文件，用来验证轨迹未丢失信息。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\lammps\in.resume` | 1.7 KB | 早期版本的续算脚本。注意：旧版此文件曾直接复制主脚本，是已废弃流程中的一个已知缺陷。 |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\requirements.txt` | 34 B | Python 依赖清单，用于重建运行环境。 |
| `Li_metastable_2x2_legacy_20260819\legacy\lammps_v2\in.bulk_energy` | 990 B | 体相能量参考计算：给出每原子势能基准，用于扣除体相贡献、得到晶界过剩量。 |
| `Li_metastable_2x2_legacy_20260819\legacy\lammps_v2\in.dump2data` | 546 B | dump→data 重建脚本：把轨迹转回结构文件，用来验证轨迹未丢失信息。 |
| `Li_metastable_2x2_legacy_20260819\legacy\lammps_v2\stage0_test\in.0_check` | 662 B | Stage 0 输入脚本：核对模型几何、原子数与势函数是否被正确加载。 |
| `Li_metastable_2x2_legacy_20260819\legacy\lammps_v2\stage1_minimize\in.1_minimize` | 1.5 KB | Stage 1 输入脚本：0 K 极小化，仅沿晶界法向做零压弛豫（y/z 固定在零压 CSL 周期），消除初始残余力。 |
| `Li_metastable_2x2_legacy_20260819\legacy\lammps_v2\stage2_anneal\in.2_anneal` | 953 B | Stage 2 输入脚本：363 K 等压（NPT）退火，是诱发晶界热重构的关键一步。 |
| `Li_metastable_2x2_legacy_20260819\legacy\lammps_v2\stage2_anneal\in.2_anneal_resume` | 1.1 KB | Stage 2 的续算版本：从中断的检查点继续退火，不重跑已完成的时段。 |
| `Li_metastable_2x2_legacy_20260819\legacy\lammps_v2\stage3_quench\in.3_quench` | 964 B | Stage 3 输入脚本：以 10 K/ps 从 363 K 淬火到 10 K，把高温无序结构冻结下来。 |
| `Li_metastable_2x2_legacy_20260819\legacy\lammps_v2\stage3_quench\in.3_quench_resume` | 1.0 KB | Stage 3 的续算版本。 |
| `Li_metastable_2x2_legacy_20260819\legacy\lammps_v2\stage4_relax\in.4_final` | 1.6 KB | Stage 4 输入脚本：最终弛豫与保持，产出用于结构分析的末态。 |
| `Li_metastable_2x2_legacy_20260819\legacy\requirements_v2.txt` | 34 B | 纯文本结果或检查报告。 |

### 文档（52 个，217.7 KB）

| 文件 | 大小 | 用途 | 内容分析 |
|---|---|---|---|
| `GB_project_20260819\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\docs\README.md` | 10.9 KB | 仓库总入口：项目组成、核心结论、文档索引、环境与日常操作。 |  |
| `GB_project_20260819\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\docs\notes\README.md` | 18 B | 仓库总入口：项目组成、核心结论、文档索引、环境与日常操作。 |  |
| `GB_project_20260819\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\docs\reports\FINAL_MIGRATION_PLAN.md` | 4.7 KB | Markdown 文档。 |  |
| `GB_project_20260819\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\docs\reports\FINAL_STATUS.md` | 426 B | Markdown 文档。 |  |
| `GB_project_20260819\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\docs\reports\PROJECT_CLEANUP_REPORT.md` | 5.2 KB | Markdown 文档。 |  |
| `GB_project_20260819\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\docs\reports\build_report.md` | 2.1 KB | Markdown 文档。 |  |
| `GB_project_20260819\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\docs\reports\conflict_report.md` | 5.3 KB | Markdown 文档。 |  |
| `GB_project_20260819\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\docs\reports\file_inventory.md` | 10.7 KB | Markdown 文档。 |  |
| `GB_project_20260819\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\docs\reports\migration_inventory.md` | 11.5 KB | Markdown 文档。 |  |
| `GB_project_20260819\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production\docs\reports\verification_report.md` | 1.8 KB | Markdown 文档。 |  |
| `GB_project_20260819\02_MD_Simulation\README.md` | 18 B | 仓库总入口：项目组成、核心结论、文档索引、环境与日常操作。 |  |
| `GB_project_20260819\03_Analysis\README.md` | 18 B | 仓库总入口：项目组成、核心结论、文档索引、环境与日常操作。 |  |
| `GB_project_20260819\04_Potential_Validation\README.md` | 18 B | 仓库总入口：项目组成、核心结论、文档索引、环境与日常操作。 |  |
| `GB_project_20260819\05_DFT_Validation\README.md` | 18 B | 仓库总入口：项目组成、核心结论、文档索引、环境与日常操作。 |  |
| `GB_project_20260819\06_Publication\README.md` | 18 B | 仓库总入口：项目组成、核心结论、文档索引、环境与日常操作。 |  |
| `GB_project_20260819\07_Documentation\PROJECT_MIGRATION_REPORT.md` | 2.1 KB | Markdown 文档。 |  |
| `GB_project_20260819\README.md` | 792 B | 仓库总入口：项目组成、核心结论、文档索引、环境与日常操作。 |  |
| `GB_project_20260819\archive\debugging\cs5_310_ref_result.json` | 33.2 KB | 参数与元数据。 | 字段: schema, description, parameters, orientation, counts, distances, gb_equivalence, bulk_bcc |
| `GB_project_20260819\archive\legacy\li_gb_project\README.md` | 5.4 KB | 仓库总入口：项目组成、核心结论、文档索引、环境与日常操作。 |  |
| `GB_project_20260819\archive\legacy\li_gb_project\SETUP.md` | 2.7 KB | Markdown 文档。 |  |
| `GB_project_20260819\archive\legacy\li_gb_project\lammps_test\README.md` | 3.2 KB | 仓库总入口：项目组成、核心结论、文档索引、环境与日常操作。 |  |
| `GB_project_20260819\archive\not_migrated\NOTES.md` | 627 B | Markdown 文档。 |  |
| `GB_project_20260819\archive\old_versions\legacy_generator\README.md` | 6.0 KB | 仓库总入口：项目组成、核心结论、文档索引、环境与日常操作。 |  |
| `GB_project_20260819\archive\old_versions\legacy_generator\lammps_workflow\anneal_save.in` | 365 B | 早期退火脚本，含中途保存设置。 |  |
| `GB_project_20260819\archive\old_versions\legacy_generator\lammps_workflow\continue.in` | 609 B | 早期续算脚本。 |  |
| `GB_project_20260819\archive\old_versions\legacy_generator\lammps_workflow\dump2data.in` | 492 B | 早期 dump→data 转换脚本。 |  |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\README.md` | 7.5 KB | 仓库总入口：项目组成、核心结论、文档索引、环境与日常操作。 |  |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\pytest.ini` | 40 B | 测试框架配置。 |  |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\reports\build_report.md` | 4.6 KB | Markdown 文档。 |  |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\reports\reference_energy_bulk.md` | 1.2 KB | Markdown 文档。 |  |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\reports\stage_io_table.md` | 2.1 KB | Markdown 文档。 |  |
| `GB_project_20260819\archive\old_versions\v2_obsolete_20260802_150439\reports\verification_report.md` | 2.6 KB | Markdown 文档。 |  |
| `Li_metastable_2x2_legacy_20260819\PROMPT_FOR_CODEX.md` | 19.1 KB | Markdown 文档。 |  |
| `Li_metastable_2x2_legacy_20260819\PROTOCOL.md` | 5.0 KB | Markdown 文档。 |  |
| `Li_metastable_2x2_legacy_20260819\PROVENANCE.md` | 2.2 KB | Markdown 文档。 |  |
| `Li_metastable_2x2_legacy_20260819\README.md` | 2.7 KB | 仓库总入口：项目组成、核心结论、文档索引、环境与日常操作。 |  |
| `Li_metastable_2x2_legacy_20260819\RUN_STATUS.md` | 2.3 KB | Markdown 文档。 |  |
| `Li_metastable_2x2_legacy_20260819\legacy\MIGRATION_MAP.md` | 1.1 KB | Markdown 文档。 |  |
| `Li_metastable_2x2_legacy_20260819\legacy\README_LEGACY_CORE.md` | 712 B | Markdown 文档。 |  |
| `Li_metastable_2x2_legacy_20260819\legacy\README_v2_final.md` | 10.9 KB | Markdown 文档。 |  |
| `Li_metastable_2x2_legacy_20260819\legacy\reports_v2\FINAL_MIGRATION_PLAN.md` | 4.7 KB | Markdown 文档。 |  |
| `Li_metastable_2x2_legacy_20260819\legacy\reports_v2\FINAL_STATUS.md` | 426 B | Markdown 文档。 |  |
| `Li_metastable_2x2_legacy_20260819\legacy\reports_v2\PROJECT_CLEANUP_REPORT.md` | 5.2 KB | Markdown 文档。 |  |
| `Li_metastable_2x2_legacy_20260819\legacy\reports_v2\build_report.md` | 2.1 KB | Markdown 文档。 |  |
| `Li_metastable_2x2_legacy_20260819\legacy\reports_v2\conflict_report.md` | 5.3 KB | Markdown 文档。 |  |
| `Li_metastable_2x2_legacy_20260819\legacy\reports_v2\file_inventory.md` | 10.7 KB | Markdown 文档。 |  |
| `Li_metastable_2x2_legacy_20260819\legacy\reports_v2\migration_inventory.md` | 11.5 KB | Markdown 文档。 |  |
| `Li_metastable_2x2_legacy_20260819\legacy\reports_v2\verification_report.md` | 1.8 KB | Markdown 文档。 |  |
| `Li_metastable_2x2_legacy_20260819\legacy\tests_v2\pytest.ini` | 40 B | 测试框架配置。 |  |
| `Li_metastable_2x2_legacy_20260819\stage0_build_validate\build_report.md` | 2.3 KB | Markdown 文档。 |  |
| `Li_metastable_2x2_legacy_20260819\stage0_build_validate\validation_report.md` | 1.8 KB | Markdown 文档。 |  |
| `Li_metastable_2x2_legacy_20260819\stage0_build_validate\verification_report.md` | 1.8 KB | Markdown 文档。 |  |

## 项目根与工具

仓库根目录与 `tools/` — 13 个文件，127.6 KB

### 结论（1 个，16.8 KB）

| 文件 | 大小 | 用途 | 内容分析 |
|---|---|---|---|
| `DATA_MANIFEST.csv` | 16.8 KB | 数据表。 | 60 行 × 7 列 [﻿相对路径, 文件名, 扩展名, 层次, 大小(字节), 用途说明 …] |

### 代码（4 个，67.3 KB）

| 文件 | 大小 | 用途 |
|---|---|---|
| `tools\artifact_ledger.py` | 28.0 KB | Python 脚本。 |
| `tools\build_artifact_catalog.py` | 30.3 KB | Python 脚本。 |
| `tools\compute_gb_energy.py` | 6.3 KB | Python 脚本。 |
| `tools\install_git_hook.ps1` | 2.7 KB | PowerShell 脚本。 |

### 输入（2 个，2.1 KB）

| 文件 | 大小 | 用途 |
|---|---|---|
| `.gitignore` | 1.5 KB | Git 忽略规则，指定哪些文件不纳入版本控制。 |
| `refresh_ledger.cmd` | 657 B | （未分类） |

### 文档（6 个，41.4 KB）

| 文件 | 大小 | 用途 |
|---|---|---|
| `.gitattributes` | 748 B | （未分类） |
| `DATA_MANIFEST.md` | 12.4 KB | Markdown 文档。 |
| `PROJECT_REVIEW.md` | 10.6 KB | Markdown 文档。 |
| `README.md` | 4.1 KB | 仓库总入口：项目组成、核心结论、文档索引、环境与日常操作。 |
| `STUDY_GUIDE.md` | 11.0 KB | Markdown 文档。 |
| `tools\README.md` | 2.7 KB | 仓库总入口：项目组成、核心结论、文档索引、环境与日常操作。 |
