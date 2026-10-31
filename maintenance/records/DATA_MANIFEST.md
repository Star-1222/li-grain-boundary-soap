# 未纳入版本控制的数据清单

> 本文件由 `maintenance/tools/artifact_ledger.py --data-manifest` 生成。
> `ARTIFACT_LEDGER.md` 只覆盖仓库内容；本文件登记的是被 `.gitignore`
> 排除在仓库之外、但**仍然重要**的本地数据。
> 全新克隆不含这些文件 —— 本清单即为还原它们时的核对表。

- 条目数：**60**
- 合计体积：**466.24 MB**

## 为什么不在仓库里

| 层次 | 文件数 | 体积 | 原因 |
|---|---:|---:|---|
| **中间** | 14 | 283.67 MB | 可由上游输入与脚本重算，代价是算力 |
| **追溯** | 1 | 8.8 KB | 运行日志，体积偏大 |
| **文献** | 43 | 182.49 MB | 第三方论文，体积大且有版权 |
| **文档** | 2 | 70.6 KB | 文档 |

## 清单

| 文件 | 大小 | 层次 | SHA256（前 16 位） | 是什么 |
|---|---:|---|---|---|
| `soap\Stage05_物理空间与方向\03_输出\stage5b_matched_direction_energy_pairs.csv` | 142.87 MB | 中间 | `9c286d7b1a370e65` | ★ 配对重采样中间表（142.9 MB，占本项目近一半体积）。最终结论只用到由它算出的相关系数，本体属可再生的中间产物。 |
| `soap\Stage02_SOAP表示与审计\03_输出\stage2_75_pairwise_sample.csv` | 48.78 MB | 中间 | `f71657f1f6fa104a` | 原子对距离样本表（48.8 MB）。支撑“原始距离里 42–46% 来自范数差”这一结论，属中间表。 |
| `soap\Stage02_SOAP表示与审计\03_输出\stage2_75_soap_vectors.npy` | 36.04 MB | 中间 | `bfa4378375dd7e9e` | ★ SOAP 主矩阵（4659 原子 × 1014 维，36 MB）。后续全部分析都建立在这个矩阵上；重算耗时较长。 |
| `grain_boundary\00_Literature\soap\2018_Priedeman_QuantifyingConnectingGBStructure_长版54页.pdf` | 17.42 MB | 文献 | `ccd4e8caf1b5872b` | PDF 文档。 |
| `soap\05_文献\2018_Priedeman_QuantifyingConnectingGBStructure_长版54页.pdf` | 17.42 MB | 文献 | `ccd4e8caf1b5872b` | PDF 文档。 |
| `grain_boundary\history\experiments\phase1_2026-09\Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\4_淬火10K\p1_quench.dump` | 15.30 MB | 中间 | `5e6c0aefe31bbdee` | 淬火过程轨迹（35.3 ps）。用于核实降温速率是否严格为 10 K/ps。 |
| `grain_boundary\history\experiments\phase1_2026-09\Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\2_退火363K零压\p1_anneal.dump` | 11.05 MB | 中间 | `267d2a85dd8aa74c` | 退火 25 ps 的过程轨迹。报告里“末 16 ps 时均 T=362.8 K”就是从这里统计出来的。 |
| `grain_boundary\00_Literature\soap\2024_Mishra_PRM_LearningFromMetastableSymmetricTiltGBs_PhysicsBasedDescriptors.pdf` | 10.41 MB | 文献 | `aff947cec0afc3ac` | PDF 文档。 |
| `soap\05_文献\2024_Mishra_PRM_LearningFromMetastableSymmetricTiltGBs_PhysicsBasedDescriptors.pdf` | 10.41 MB | 文献 | `aff947cec0afc3ac` | PDF 文档。 |
| `grain_boundary\history\construction_baseline_3600\Sigma5_310_Metastable_GB\v2_production\data\quenched\stage3_quench.dump` | 9.89 MB | 中间 | `179a17eebc4cf6c9` | 淬火过程轨迹。用于核对降温速率与系统响应是否正常。 |
| `grain_boundary\00_Literature\soap\2018_Priedeman_ActaMater_QuantifyingConnectingAtomicCrystallographicGBStructure.pdf` | 8.02 MB | 文献 | `ac036446151255f1` | PDF 文档。 |
| `soap\05_文献\2018_Priedeman_ActaMater_QuantifyingConnectingAtomicCrystallographicGBStructure.pdf` | 8.02 MB | 文献 | `ac036446151255f1` | PDF 文档。 |
| `grain_boundary\history\construction_baseline_3600\Sigma5_310_Metastable_GB\v2_production\data\annealed\stage2_anneal.dump` | 7.90 MB | 中间 | `18d72049a0bc68b8` | 退火过程轨迹。用来观察晶界重构是怎样逐步发生的，而不只是看首末两帧。 |
| `grain_boundary\00_Literature\soap\2020_Fujii_NatCommun_GBThermalConductivity_LocalAtomicEnvironments_作者版.pdf` | 7.00 MB | 文献 | `62172a7dcfcf808c` | PDF 文档。 |
| `soap\05_文献\2020_Fujii_NatCommun_GBThermalConductivity_LocalAtomicEnvironments_作者版.pdf` | 7.00 MB | 文献 | `62172a7dcfcf808c` | PDF 文档。 |
| `grain_boundary\00_Literature\soap\2013_Bartok_PRB_OnRepresentingChemicalEnvironments_SOAP原始论文.pdf` | 6.61 MB | 文献 | `51c56e5bb24ec03d` | PDF 文档。 |
| `soap\05_文献\2013_Bartok_PRB_OnRepresentingChemicalEnvironments_SOAP原始论文.pdf` | 6.61 MB | 文献 | `51c56e5bb24ec03d` | PDF 文档。 |
| `grain_boundary\00_Literature\lithium_metal\mechanistic-insights-into-li-transport-enabled-by-isolated-sulfur-species-in-li3ps4-glasses.pdf` | 6.36 MB | 文献 | `fa9ce638d7d76172` | PDF 文档。 |
| `grain_boundary\00_Literature\grain_boundary\2024 --Metastable grain boundaries the roles of structural and chemical disorders in their energetics, non-equilibrium kinetic evolution, and mechanical behaviors.pdf` | 5.46 MB | 文献 | `ce4cba3d2a8062f7` | PDF 文档。 |
| `grain_boundary\00_Literature\soap\2020_Fujii_NatCommun_GBThermalConductivity_LocalAtomicEnvironments_正式版.pdf` | 4.96 MB | 文献 | `355cba548cd681e2` | PDF 文档。 |
| `soap\05_文献\2020_Fujii_NatCommun_GBThermalConductivity_LocalAtomicEnvironments_正式版.pdf` | 4.96 MB | 文献 | `355cba548cd681e2` | PDF 文档。 |
| `grain_boundary\history\legacy\GB_project_20260819\archive\old_versions\legacy_generator\lammps_workflow\traj_quench.dump` | 4.92 MB | 中间 | `64967a2841e6cdb1` | LAMMPS 轨迹：该阶段的过程演化，多帧记录。 |
| `grain_boundary\history\experiments\phase1_2026-09\Li_Sigma5_310_2x2_阶段一零压结果_2026-09-03\5_保持10K零压\p1_hold.dump` | 4.68 MB | 中间 | `e5a86bae9de6124f` | 10 K 保持段的过程轨迹。用于确认无序结构在整个保持期内持续存在，而非中途塌陷。 |
| `grain_boundary\00_Literature\lithium_metal\Diﬀusion and Deformation Mechanism Maps in Li.pdf` | 4.63 MB | 文献 | `05422c659103a232` | PDF 文档。 |
| `grain_boundary\00_Literature\soap\2018_Priedeman_BYU硕士论文_QuantifyingGBAtomicStructures_SOAP.pdf` | 4.35 MB | 文献 | `a73bc6e1a2f67fe1` | PDF 文档。 |
| `soap\05_文献\2018_Priedeman_BYU硕士论文_QuantifyingGBAtomicStructures_SOAP.pdf` | 4.35 MB | 文献 | `a73bc6e1a2f67fe1` | PDF 文档。 |
| `grain_boundary\00_Literature\soap\Ni SOAP报告.pdf` | 4.19 MB | 文献 | `5efaf15dd64a5b88` | PDF 文档。 |
| `grain_boundary\00_Literature\soap\2022_SongDeng_PRM_AtomicEnergyInGrainBoundaries_MachineLearning.pdf` | 3.31 MB | 文献 | `790c6490db7446dd` | PDF 文档。 |
| `soap\05_文献\2022_SongDeng_PRM_AtomicEnergyInGrainBoundaries_MachineLearning.pdf` | 3.31 MB | 文献 | `790c6490db7446dd` | PDF 文档。 |
| `grain_boundary\00_Literature\soap\2017_Rosenbrock_npjComputMater_DiscoveringBuildingBlocks_ML_GBs.pdf` | 3.00 MB | 文献 | `46aaf24b8165707e` | PDF 文档。 |
| `soap\05_文献\2017_Rosenbrock_npjComputMater_DiscoveringBuildingBlocks_ML_GBs.pdf` | 3.00 MB | 文献 | `46aaf24b8165707e` | PDF 文档。 |
| `grain_boundary\00_Literature\grain_boundary\Grain Boundary Softening A Potential Mechanism for Lithium Metal.pdf` | 2.95 MB | 文献 | `65a2a36a026160fe` | PDF 文档。 |
| `grain_boundary\00_Literature\soap\Feature engineering descriptors, transforms, and machine learning for grain boundaries and variable-sized atom clusters.pdf` | 2.77 MB | 文献 | `a928ea2c90b10937` | PDF 文档。 |
| `soap\05_文献\Feature engineering descriptors, transforms, and machine learning for grain boundaries and variable-sized atom clusters.pdf` | 2.77 MB | 文献 | `a928ea2c90b10937` | PDF 文档。 |
| `grain_boundary\00_Literature\soap\1-s2.0-S1359645424010450-main.pdf` | 2.67 MB | 文献 | `0c69d7126e9df3a7` | PDF 文档。 |
| `soap\05_文献\1-s2.0-S1359645424010450-main.pdf` | 2.67 MB | 文献 | `0c69d7126e9df3a7` | PDF 文档。 |
| `grain_boundary\00_Literature\soap\2017_Rosenbrock_npjComputMater_SI_SOAP参数选择_补充材料.pdf` | 2.25 MB | 文献 | `e5bc0c22e26e7531` | PDF 文档。 |
| `soap\05_文献\2017_Rosenbrock_npjComputMater_SI_SOAP参数选择_补充材料.pdf` | 2.25 MB | 文献 | `e5bc0c22e26e7531` | PDF 文档。 |
| `grain_boundary\00_Literature\grain_boundary\Atomistic modeling of the mechanical.pdf` | 2.21 MB | 文献 | `1f37591627c08f67` | PDF 文档。 |
| `grain_boundary\00_Literature\lithium_metal\Elastic and plastic mechanical properties of lithium measured .pdf` | 2.08 MB | 文献 | `5ce23d3ad6e386a7` | PDF 文档。 |
| `grain_boundary\history\legacy\GB_project_20260819\archive\old_versions\legacy_generator\lammps_workflow\traj_anneal.dump` | 1.89 MB | 中间 | `d4c9ee4393511f28` | LAMMPS 轨迹：该阶段的过程演化，多帧记录。 |
| `grain_boundary\00_Literature\soap\1-s2.0-S1359645415300860-main.pdf` | 1.75 MB | 文献 | `91adda1828777ead` | PDF 文档。 |
| `soap\05_文献\1-s2.0-S1359645415300860-main.pdf` | 1.75 MB | 文献 | `91adda1828777ead` | PDF 文档。 |
| `grain_boundary\00_Literature\soap\2020_Wagih_NatCommun_GBSegregationEnergySpectra_Polycrystals.pdf` | 1.28 MB | 文献 | `caf2af0725887b63` | PDF 文档。 |
| `soap\05_文献\2020_Wagih_NatCommun_GBSegregationEnergySpectra_Polycrystals.pdf` | 1.28 MB | 文献 | `caf2af0725887b63` | PDF 文档。 |
| `grain_boundary\00_Literature\lithium_metal\The growth mechanism of lithium dendrites and its.pdf` | 756.1 KB | 文献 | `b082f9a80de4dfc6` | PDF 文档。 |
| `grain_boundary\history\presentations\presentations\Li_metal_metastable_atomic_structure_mechanical_analysis.pdf` | 678.0 KB | 文献 | `23bd55d59dc7b96a` | PDF 文档。 |
| `grain_boundary\00_Literature\soap\笔记_SOAP原理_结构有序无序度表征.pdf` | 421.7 KB | 文献 | `093369c5d22e2e25` | PDF 文档。 |
| `soap\05_文献\笔记_SOAP原理_结构有序无序度表征.pdf` | 421.7 KB | 文献 | `093369c5d22e2e25` | PDF 文档。 |
| `grain_boundary\00_Literature\soap\笔记_SOAP原理_信息来源与参考文献.pdf` | 358.5 KB | 文献 | `aff3850764f097b2` | PDF 文档。 |
| `soap\05_文献\笔记_SOAP原理_信息来源与参考文献.pdf` | 358.5 KB | 文献 | `aff3850764f097b2` | PDF 文档。 |
| `soap\Stage02_SOAP表示与审计\03_输出\stage2_5_raw_soap_vectors.npy` | 261.5 KB | 中间 | `e91756c89f16dfeb` | Stage02.5 表示审计产物。 |
| `grain_boundary\00_Literature\soap\NI_SOAP_PROOF_FINAL_BI.html` | 67.0 KB | 文献 | `e9b873ad47ffbd61` | HTML 报告（含公式渲染）。 |
| `soap\Stage01_数据准备与SOAP初测\03_输出\soap_vectors_metastable.npz` | 39.5 KB | 中间 | `195920ac2d0414b4` | NumPy 压缩数组（中间计算量）。 |
| `soap\Stage01_数据准备与SOAP初测\03_输出\soap_vectors_stable.npz` | 39.5 KB | 中间 | `beaebc8598617932` | NumPy 压缩数组（中间计算量）。 |
| `tools\__pycache__\build_artifact_catalog.cpython-313.pyc` | 35.5 KB | 文档 | `c4386a87c7da6294` | （未分类） |
| `tools\__pycache__\build_artifact_catalog.cpython-312.pyc` | 35.1 KB | 文档 | `5f1438b0afe15528` | （未分类） |
| `grain_boundary\00_Literature\soap\NI_SOAP_PROOF_FINAL_BI_files\katex.min.css` | 23.3 KB | 文献 | `0289a02cf451a44d` | 网页样式表（随 HTML 报告）。 |
| `soap\Stage02_SOAP表示与审计\03_输出\stage2_75_bulk_reference_vectors.npy` | 16.0 KB | 中间 | `81a636b05b9c1251` | Stage02.75 全人口 SOAP 与距离审计产物。 |
| `grain_boundary\history\legacy\GB_project_20260819\archive\legacy\li_gb_project\lammps_test\log.lammps` | 8.8 KB | 追溯 | `ee4069064f9d2f57` | LAMMPS 运行日志，含各阶段数值历史与断言结果。 |

## 还原步骤

1. 从备份介质取回上述文件，保持相对路径不变
2. 核对校验和：`certutil -hashfile <文件> SHA256`，比对上表
3. 重新生成台账：双击 `refresh_ledger.cmd`
