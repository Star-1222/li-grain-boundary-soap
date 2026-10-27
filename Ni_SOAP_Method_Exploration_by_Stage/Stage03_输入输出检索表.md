# Stage 03：LAE 阈值与 LAE/LER 构建　输入・输出检索表

> 仓库 [`Star-1222/li-grain-boundary-soap`](https://github.com/Star-1222/li-grain-boundary-soap)　·　目录 `Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE阈值与LAE_LER构建/`
> 本文件把 Stage 03 单独拎出来，逐个文件给出作用。链接指向 `main` 分支。

## 这一阶段在做什么

在原始 SOAP 距离上扫 115 个 epsilon，用贪心算法建局部原子环境（LAE）类别，检查它对原子顺序有多敏感、邻域图是否连通，最后建 LER 人口指纹。结论是否定：离散分类不可用。

## 数据流

```
stage2_75_soap_vectors.npy ─┐
                                ├─→ stage3a_lae_threshold.py ─→ epsilon 扫描 115 × 30 次
stage2_75_atom_metadata.csv ────┘            ├─→ partition stability / ARI
                                             └─→ 候选区间 7.224–8.309
                                                        │（人工转录成常数）
                                                        ▼
                              stage3b_lae_ler.py ─→ 5 个 epsilon 的 LAE + LER
```

## 一、输入

本阶段没有外部数据，吃的是上游产物。

| 文件 | 从哪来 | 作用 |
|---|---|---|
| `stage2_75_soap_vectors.npy` | Stage02（未入库） | 主矩阵 4659 × 1014，`stage3a` 第 100 行 `np.load` 读它。 |
| [`stage2_75_atom_metadata.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/stage2_75_atom_metadata.csv) | Stage02 | 行对齐的元数据，提供坐标与标签。 |
| [`stage2_75_population_summary.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/stage2_75_population_summary.csv) | Stage02 | 人口分位数，用来判断类别划分是否合理。 |

> 注意：这些脚本的 `OUTPUT_DIR` 都是「脚本所在目录 / outputs」，而且彼此用裸模块名 `import`，
> 运行时所有脚本与数据必须放在同一个目录里。归档时按阶段拆成了 `03_输出/`，
> 所以重跑之前要先把上游文件复制到本阶段的 `02_脚本/` 下。

## 二、处理代码

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`stage3a_lae_threshold.py`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/02_%E8%84%9A%E6%9C%AC/stage3a_lae_threshold.py) | 56.6 KB | epsilon 扫描主体：115 个阈值 × 30 次随机原子顺序，统计类数、单例比例、最大类占比、连通分量。同时提供 `load_and_validate()`，是 Stage04–06 共用的数据入口。 |
| [`stage3b_lae_ler.py`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/02_%E8%84%9A%E6%9C%AC/stage3b_lae_ler.py) | 39.4 KB | 用 Stage3A 定下的 5 个 epsilon 建 LAE 与 LER，算运行间 ARI、类紧致度、两状态的 LER 距离。 |

## 三、输出：数据表

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`stage3a_epsilon_grid.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3a_epsilon_grid.csv) | 5.0 KB | 115 个 epsilon 的网格值与来源（粗扫还是精扫）。 |
| [`stage3a_epsilon_sweep.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3a_epsilon_sweep.csv) | 89.4 KB | 115 行 × 55 列。每个 epsilon 在 30 次随机顺序下的统计，`_mean` / `_std` 成对出现。`n_classes_cv` 直接度量顺序依赖。 |
| [`stage3a_nearest_neighbor_distances.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3a_nearest_neighbor_distances.csv) | 331.0 KB | 4659 个原子各自到最近邻原子的距离，是判断 epsilon 该取多大尺度的依据。 |
| [`stage3a_nearest_neighbor_summary.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3a_nearest_neighbor_summary.csv) | 2.0 KB | 上面那列的分布（按结构 / 区域分组）。 |
| [`stage3a_partition_stability.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3a_partition_stability.csv) | 2.4 KB | 每个 epsilon 的 ARI 与连通分量统计。**`connected_component_count` 恒为 1 就是 chaining 的直接证据**。 |
| [`stage3a_refined_epsilon_sweep.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3a_refined_epsilon_sweep.csv) | 28.1 KB | 候选区间 7.224–8.309 内的精细扫描（35 个点）。 |
| [`stage3a_threshold_candidates.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3a_threshold_candidates.csv) | 281 B | 只挑出的那 1 行候选区间，以及该区间满足的各项门槛。 |
| [`stage3b_atomwise_mid_epsilon.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3b_atomwise_mid_epsilon.csv) | 611.5 KB | 4659 个原子在 mid epsilon 下被分到哪个 `class_id`，以及到该类的代表原子的距离。 |
| [`stage3b_canonical_classes.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3b_canonical_classes.csv) | 2.6 KB | 每个类别在四组（stable / metastable × interface / bulk）里的人口与占比。 |
| [`stage3b_class_compactness.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3b_class_compactness.csv) | 3.0 KB | 类内直径、p95 内距、代表原子半径、`non_mutual_equivalence_fraction`（这个数越高越说明类别边界是假的）。 |
| [`stage3b_class_statistics.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3b_class_statistics.csv) | 33.8 KB | 长格式表：`epsilon, class_id, stat, group, value`，一个类占多行。查某个类的 PEATOM 中位数就翻这里。 |
| [`stage3b_epsilon_set.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3b_epsilon_set.csv) | 544 B | 最终采用的 5 个 epsilon 及角色（below / low / mid / high / above）与理由。 |
| [`stage3b_ler_by_run.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3b_ler_by_run.csv) | 42.1 KB | 每次运行下每个类的人口（occupancy）与两个状态各自的 LER。 |
| [`stage3b_ler_comparison.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3b_ler_comparison.csv) | 13.1 KB | stable 与 metastable 的 LER 距离，三种度量：L1、JSD、余弦。 |
| [`stage3b_ler_summary.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3b_ler_summary.csv) | 3.5 KB | LER 的跨运行汇总：中位数、均值、标准差。 |
| [`stage3b_partition_ari.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3b_partition_ari.csv) | 67.3 KB | 2175 行 = 每个 epsilon 下 30 次运行的 435 个两两配对 ARI。 |
| [`stage3b_representative_atoms.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3b_representative_atoms.csv) | 6.5 KB | 每个类的代表原子（medoid）及其坐标、标签、PEATOM。 |
| [`stage3b_run_summary.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3b_run_summary.csv) | 7.4 KB | 150 行 = 5 个 epsilon × 30 次运行。每次运行的类数、与其他运行的平均 ARI、是否为 medoid。 |

## 四、输出：图

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`stage3a_bulk_behavior_vs_epsilon.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3a_bulk_behavior_vs_epsilon.png) | 32.0 KB | 体相原子落入主导类的比例随 epsilon 的变化。 |
| [`stage3a_class_count_vs_epsilon.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3a_class_count_vs_epsilon.png) | 21.0 KB | 类数随 epsilon 的变化（中位数 + 波动带）。看有没有平台区。 |
| [`stage3a_nn_distributions.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3a_nn_distributions.png) | 25.1 KB | 最近邻距离的分布图，定 epsilon 的量级参考。 |
| [`stage3a_order_dependence_vs_epsilon.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3a_order_dependence_vs_epsilon.png) | 22.4 KB | 顺序依赖随 epsilon 的变化。曲线越抖说明越不可靠。 |
| [`stage3a_partition_stability_vs_epsilon.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3a_partition_stability_vs_epsilon.png) | 24.8 KB | 运行间 ARI 随 epsilon 的变化。 |
| [`stage3a_refined_threshold_plot.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3a_refined_threshold_plot.png) | 25.9 KB | 候选区间的精扫结果，用来看阈值是否真的落在平台上。 |
| [`stage3a_singleton_fraction_vs_epsilon.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3a_singleton_fraction_vs_epsilon.png) | 19.2 KB | 单例类比例随 epsilon 的变化。 |
| [`stage3b_class_compactness.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3b_class_compactness.png) | 18.2 KB | 类紧致度（直径、p95 内距）随 epsilon 的变化。 |
| [`stage3b_class_population_vs_epsilon.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3b_class_population_vs_epsilon.png) | 18.3 KB | 各类人口随 epsilon 的变化。 |
| [`stage3b_lae_vs_z_metastable.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3b_lae_vs_z_metastable.png) | 13.1 KB | 亚稳结构版本。 |
| [`stage3b_lae_vs_z_stable.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3b_lae_vs_z_stable.png) | 13.2 KB | 稳定结构里 LAE 类别随 z 的分布，看类别是否按界面分层。 |
| [`stage3b_ler_stable_vs_metastable.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3b_ler_stable_vs_metastable.png) | 21.6 KB | 两个状态的 LER 并排对照。L1 约 0.37–0.45 就来自这张图。 |
| [`stage3b_ler_vs_epsilon.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3b_ler_vs_epsilon.png) | 27.1 KB | LER 排名曲线随 epsilon 的变化。 |

## 五、输出：结构文件

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`stage3b_mid_epsilon_lae.xyz`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/03_%E8%BE%93%E5%87%BA/stage3b_mid_epsilon_lae.xyz) | 543.2 KB | 4659 个原子按类别着色的结构文件，OVITO 直接打开。 |

## 六、未入库的产物

本阶段全部产物都已入库。

> 未入库项登记在 `DATA_MANIFEST.md`，附 SHA256 前 16 位，用于换机器后核对。

## 七、报告

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`Stage03_确定性结论与待讨论.md`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/01_%E6%8A%A5%E5%91%8A/Stage03_%E7%A1%AE%E5%AE%9A%E6%80%A7%E7%BB%93%E8%AE%BA%E4%B8%8E%E5%BE%85%E8%AE%A8%E8%AE%BA.md) | 3.1 KB | epsilon 扫描、候选区间、排序依赖、chaining、LER 结果五节，每节都给出「能确定」和「不能确定」。结论是离散 LAE 不采用。 |

## 八、本阶段的关键数字

| 项 | 值 | 出处 |
|---|---|---|
| epsilon 网格 | 115 个点，0.0038 → 13.08 | `stage3a_epsilon_grid.csv` |
| 随机顺序重复 | 30 次 / 每个 epsilon | 用来度量顺序依赖 |
| 运行间 ARI | 0.56–0.71 | 换一次随机顺序，分类就变了 |
| 类数范围 | 4–8 | 同一个 epsilon 下都不稳定 |
| 连通分量 | 1（而类数约 5） | **chaining**：阈值图完全连通，类别边界是假的 |
| 候选区间 | 7.224 → 8.309 | 类数中位数稳定为 5，单例比例为 0 |
| LER 距离 | L1 0.37–0.45 / JSD 0.046–0.073 / 余弦 0.876–0.902 | 能区分两个状态，但只是人口指纹 |

## 九、本阶段交给下游的东西

- **主分析取消离散 LAE**：后续阶段不再继承「5 类物理基元」的解释。
- **LER 只作描述性指纹**：不用于因果推断。
- **一个副作用**：`stage3a_lae_threshold.py` 里的 `load_and_validate()` 成了 Stage04–06 共用的数据入口，是整条链的枢纽代码。

## 十、归档缺口

`stage3b_lae_ler.py` 不读 Stage3A 的输出 CSV，五个 epsilon 值（7.224 / 8.309 / 6.790452 / 8.525591 / 中点）是**硬编码常数**，靠人工从 3A 结果转录。两个脚本还会各写 `stage3a_report.md` / `stage3b_report.md`，归档里没有。
