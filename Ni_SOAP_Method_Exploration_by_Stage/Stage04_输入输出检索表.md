# Stage 04：连续无序度与能量关联　输入・输出检索表

> 仓库 [`Star-1222/li-grain-boundary-soap`](https://github.com/Star-1222/li-grain-boundary-soap)　·　目录 `Ni_SOAP_Method_Exploration_by_Stage/Stage04_无序度与能量关联/`
> 本文件把 Stage 04 单独拎出来，逐个文件给出作用。链接指向 `main` 分支。

## 这一阶段在做什么

放弃离散类别，改用连续坐标 D(k)：一个原子到它最近的 k 个「清洁体相参考原子」的平均距离。先检查 D 对 k、对参考体制是否稳健，再拿 `c_PEATOM` 当靶子，看 D 到底携带多少信息。

## 数据流

```
stage2_75_soap_vectors.npy ─┐
stage2_75_atom_metadata.csv ─┴─→ stage3a.load_and_validate() ─→ X, meta
                                        │
                                        ├─→ stage4a_continuous_disorder.py
                                        │      └─→ D(k=1,5,10,20) + 964 个清洁体相参考
                                        │      └─→ k 稳健性 / 参考敏感性 / bootstrap / PCA
                                        │
                                        └─→ stage4b_disorder_energy.py
                                               └─→ 三种基线的过剩能 + 分箱趋势 + 共同支撑
```

## 一、输入

本阶段没有外部数据，吃的是上游产物。

| 文件 | 从哪来 | 作用 |
|---|---|---|
| `stage2_75_soap_vectors.npy` | Stage02（未入库） | 通过 `stage3a.load_and_validate()` 读入。 |
| [`stage2_75_atom_metadata.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/stage2_75_atom_metadata.csv) | Stage02 | 同上，由共用 loader 一起读。 |
| [`stage3a_lae_threshold.py`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/02_%E8%84%9A%E6%9C%AC/stage3a_lae_threshold.py) | Stage03 | **代码依赖**：`import stage3a_lae_threshold as s3a`。4A 与 4B 都靠它取数据。 |

> 注意：这些脚本的 `OUTPUT_DIR` 都是「脚本所在目录 / outputs」，而且彼此用裸模块名 `import`，
> 运行时所有脚本与数据必须放在同一个目录里。归档时按阶段拆成了 `03_输出/`，
> 所以重跑之前要先把上游文件复制到本阶段的 `02_脚本/` 下。

## 二、处理代码

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`stage4a_continuous_disorder.py`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/02_%E8%84%9A%E6%9C%AC/stage4a_continuous_disorder.py) | 27.6 KB | 定义连续坐标 D(k)：到 k 个最近清洁体相参考原子的平均 `d_raw`。做 k 稳健性、参考体制敏感性、bootstrap 与 PCA。 |
| [`stage4b_disorder_energy.py`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/02_%E8%84%9A%E6%9C%AC/stage4b_disorder_energy.py) | 33.2 KB | 把 D 与 `c_PEATOM` 对上：三种基线、分箱趋势、共同支撑比较、固定 D 区间的能量离散。 |

## 三、输出：数据表

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`stage4a_atomwise_disorder.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4a_atomwise_disorder.csv) | 997.6 KB | 4659 行 × 17 列。`disorder_k1/k5/k10/k20`、三个主成分、以及 `in_clean_bulk_reference` 标记哪 964 个原子被用作参考。 |
| [`stage4a_bootstrap_summary.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4a_bootstrap_summary.csv) | 942 B | 关键量的 bootstrap 区间。**界面 D 差值的 95% 区间含零就在这里**。 |
| [`stage4a_group_statistics.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4a_group_statistics.csv) | 3.4 KB | 四组（state × region）× 四个 k 的均值、中位数、分位数。 |
| [`stage4a_k_robustness.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4a_k_robustness.csv) | 1.5 KB | k 两两之间的 Pearson / Spearman 与绝对差。k=5/10/20 相关高于 0.99。 |
| [`stage4a_pca_variance.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4a_pca_variance.csv) | 459 B | 前几个主成分的解释率与累计解释率。 |
| [`stage4a_reference_sensitivity.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4a_reference_sensitivity.csv) | 1.2 KB | pooled 参考与分状态参考的对比。排序相关约 0.999。 |
| [`stage4b_atomwise_energy_disorder.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_atomwise_energy_disorder.csv) | 1.26 MB | 4659 行 × 21 列。在 4A 基础上加了三个基线的过剩能，是 Stage05 与 Stage06 的主要输入。 |
| [`stage4b_binned_interface_trend.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_binned_interface_trend.csv) | 2.6 KB | 按 D 分箱后的能量中位数与四分位。**看是否单调**（结论：高 D 尾部下降）。 |
| [`stage4b_binned_metastable_trend.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_binned_metastable_trend.csv) | 1.7 KB | 亚稳状态单独分箱的趋势。 |
| [`stage4b_binned_stable_trend.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_binned_stable_trend.csv) | 1.7 KB | 稳定状态单独分箱的趋势。 |
| [`stage4b_bootstrap_summary.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_bootstrap_summary.csv) | 3.2 KB | 带 `interval_type` 列的区间估计，比 4A 更完整。 |
| [`stage4b_common_support_comparison.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_common_support_comparison.csv) | 1.1 KB | 只在两状态 D 重叠区间内比较，避免外推。 |
| [`stage4b_correlations.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_correlations.csv) | 1.2 KB | 18 行 = 基线 × 子集的 D–能量相关。界面 Spearman 约 0.28。 |
| [`stage4b_energy_baselines.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_energy_baselines.csv) | 90 B | 三个基线各自的原子数与 PEATOM 中位数（pooled 基线 −4.450650 eV/atom）。 |
| [`stage4b_fixed_disorder_energy_spread.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_fixed_disorder_energy_spread.csv) | 761 B | 固定 D 区间内能量的离散程度，说明「同样无序度可以对应多大能量差」。 |
| [`stage4b_group_energy_statistics.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_group_energy_statistics.csv) | 877 B | 四组过剩能的均值、分位数、IQR。 |
| [`stage4b_k_robustness.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_k_robustness.csv) | 601 B | 换 k 后 D 与能量的相关是否稳定。 |

## 四、输出：图

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`stage4a_bulk_disorder_distribution.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4a_bulk_disorder_distribution.png) | 19.7 KB | 体相原子的 D 分布。中位数约 0.03，说明参考定义成立。 |
| [`stage4a_four_group_disorder.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4a_four_group_disorder.png) | 23.8 KB | 四组 D 的箱线对照。 |
| [`stage4a_interface_disorder_distribution.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4a_interface_disorder_distribution.png) | 22.8 KB | 界面原子的 D 分布（两个状态叠在一起）。 |
| [`stage4a_k_sensitivity.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4a_k_sensitivity.png) | 17.2 KB | 不同 k 下 D 的分布与排序，用来判断 k=10 是否只是数值平滑。 |
| [`stage4a_pca_disorder.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4a_pca_disorder.png) | 61.5 KB | PCA 与 D 的对照。 |
| [`stage4a_pca_groups.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4a_pca_groups.png) | 31.8 KB | PCA 前两维按四组着色。 |
| [`stage4a_pca_peatom.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4a_pca_peatom.png) | 60.9 KB | PCA 与 PEATOM 的对照。 |
| [`stage4b_all_atoms_scatter.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_all_atoms_scatter.png) | 30.6 KB | 全体原子的 D 与过剩能散点。 |
| [`stage4b_conditional_state_difference.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_conditional_state_difference.png) | 21.3 KB | 共同 D 区间内两状态的能量差。 |
| [`stage4b_interface_binned_trend.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_interface_binned_trend.png) | 22.5 KB | 界面原子分箱后的能量中位数曲线。 |
| [`stage4b_interface_scatter.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_interface_scatter.png) | 28.8 KB | 只取界面原子的散点。相关 0.258 / 0.277 来自这里。 |
| [`stage4b_metastable_interface_scatter.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_metastable_interface_scatter.png) | 27.3 KB | 亚稳结构界面原子的散点。 |
| [`stage4b_pca_disorder.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_pca_disorder.png) | 61.5 KB | PCA 与 D 的对照。 |
| [`stage4b_pca_energy.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_pca_energy.png) | 62.3 KB | PCA 与过剩能的对照。 |
| [`stage4b_pca_interface_energy.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_pca_interface_energy.png) | 60.4 KB | 只取界面原子的 PCA 与能量对照。 |
| [`stage4b_stable_interface_scatter.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_stable_interface_scatter.png) | 11.8 KB | 稳定结构界面原子的散点。 |
| [`stage4b_state_binned_trends.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_state_binned_trends.png) | 23.8 KB | 两个状态的分箱趋势并排，看形状是否一致。 |

## 五、输出：结构文件

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`stage4a_metastable_disorder.xyz`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4a_metastable_disorder.xyz) | 136.2 KB | 亚稳结构版本。 |
| [`stage4a_stable_disorder.xyz`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4a_stable_disorder.xyz) | 111.7 KB | 稳定结构按 D 着色的结构文件，OVITO 打开看空间分布。 |

## 六、未入库的产物

本阶段全部产物都已入库。

> 未入库项登记在 `DATA_MANIFEST.md`，附 SHA256 前 16 位，用于换机器后核对。

## 七、报告

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`Stage04_确定性结论与待讨论.md`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/01_%E6%8A%A5%E5%91%8A/Stage04_%E7%A1%AE%E5%AE%9A%E6%80%A7%E7%BB%93%E8%AE%BA%E4%B8%8E%E5%BE%85%E8%AE%A8%E8%AE%BA.md) | 3.7 KB | 体相参考、连续 D、D 的数值组成、辅助标量关联四节。关键一句是「D 不是热力学无序度」。 |

## 八、本阶段的关键数字

| 项 | 值 | 出处 |
|---|---|---|
| 清洁体相参考 | 964 个原子（stable 476 + metastable 488） | 条件：`region=bulk` + PTM=1 + CNA=1 + PEATOM 在 Q1−3IQR 内 |
| 主指标 | D(k=10) | 到 10 个最近参考原子的平均 `d_raw` |
| k 稳健性 | k = 5/10/20 之间相关 > 0.99 | 所以 k=10 不是关键选择 |
| 参考体制 | pooled vs 分状态，排序相关约 0.999 | 两种参考给出几乎相同的排序 |
| 界面 D 中位数 | stable 6.744 / metastable 5.972，差 0.772 | 看起来 stable 更无序 |
| 该差值的 95% 区间 | −0.129 ~ 1.601，**含零** | 所以不能声称状态间的无序度排序 |
| 体相 D 中位数 | 0.03 | 参考原子自身的 D 接近零，说明定义成功 |
| PEATOM 基线 | pooled −4.450650 eV/atom | 后面所有过剩能的零点 |
| D 与能量相关 | 全体 0.449 / 界面 0.258（Pearson） | 界面 Spearman 0.277 |

## 九、本阶段交给下游的东西

- **`stage4b_atomwise_energy_disorder.csv`** 是 Stage05、Stage06 的直接输入，里面带着三种口径的过剩能。
- **一条禁令**：D 只能叫描述性连续坐标，不能叫热力学无序度。它用的仍是 `d_raw`，仍混着 42–46% 的范数效应。
- **一个警示**：界面 D 的状态差区间含零，所以后续任何「亚稳更无序」的说法都缺依据。

## 十、归档缺口

`stage4a` 与 `stage4b` 各写一份 `stage4a_report.md` / `stage4b_report.md`，归档里没有。两个脚本都不直接 `np.load`，而是通过 import 的 `s3a.load_and_validate()` 取数据，所以运行时 `stage3a_lae_threshold.py` 必须和它们在同一目录。
