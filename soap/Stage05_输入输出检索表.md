# Stage 05：物理描述符与 SOAP 方向　输入・输出检索表

> 仓库 [`Star-1222/li-grain-boundary-soap`](https://github.com/Star-1222/li-grain-boundary-soap)　·　目录 `Ni_SOAP_Method_Exploration_by_Stage/Stage05_物理空间与方向/`
> 本文件把 Stage 05 单独拎出来，逐个文件给出作用。链接指向 `main` 分支。

## 这一阶段在做什么

把「结构无序」拆成可计算的物理量（配位数、键应变、Voronoi 体积、CSP），用它们对照 D 解释不掉的那部分能量；再给 SOAP 向量定义「位移方向」，在匹配 D 的条件下看方向能不能解释能量差。

## 数据流

```
stage4b_atomwise_energy_disorder.csv ─┐
ni_stable.cfg / ni_metastable.cfg ────────┼─→ stage5a_physical_spatial_residual.py
                                          │      └─→ 10 个描述符 + 残差能 + 空间自相关
                                          │
stage2_75_soap_vectors.npy ───────────────┴─→ stage5b_soap_direction.py
                                                 └─→ delta_s / u + 匹配对 + 置换检验
```

## 一、输入

本阶段没有外部数据，吃的是上游产物。

| 文件 | 从哪来 | 作用 |
|---|---|---|
| [`stage4b_atomwise_energy_disorder.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_atomwise_energy_disorder.csv) | Stage04 | 过剩能与 D，是残差能定义的基础。 |
| [`stage2_75_atom_metadata.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/stage2_75_atom_metadata.csv) | Stage02 | 5B 直接 `pd.read_csv` 读它。 |
| `stage2_75_soap_vectors.npy` | Stage02（未入库） | 5B 第 70 行 `np.load`，用于算方向。 |
| [`ni_stable.cfg`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage01_%E6%95%B0%E6%8D%AE%E5%87%86%E5%A4%87%E4%B8%8ESOAP%E5%88%9D%E6%B5%8B/02_%E8%84%9A%E6%9C%AC/ni_stable.cfg) | Stage01 | 5A 要重新读结构算 CN、Voronoi 等几何量。 |
| [`ni_metastable.cfg`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage01_%E6%95%B0%E6%8D%AE%E5%87%86%E5%A4%87%E4%B8%8ESOAP%E5%88%9D%E6%B5%8B/02_%E8%84%9A%E6%9C%AC/ni_metastable.cfg) | Stage01 | 同上。 |
| [`stage3a_lae_threshold.py`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/02_%E8%84%9A%E6%9C%AC/stage3a_lae_threshold.py) | Stage03 | **代码依赖**：5A、5B 都 import 它。 |

> 注意：这些脚本的 `OUTPUT_DIR` 都是「脚本所在目录 / outputs」，而且彼此用裸模块名 `import`，
> 运行时所有脚本与数据必须放在同一个目录里。归档时按阶段拆成了 `03_输出/`，
> 所以重跑之前要先把上游文件复制到本阶段的 `02_脚本/` 下。

## 二、处理代码

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`stage5a_physical_spatial_residual.py`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/02_%E8%84%9A%E6%9C%AC/stage5a_physical_spatial_residual.py) | 49.9 KB | 算 10 个物理描述符（CN、近邻距离统计、键应变、CSP、Voronoi 体积、超额体积），定义残差能 = 过剩能 − 同 D 分箱中位数，并做空间自相关与 patchiness。 |
| [`stage5b_soap_direction.py`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/02_%E8%84%9A%E6%9C%AC/stage5b_soap_direction.py) | 46.7 KB | 构造 SOAP 位移方向 `delta_s` 与单位方向 `u`，在匹配 D 的条件下比较方向距离与能量差，做置换检验与 PCA。 |

## 三、输出：数据表

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`stage5a_atomwise_physical_descriptors.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_atomwise_physical_descriptors.csv) | 1.44 MB | 4659 行 × 23 列，每个原子的全部物理描述符与残差能。**Stage05 的中心表**，Stage06 也读它。 |
| [`stage5a_descriptor_correlations.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_descriptor_correlations.csv) | 2.3 KB | 30 行 = 子集 × 描述符，与残差能的 Pearson / Spearman。 |
| [`stage5a_descriptor_group_statistics.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_descriptor_group_statistics.csv) | 5.2 KB | 按组统计各描述符。 |
| [`stage5a_disorder_bin_conditional_stats.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_disorder_bin_conditional_stats.csv) | 1.5 KB | 在每个 D 箱内单独算描述符与残差的相关，用来确认不是 D 在作怪。 |
| [`stage5a_matched_pairs.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_matched_pairs.csv) | 32.7 KB | D 匹配的原子对（容差 0.10），72 行。带两侧的全部描述符。 |
| [`stage5a_residual_definition.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_residual_definition.csv) | 251 B | 残差定义的稳健性：18 箱为主，12 / 15 / 20 箱为敏感性，给出与主定义的 Spearman。 |
| [`stage5a_spatial_autocorrelation.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_spatial_autocorrelation.csv) | 16.3 KB | 面内自相关：每个状态、每条界面、每个场、每个半径。 |
| [`stage5a_spatial_patchiness.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_spatial_patchiness.csv) | 982 B | top-10% 热点原子的近邻距离与 3 / 5 / 7 Å 内共现比例。比例低说明热点不成块。 |
| [`stage5b_atomwise_soap_direction.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5b_atomwise_soap_direction.csv) | 1.41 MB | 4659 行 × 25 列。`displacement_norm`（幅度）、`direction_norm`、`disp_PC1-3`、`dir_PC1-3`，以及全部物理描述符。 |
| [`stage5b_direction_autocorrelation.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5b_direction_autocorrelation.csv) | 2.4 KB | 方向的空间自相关：每个状态、每条界面、每个半径的平均余弦。用来判断方向是否成块，还是随机散布。 |
| [`stage5b_directional_energy_statistics.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5b_directional_energy_statistics.csv) | 3.4 KB | 按容差与类别统计方向与能量差的关系，含方向五分位下的中位能量差。 |
| [`stage5b_directional_statistics.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5b_directional_statistics.csv) | 652 B | 三个类别（stable-stable / meta-meta / stable-meta）的余弦统计。 |
| [`stage5b_pca_summary.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5b_pca_summary.csv) | 1.1 KB | 位移 PCA 与方向 PCA 的解释率。 |
| [`stage5b_permutation_null.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5b_permutation_null.csv) | 227 B | 只 1 行，但决定 0.83 算不算真信号：实测 Spearman vs 置换零模型的均值与分位。 |
| [`stage5b_representative_cases.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5b_representative_cases.csv) | 3.9 KB | 8 个代表案例（方向相似但能量差大的对）。 |

## 四、输出：图

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`stage5a_matched_pair_energy_vs_disorder.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_matched_pair_energy_vs_disorder.png) | 27.2 KB | 匹配对的能量差 vs D 差。 |
| [`stage5a_metastable_bond_strain_map.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_metastable_bond_strain_map.png) | 62.3 KB | 亚稳结构键应变空间图。 |
| [`stage5a_metastable_disorder_map.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_metastable_disorder_map.png) | 62.5 KB | 亚稳结构 D 的空间图。 |
| [`stage5a_metastable_physical_descriptors.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_metastable_physical_descriptors.png) | 56.3 KB | 亚稳结构版本。 |
| [`stage5a_metastable_residual_energy_map.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_metastable_residual_energy_map.png) | 60.0 KB | 亚稳结构残差能空间图。 |
| [`stage5a_residual_vs_bond_strain.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_residual_vs_bond_strain.png) | 27.5 KB | 残差能 vs 键应变（−0.503，最强的单变量）。 |
| [`stage5a_residual_vs_cn.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_residual_vs_cn.png) | 14.9 KB | 残差能 vs 配位数（Spearman −0.392）。 |
| [`stage5a_residual_vs_csp.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_residual_vs_csp.png) | 23.9 KB | 残差能 vs 中心对称参数（+0.202，很弱）。 |
| [`stage5a_residual_vs_disorder.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_residual_vs_disorder.png) | 26.7 KB | 残差能 vs D，用来确认残差确实与 D 解耦。 |
| [`stage5a_residual_vs_excess_volume.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_residual_vs_excess_volume.png) | 25.8 KB | 残差能 vs 超额 Voronoi 体积（−0.360）。 |
| [`stage5a_spatial_autocorrelation.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_spatial_autocorrelation.png) | 59.2 KB | 自相关曲线：相关随半径衰减得多快，反映空间组织尺度。 |
| [`stage5a_stable_bond_strain_map.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_stable_bond_strain_map.png) | 49.5 KB | 稳定结构键应变空间图。 |
| [`stage5a_stable_disorder_map.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_stable_disorder_map.png) | 46.7 KB | 稳定结构 D 的空间图。 |
| [`stage5a_stable_physical_descriptors.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_stable_physical_descriptors.png) | 55.1 KB | 稳定结构各描述符的面内分布面板图。 |
| [`stage5a_stable_residual_energy_map.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_stable_residual_energy_map.png) | 44.9 KB | 控制 D 之后的残差能空间图。看热点在哪。 |
| [`stage5b_direction_similarity_distributions.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5b_direction_similarity_distributions.png) | 22.0 KB | 方向余弦的分布，按三个类别分开。 |
| [`stage5b_direction_vs_disorder.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5b_direction_vs_disorder.png) | 123.3 KB | 方向距离 vs D。 |
| [`stage5b_direction_vs_residual_energy.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5b_direction_vs_residual_energy.png) | 60.5 KB | 方向距离 vs 残差能。 |
| [`stage5b_direction_vs_stage5a_descriptors.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5b_direction_vs_stage5a_descriptors.png) | 49.4 KB | 方向距离 vs Stage5A 各描述符。**看方向是否只是描述符的重新组合**。 |
| [`stage5b_matched_energy_difference_vs_direction.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5b_matched_energy_difference_vs_direction.png) | 45.6 KB | 匹配对的能量差 vs 方向距离。 |
| [`stage5b_metastable_direction_autocorrelation.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5b_metastable_direction_autocorrelation.png) | 25.0 KB | 亚稳状态的方向自相关。 |
| [`stage5b_pca_direction.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5b_pca_direction.png) | 33.0 KB | 单位方向的 PCA 投影。PC1 / PC2 分别约 44.6% / 36.7%。 |
| [`stage5b_pca_displacement.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5b_pca_displacement.png) | 30.9 KB | 原始位移的 PCA 投影。 |
| [`stage5b_pca_displacement_3d.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5b_pca_displacement_3d.png) | 37.4 KB | 同上，三维版本。 |
| [`stage5b_representative_cases.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5b_representative_cases.png) | 35.8 KB | 8 个代表案例的图示。 |
| [`stage5b_stable_direction_autocorrelation.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5b_stable_direction_autocorrelation.png) | 25.0 KB | 稳定状态的方向自相关。 |

## 五、输出：结构文件

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`stage5a_matched_neighborhoods.xyz`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5a_matched_neighborhoods.xyz) | 1.64 MB | 匹配对的局域环境导出（63 个原子），OVITO 逐对查看。 |
| [`stage5b_representative_direction_cases.xyz`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5b_representative_direction_cases.xyz) | 117.4 KB | 代表案例的结构导出，用于人工核对。 |

## 六、未入库的产物

| 文件 | 大小 | 作用 |
|---|---:|---|
| `stage5b_matched_direction_energy_pairs.csv` | 142.87 MB | 全部匹配对，33 万行，142.9 MB，占整个项目近一半体积。最终结论只用由它算出的相关系数。 |

> 未入库项登记在 `DATA_MANIFEST.md`，附 SHA256 前 16 位，用于换机器后核对。

## 七、报告

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`Stage05_确定性结论与待讨论.md`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/01_%E6%8A%A5%E5%91%8A/Stage05_%E7%A1%AE%E5%AE%9A%E6%80%A7%E7%BB%93%E8%AE%BA%E4%B8%8E%E5%BE%85%E8%AE%A8%E8%AE%BA.md) | 3.8 KB | 物理描述符、残差关系、空间组织、方向定义与结果四节，附各描述符与残差能的 Spearman。 |

## 八、本阶段的关键数字

| 项 | 值 | 出处 |
|---|---|---|
| 物理截断 | 2.9924 Å（第一近邻） | 体相平均近邻距离 2.4793 Å，Voronoi 参考体积 10.7757 Å³ |
| 残差能定义 | `dE − 同 D 分箱中位数` | 主 18 箱，敏感 12 / 15 / 20 箱（Spearman 0.823 / 0.808 / 0.858） |
| 描述符 vs 残差 | 近邻距离 −0.503 / 键应变 −0.503 / CN −0.392 / 超额体积 −0.360 / CSP +0.202 | 界面原子 Spearman |
| 匹配相关 | 全体 0.8317 | stable-stable 0.2591 / meta-meta 0.6045 / stable-meta 0.7530 |
| 置换零模型 | 均值 0.0874 | 实测远高于零模型，说明方向信号非随机 |
| 方向 PCA | PC1 44.6% / PC2 36.7% | 描述性投影，不是独立验证 |
| 方向自相关 | meta 0.5235 / stable 0.6259 | 状态内的平均方向余弦 |

## 九、本阶段交给下游的东西

- **物理描述符基线**：Stage06 的每个模型都要拿它当对照，任何方向增益必须跑赢它才算数。
- **`stage5b_atomwise_soap_direction.csv`** 是 Stage06 的主变量表。
- **一条边界**：当前方向由原始位移构造，混着范数效应；PCA 在全体数据上拟合，属于预变换泄漏。

## 十、归档缺口

`stage5a` 与 `stage5b` 各写一份 `stage*_report.md`，归档里没有。5A 会重新读两个 CFG 算几何量，5B 依赖 5A 的输出表。
