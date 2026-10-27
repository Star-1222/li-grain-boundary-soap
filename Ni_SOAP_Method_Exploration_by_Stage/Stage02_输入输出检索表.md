# Stage 02：SOAP 表示与度量审计　输入・输出检索表

> 仓库 [`Star-1222/li-grain-boundary-soap`](https://github.com/Star-1222/li-grain-boundary-soap)　·　目录 `Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP表示与审计/`
> 本文件把 Stage 02 单独拎出来，逐个文件给出作用。链接指向 `main` 分支。

## 这一阶段在做什么

沿三个轴各切 200 箱找界面在哪里，取 33 个代表原子把四种距离的关系厘清，再把全人口 4659 个原子的 SOAP 矩阵算出来，量化「原始距离里有多少来自范数差」。

## 数据流

```
ni_stable.cfg ─┐
               ├─→ stage2_profiles.py ─→ spatial_profile_summary.txt + 2 图（定界面法向 = z）
ni_metastable.cfg ─┘        │
                            ├─→ stage2_soap.py ─→ 33 原子 ─→ soap_similarity_matrix.csv
                            │                        └─→ stage2_5_audit.py ─→ 528 对四度量
                            └─→ stage2_75_full.py ─→ 4659 × 1014 主矩阵 + 元数据 + 范数分解
```

## 一、输入

本阶段没有外部数据，吃的是上游产物。

| 文件 | 从哪来 | 作用 |
|---|---|---|
| [`ni_stable.cfg`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage01_%E6%95%B0%E6%8D%AE%E5%87%86%E5%A4%87%E4%B8%8ESOAP%E5%88%9D%E6%B5%8B/02_%E8%84%9A%E6%9C%AC/ni_stable.cfg) | Stage01 | 稳定态结构，本阶段的全部结构信息都来自它。 |
| [`ni_metastable.cfg`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage01_%E6%95%B0%E6%8D%AE%E5%87%86%E5%A4%87%E4%B8%8ESOAP%E5%88%9D%E6%B5%8B/02_%E8%84%9A%E6%9C%AC/ni_metastable.cfg) | Stage01 | 亚稳态结构，与上面共用同一套解析与抽样逻辑。 |
| [`cfg_reader.py`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage01_%E6%95%B0%E6%8D%AE%E5%87%86%E5%A4%87%E4%B8%8ESOAP%E5%88%9D%E6%B5%8B/02_%E8%84%9A%E6%9C%AC/cfg_reader.py) | Stage01 | 四个脚本都 `import cfg_reader`，用它读 CFG。 |

> 注意：这些脚本的 `OUTPUT_DIR` 都是「脚本所在目录 / outputs」，而且彼此用裸模块名 `import`，
> 运行时所有脚本与数据必须放在同一个目录里。归档时按阶段拆成了 `03_输出/`，
> 所以重跑之前要先把上游文件复制到本阶段的 `02_脚本/` 下。

## 二、处理代码

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`stage2_5_audit.py`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/02_%E8%84%9A%E6%9C%AC/stage2_5_audit.py) | 26.4 KB | 度量审计的核心：33 个原子共 528 对，同时算 `d_raw`、`s_raw`、`cos`、`d_norm`，并做 4 组 SOAP 参数敏感性。 |
| [`stage2_75_full.py`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/02_%E8%84%9A%E6%9C%AC/stage2_75_full.py) | 31.4 KB | 全人口 4659 个原子的 SOAP 主矩阵（4659 × 1014）与派生量，并把原始距离拆成范数部分和角度部分。 |
| [`stage2_profiles.py`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/02_%E8%84%9A%E6%9C%AC/stage2_profiles.py) | 8.9 KB | 沿 x / y / z 各切 200 箱做结构剖面，定位晶界法向。产出 `spatial_profile_summary.txt` 与两张分布图。 |
| [`stage2_soap.py`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/02_%E8%84%9A%E6%9C%AC/stage2_soap.py) | 14.1 KB | 从两个 CFG 里挑 33 个代表原子（体相、非主导配位、高能界面、界面成对），两两算 SOAP 余弦。产出相似度矩阵与热图。 |

## 三、输出：数据表

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`representative_atoms.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/representative_atoms.csv) | 4.8 KB | 33 个代表原子清单。`selection_type` 分四类，`selection_reason` 写明为什么选它。 |
| [`soap_similarity_matrix.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/soap_similarity_matrix.csv) | 11.4 KB | 33 × 33 余弦相似度矩阵，行列标签是「结构 / 原子id / 类别」。对角恒为 1，非对角都在 0.99 以上。 |
| [`spatial_profile_summary.txt`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/spatial_profile_summary.txt) | 1.3 KB | 200 箱剖面结果。**先看三处**：三个轴各识别出几个界面箱（只有 z 有）、界面原子数 1600 / 2059、界面与体相的 PEATOM 中位数差。 |
| [`stage2_5_atom_metadata.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/stage2_5_atom_metadata.csv) | 4.4 KB | 33 个代表原子的坐标、PTM / CNA 标签、PEATOM 与选取类别。 |
| [`stage2_5_pair_category_summary.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/stage2_5_pair_category_summary.csv) | 4.0 KB | 把 528 对按类别与度量分组后的 min / median / mean / std。 |
| [`stage2_5_pairwise_metrics.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/stage2_5_pairwise_metrics.csv) | 104.2 KB | 528 个原子对的四种距离。**四列的关系是本阶段的核心**：`d_norm² = 2 − 2·cos`，`s_raw ≠ cos`，`consistency_residual` 是等式残差。 |
| [`stage2_5_parameter_sensitivity.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/stage2_5_parameter_sensitivity.csv) | 254 B | 换 `n_max/l_max`（6/6、8/8、10/10）与 12/12 比，给出原始距离排序的 Spearman。**这就是 12/12 病态警告的后果量化**。 |
| [`stage2_75_atom_metadata.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/stage2_75_atom_metadata.csv) | 1.17 MB | 全人口 4659 行 × 19 列的主元数据。派生列 `soap_norm`、`d_mag_sq`、`d_ang_sq`、`f_mag`、`f_ang` 是全项目最该看懂的一组。 |
| [`stage2_75_pairwise_summary.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/stage2_75_pairwise_summary.csv) | 4.5 KB | 原子对距离的分位数汇总，支撑「原始距离约 42–46% 来自范数差」。 |
| [`stage2_75_population_summary.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/stage2_75_population_summary.csv) | 7.7 KB | 按 `category × metric` 汇总的分位数表（min / p5 / p25 / median / p75 / p95 / max）。 |
| [`stage2_75_ptm_cna_summary.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/stage2_75_ptm_cna_summary.csv) | 2.8 KB | 按 PTM / CNA 标签分组的 SOAP 范数、体相距离、PEATOM 中位数与 IQR。 |

## 四、输出：图

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`soap_similarity_heatmap.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/soap_similarity_heatmap.png) | 181.7 KB | 33 × 33 余弦矩阵的热图。看有没有按类别分块。 |
| [`spatial_profile_metastable.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/spatial_profile_metastable.png) | 121.9 KB | 同上，亚稳结构版本。两个结构都要单独看，不能只看一张。 |
| [`spatial_profile_stable.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/spatial_profile_stable.png) | 109.9 KB | 稳定结构的 x / y / z 三轴剖面。只有 z 出现尖峰，所以法向是 z。 |
| [`stage2_5_metric_heatmaps.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/stage2_5_metric_heatmaps.png) | 218.3 KB | 四种距离两两对照的热图，用来看它们何时一致、何时分叉。 |
| [`stage2_75_bulk_distance_distributions.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/stage2_75_bulk_distance_distributions.png) | 35.4 KB | 体相参考内部的原子间距离分布，用来判断 500 个抽样是否够代表。 |
| [`stage2_75_distance_decomposition.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/stage2_75_distance_decomposition.png) | 161.5 KB | **本阶段最重要的一张图**：平方距离里范数贡献与角度贡献的分解。 |
| [`stage2_75_distance_vs_PEATOM.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/stage2_75_distance_vs_PEATOM.png) | 149.4 KB | 到体相的距离与局域能量的散点，看有没有单调关系。 |
| [`stage2_75_distance_vs_z.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/stage2_75_distance_vs_z.png) | 111.7 KB | 距离随 z 坐标的分布，用来确认界面位置与无序带宽度。 |
| [`stage2_75_norm_distributions.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/03_%E8%BE%93%E5%87%BA/stage2_75_norm_distributions.png) | 29.2 KB | SOAP 向量自身的范数分布。体相约 77.4，界面约 72。**范数低不等于无序**，这一点必须记住。 |

## 五、输出：结构文件

本阶段没有结构文件输出。

## 六、未入库的产物

| 文件 | 大小 | 作用 |
|---|---:|---|
| `stage2_5_raw_soap_vectors.npy` | 261.5 KB | 33 个代表原子的原始 SOAP 向量，261 KB，是 528 对度量的原料。 |
| `stage2_75_bulk_reference_vectors.npy` | 16.0 KB | 两个结构各自的体相平均向量（16 KB），用来算 `d_to_bulk_mean` 与 `cos_to_bulk_mean`。 |
| `stage2_75_pairwise_sample.csv` | 48.78 MB | 70 万对原子对的抽样表，48.8 MB。支撑「42–46% 来自范数差」这个结论，本体可重算。 |
| `stage2_75_soap_vectors.npy` | 36.04 MB | SOAP 主矩阵，4659 × 1014，float64，36 MB。**Stage03–06 全部建立在这个文件上**，没它就是断链。 |

> 未入库项登记在 `DATA_MANIFEST.md`，附 SHA256 前 16 位，用于换机器后核对。

## 七、报告

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`Stage02_确定性结论与待讨论.md`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage02_SOAP%E8%A1%A8%E7%A4%BA%E4%B8%8E%E5%AE%A1%E8%AE%A1/01_%E6%8A%A5%E5%91%8A/Stage02_%E7%A1%AE%E5%AE%9A%E6%80%A7%E7%BB%93%E8%AE%BA%E4%B8%8E%E5%BE%85%E8%AE%A8%E8%AE%BA.md) | 4.0 KB | 本阶段能确定的 6 条（界面定位、代表原子、距离数学关系、参数敏感性、全人口产物、原始距离分解）、不能确定 8 条、待讨论 6 个问题、当前处置。 |

## 八、本阶段的关键数字

| 项 | 值 | 出处 |
|---|---|---|
| SOAP 人口 | 4659 = 1600 + 2059 + 500 + 500 | 界面原子按 z 掩膜取，体相各抽 500 |
| 界面掩膜 | 非主导 CNA/PTM 比例 > 0.1 且箱内 ≥ 50 原子 | 三个轴都试过，只有 z 出现界面箱 |
| SOAP 维度 | 1014 | 12×13/2 个径向对 × 13 个角层 |
| 范数差贡献 | 42–46% | `f_mag` 中位数，界面平方距离里来自范数的部分 |
| SOAP 范数中位数 | 体相 77.405 / 77.417，界面 72.006 / 72.737 | 界面范数更低，但**不能当成更无序** |
| 距离等式 | `d_norm² = 2 − 2·cos`，最大残差 1.5e-15 | `stage2_5_pairwise_metrics.csv` 的 `consistency_residual` |
| 参数敏感性 | 6/6 → 0.923，8/8 → 0.932，10/10 → 0.979 | 与 12/12 的距离排序 Spearman |

## 九、本阶段交给下游的东西

- **主矩阵**：`stage2_75_soap_vectors.npy`（未入库）是 Stage03–06 的共同地基，丢了要重算。
- **元数据**：4659 行逐原子表，后续每个阶段都要按 `row_index` 对齐。
- **一条结论**：`d_raw` 和 `d_norm` 是不同量，报告里必须分开命名。这条一直管到 Stage06。

## 十、归档缺口

四个脚本里的 `stage2_5_audit.py` 与 `stage2_75_full.py` 还会各写一份 `stage2_5_report.md` / `stage2_75_report.md`，但归档里没有这两个文件。
