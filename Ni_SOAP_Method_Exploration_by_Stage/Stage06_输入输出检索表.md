# Stage 06：最终评估与独立性检查　输入・输出检索表

> 仓库 [`Star-1222/li-grain-boundary-soap`](https://github.com/Star-1222/li-grain-boundary-soap)　·　目录 `Ni_SOAP_Method_Exploration_by_Stage/Stage06_最终评估/`
> 本文件把 Stage 06 单独拎出来，逐个文件给出作用。链接指向 `main` 分支。

## 这一阶段在做什么

把方向信号放进最严格的控制里：匹配 D、控制物理描述符、分状态、留一界面交叉验证、排除极端值。看它还剩多少。结论是方向没有独立增量。

## 数据流

```
stage4b_atomwise_energy_disorder.csv ─┐
stage5b_atomwise_soap_direction.csv ───┼─→ stage6_final_analysis.py
stage2_75_soap_vectors.npy ────────────┤      ├─→ 相关矩阵 + 偏相关（324 行）
stage5b_directional_statistics.csv ────┘      ├─→ 分组 CV 能量模型 / 状态分类
                                              ├─→ 极端值稳健性
                                              └─→ 代表案例打分
```

## 一、输入

本阶段没有外部数据，吃的是上游产物。

| 文件 | 从哪来 | 作用 |
|---|---|---|
| [`stage4b_atomwise_energy_disorder.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage04_%E6%97%A0%E5%BA%8F%E5%BA%A6%E4%B8%8E%E8%83%BD%E9%87%8F%E5%85%B3%E8%81%94/03_%E8%BE%93%E5%87%BA/stage4b_atomwise_energy_disorder.csv) | Stage04 | D 与过剩能。 |
| [`stage5b_atomwise_soap_direction.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5b_atomwise_soap_direction.csv) | Stage05 | 方向与位移，Stage06 的主变量。 |
| [`stage5b_directional_statistics.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/03_%E8%BE%93%E5%87%BA/stage5b_directional_statistics.csv) | Stage05 | 画图用的方向余弦统计。 |
| `stage2_75_soap_vectors.npy` | Stage02（未入库） | 第 296 行 `np.load`，分组 CV 需要原始向量做预变换。 |
| [`stage3a_lae_threshold.py`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage03_LAE%E9%98%88%E5%80%BC%E4%B8%8ELAE_LER%E6%9E%84%E5%BB%BA/02_%E8%84%9A%E6%9C%AC/stage3a_lae_threshold.py) | Stage03 | **代码依赖**：`import stage3a_lae_threshold as s3a`。 |
| [`stage5b_soap_direction.py`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage05_%E7%89%A9%E7%90%86%E7%A9%BA%E9%97%B4%E4%B8%8E%E6%96%B9%E5%90%91/02_%E8%84%9A%E6%9C%AC/stage5b_soap_direction.py) | Stage05 | **代码依赖**：`import stage5b_soap_direction as s5b`。 |

> 注意：这些脚本的 `OUTPUT_DIR` 都是「脚本所在目录 / outputs」，而且彼此用裸模块名 `import`，
> 运行时所有脚本与数据必须放在同一个目录里。归档时按阶段拆成了 `03_输出/`，
> 所以重跑之前要先把上游文件复制到本阶段的 `02_脚本/` 下。

## 二、处理代码

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`stage6_final_analysis.py`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage06_%E6%9C%80%E7%BB%88%E8%AF%84%E4%BC%B0/02_%E8%84%9A%E6%9C%AC/stage6_final_analysis.py) | 37.4 KB | 最终审计：相关矩阵、偏相关、匹配分析 + bootstrap + 置换、分组 CV 的能量模型与状态分类、极端值稳健性、空间汇总、代表案例打分。 |

## 三、输出：数据表

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`stage6_correlations.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage06_%E6%9C%80%E7%BB%88%E8%AF%84%E4%BC%B0/03_%E8%BE%93%E5%87%BA/stage6_correlations.csv) | 24.7 KB | 324 行 × 9 列。`scope / var_a / var_b / pearson / spearman / x / y / controls / partial_spearman`。**第一行就是 D 与位移范数 0.99998**。 |
| [`stage6_direction_by_disorder_bin.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage06_%E6%9C%80%E7%BB%88%E8%AF%84%E4%BC%B0/03_%E8%BE%93%E5%87%BA/stage6_direction_by_disorder_bin.csv) | 1.4 KB | 20 行 = 状态 × D 分箱，看每个箱内的方向一致性。 |
| [`stage6_energy_models.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage06_%E6%9C%80%E7%BB%88%E8%AF%84%E4%BC%B0/03_%E8%BE%93%E5%87%BA/stage6_energy_models.csv) | 732 B | 8 行 = 4 个模型 × 2 折。R² 0.0393 / 0.7497 / 0.9869 / 0.9886（两折平均）。 |
| [`stage6_figure_summary.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage06_%E6%9C%80%E7%BB%88%E8%AF%84%E4%BC%B0/03_%E8%BE%93%E5%87%BA/stage6_figure_summary.csv) | 536 B | 6 张图的文件名与内容对照，纯索引表。 |
| [`stage6_matched_analysis.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage06_%E6%9C%80%E7%BB%88%E8%AF%84%E4%BC%B0/03_%E8%BE%93%E5%87%BA/stage6_matched_analysis.csv) | 3.4 KB | 16 行。匹配相关 + bootstrap 区间 + 置换检验 + 偏相关。0.8317 → 0.2058 就出自这里。 |
| [`stage6_representative_directional_cases.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage06_%E6%9C%80%E7%BB%88%E8%AF%84%E4%BC%B0/03_%E8%BE%93%E5%87%BA/stage6_representative_directional_cases.csv) | 7.8 KB | 12 个代表案例，带 z 分数与 Score I–IV 四套打分。 |
| [`stage6_robustness.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage06_%E6%9C%80%E7%BB%88%E8%AF%84%E4%BC%B0/03_%E8%BE%93%E5%87%BA/stage6_robustness.csv) | 211 B | 排除极端 CSP / 体积 / 能量后的匹配 Spearman，0.8313 / 0.8320 / 0.8316，说明结论对极端值不敏感。 |
| [`stage6_spatial_summary.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage06_%E6%9C%80%E7%BB%88%E8%AF%84%E4%BC%B0/03_%E8%BE%93%E5%87%BA/stage6_spatial_summary.csv) | 955 B | 空间场的汇总统计与 top-10% 热点近邻距离。 |
| [`stage6_state_classification.csv`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage06_%E6%9C%80%E7%BB%88%E8%AF%84%E4%BC%B0/03_%E8%BE%93%E5%87%BA/stage6_state_classification.csv) | 1.1 KB | 10 行 = 5 个模型 × 2 折。ROC-AUC 0.606 / 0.684 / 0.703 / 0.577 / 0.689。 |

## 四、输出：图

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`stage6_directional_heterogeneity.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage06_%E6%9C%80%E7%BB%88%E8%AF%84%E4%BC%B0/03_%E8%BE%93%E5%87%BA/stage6_directional_heterogeneity.png) | 21.7 KB | 状态内与跨状态的余弦分布对照。 |
| [`stage6_energy_model_improvement.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage06_%E6%9C%80%E7%BB%88%E8%AF%84%E4%BC%B0/03_%E8%BE%93%E5%87%BA/stage6_energy_model_improvement.png) | 20.5 KB | 四个模型的 R² / RMSE / MAE 对比。**加方向只涨 0.0017 就在这张图里**。 |
| [`stage6_magnitude_matched_direction_energy.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage06_%E6%9C%80%E7%BB%88%E8%AF%84%E4%BC%B0/03_%E8%BE%93%E5%87%BA/stage6_magnitude_matched_direction_energy.png) | 45.2 KB | 幅度匹配后的方向距离 vs 能量差。 |
| [`stage6_representative_cases.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage06_%E6%9C%80%E7%BB%88%E8%AF%84%E4%BC%B0/03_%E8%BE%93%E5%87%BA/stage6_representative_cases.png) | 35.7 KB | 12 个方向性代表案例。 |
| [`stage6_scalar_disorder_vs_energy.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage06_%E6%9C%80%E7%BB%88%E8%AF%84%E4%BC%B0/03_%E8%BE%93%E5%87%BA/stage6_scalar_disorder_vs_energy.png) | 28.5 KB | 界面 D 与过剩能的关系（按状态着色）。 |
| [`stage6_spatial_direction_residual.png`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage06_%E6%9C%80%E7%BB%88%E8%AF%84%E4%BC%B0/03_%E8%BE%93%E5%87%BA/stage6_spatial_direction_residual.png) | 80.3 KB | 方向异质性与残差能的空间图并排。 |

## 五、输出：结构文件

本阶段没有结构文件输出。

## 六、未入库的产物

本阶段全部产物都已入库。

> 未入库项登记在 `DATA_MANIFEST.md`，附 SHA256 前 16 位，用于换机器后核对。

## 七、报告

| 文件 | 大小 | 作用 |
|---|---:|---|
| [`Stage06_确定性结论与待讨论.md`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/Ni_SOAP_Method_Exploration_by_Stage/Stage06_%E6%9C%80%E7%BB%88%E8%AF%84%E4%BC%B0/01_%E6%8A%A5%E5%91%8A/Stage06_%E7%A1%AE%E5%AE%9A%E6%80%A7%E7%BB%93%E8%AE%BA%E4%B8%8E%E5%BE%85%E8%AE%A8%E8%AE%BA.md) | 3.5 KB | 匹配分析的实际结果、分组 CV 的模型表现、稳健性排除、本次补充核对四节，是全项目方法结论的落点。 |

## 八、本阶段的关键数字

| 项 | 值 | 出处 |
|---|---|---|
| 匹配相关 | 0.8317 → 控制后 0.2058 | **全项目最关键的一次衰减** |
| 分状态 | 0.2591 → 0.0519 / 0.6045 → 0.1823 / 0.7530 → 0.0907 | stable-stable 几乎归零 |
| 能量模型 R² | D 0.0393 / D+方向 0.7497 / D+物理 0.9869 / 全加 0.9886 | 两折平均 |
| 状态分类 ROC-AUC | D 0.606 / 物理 0.684 / D+物理 0.703 / D+方向 0.577 / 全加 0.689 | 加方向反而变差 |
| 极端值稳健性 | 0.8313 / 0.8320 / 0.8316 | 排除极端 CSP / 体积 / 能量后 |
| 最大疑点 | D 与位移范数相关 0.99998 | `stage6_correlations.csv` 第一行 |

## 九、本阶段交给下游的东西

- **方法线的判决**：物理描述符是强基线，方向没有独立于它的预测增量。
- **可以说的**：匹配 D 时方向与能量差存在描述性关联，且非随机。
- **不能说的**：方向提供大幅预测增益、方向导致能量变化、结论适用于所有晶界。

## 十、归档缺口

`stage6_final_analysis.py` 是唯一不写 `*_report.md` 的脚本之一，但它的图表有 `stage6_figure_summary.csv` 做索引。全部结论只限当前两个 Ni 构型。
