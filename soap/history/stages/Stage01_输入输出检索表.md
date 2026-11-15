# Stage 01 输入・输出检索表

> 仓库 [`Star-1222/li-grain-boundary-soap`](https://github.com/Star-1222/li-grain-boundary-soap)　·　目录 `Ni_SOAP_Method_Exploration_by_Stage/Stage01_数据准备与SOAP初测/`
> 本文件把 Stage 01 单独拎出来，逐个文件给出作用。链接均指向 `main` 分支的文件页。

## 这一阶段在做什么

一句话：**确认两个 Ni CFG 能被可靠解析，SOAP 能在这台机器上跑通。**

它只检查数据与工具链的可信性，不评价势函数、能量来源或稳定/亚稳的物理含义。
判定通过不等于物理结论成立，只是拿到进入 Stage 02 的资格。

## 数据流

```
ni_stable.cfg ─┐
               ├─→ cfg_reader.py
ni_metastable.cfg ─┘        │
                            ├─→ inspect_structures.py ─→ inspection_summary.txt
                            ├─→ soap_test.py ─→ soap_test_summary.txt
                            │                └─→ soap_vectors_*.npz
                            └─→ stage1_qc.py ─→ stage1_qc_report.txt（A–F，Overall PASS）
```

## 一、输入

| 文件 | 大小 | 角色与作用 |
|---|---:|---|
| [`ni_stable.cfg`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage01_%E6%95%B0%E6%8D%AE%E5%87%86%E5%A4%87%E4%B8%8ESOAP%E5%88%9D%E6%B5%8B/02_%E8%84%9A%E6%9C%AC/ni_stable.cfg) | 4.21 MB | 稳定态 Ni 双晶构型。47,120 原子，`entry_count=11`（3 列坐标 + 8 个辅助字段），晶胞 101.068 × 17.5399 × 286.870 Å。外部给的原始数据，不是本项目生成。 |
| [`ni_metastable.cfg`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage01_%E6%95%B0%E6%8D%AE%E5%87%86%E5%A4%87%E4%B8%8ESOAP%E5%88%9D%E6%B5%8B/02_%E8%84%9A%E6%9C%AC/ni_metastable.cfg) | 4.21 MB | 亚稳态 Ni 双晶构型，同样 47,120 原子，晶胞 101.067 × 17.5397 × 286.866 Å。与 stable 对照，用来检验方法能否把两个状态区分开。 |

两个 CFG 都是**外部给的原始数据**，本仓库不生成它们。文件 4.2 MB 的原因是有 47,120 行原子记录。

## 二、处理代码

| 文件 | 大小 | 角色与作用 |
|---|---:|---|
| [`cfg_reader.py`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage01_%E6%95%B0%E6%8D%AE%E5%87%86%E5%A4%87%E4%B8%8ESOAP%E5%88%9D%E6%B5%8B/02_%E8%84%9A%E6%9C%AC/cfg_reader.py) | 8.7 KB | CFG 解析模块，无 `main`。按头部字段名读取，不依赖固定行号；自动判断 scaled / cartesian 坐标。被下面三个脚本导入。 |
| [`inspect_structures.py`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage01_%E6%95%B0%E6%8D%AE%E5%87%86%E5%A4%87%E4%B8%8ESOAP%E5%88%9D%E6%B5%8B/02_%E8%84%9A%E6%9C%AC/inspect_structures.py) | 4.2 KB | 逐字段统计并写出 `inspection_summary.txt`。给“数据到底长什么样”的第一手描述。 |
| [`soap_test.py`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage01_%E6%95%B0%E6%8D%AE%E5%87%86%E5%A4%87%E4%B8%8ESOAP%E5%88%9D%E6%B5%8B/02_%E8%84%9A%E6%9C%AC/soap_test.py) | 5.0 KB | 取前 5 个原子跑 DScribe SOAP，检查维度、有限性、可复现性，写出 `soap_test_summary.txt` 与两个 npz。 |
| [`stage1_qc.py`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage01_%E6%95%B0%E6%8D%AE%E5%87%86%E5%A4%87%E4%B8%8ESOAP%E5%88%9D%E6%B5%8B/02_%E8%84%9A%E6%9C%AC/stage1_qc.py) | 21.2 KB | 六项质检 A–F 并写出 `stage1_qc_report.txt`：CFG 完整性、数值完整性、坐标与晶胞一致性、PBC 一致性、SOAP 冒烟复用、PBC 邻居合理性。 |

`stage1_qc.py` 是唯一的汇总入口，它会把 `soap_test_summary.txt` 直接复用，不重跑 SOAP。

## 三、输出

| 文件 | 大小 | 角色与作用 |
|---|---:|---|
| [`inspection_summary.txt`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage01_%E6%95%B0%E6%8D%AE%E5%87%86%E5%A4%87%E4%B8%8ESOAP%E5%88%9D%E6%B5%8B/03_%E8%BE%93%E5%87%BA/inspection_summary.txt) | 2.8 KB | 两个 CFG 的解析摘要。**先看三处**：`header atom count` 与 `parsed atom count` 是否都等于 47120；8 个辅助字段的 min/max/unique；`source` 里记录的原始路径。 |
| [`soap_test_summary.txt`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage01_%E6%95%B0%E6%8D%AE%E5%87%86%E5%A4%87%E4%B8%8ESOAP%E5%88%9D%E6%B5%8B/03_%E8%BE%93%E5%87%BA/soap_test_summary.txt) | 1.4 KB | 5 原子 SOAP 冒烟测试。**看四行**：`SOAP matrix shape: (5, 1014)`、`all finite: True`、`reproducible: True`、非对角余弦区间 0.9914–1.0000。 |
| [`stage1_qc_report.txt`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage01_%E6%95%B0%E6%8D%AE%E5%87%86%E5%A4%87%E4%B8%8ESOAP%E5%88%9D%E6%B5%8B/03_%E8%BE%93%E5%87%BA/stage1_qc_report.txt) | 7.3 KB | A–F 六节质检，每节末尾 `PASS`。最后 `[Overall] overall: PASS` 就是本阶段的总判定。末尾 `[Notes]` 列了三条已知边界（病态矩阵警告、坐标约定靠启发式判断、复用旧结果）。 |

三份 txt 的分工：`inspection` 说数据长什么样，`soap_test` 说 SOAP 能不能算，
`stage1_qc` 把两者收进 A–F 六项判定。

## 四、报告

| 文件 | 大小 | 角色与作用 |
|---|---:|---|
| [`Stage01_确定性结论与待讨论.md`](https://github.com/Star-1222/li-grain-boundary-soap/blob/main/soap/Stage01_%E6%95%B0%E6%8D%AE%E5%87%86%E5%A4%87%E4%B8%8ESOAP%E5%88%9D%E6%B5%8B/01_%E6%8A%A5%E5%91%8A/Stage01_%E7%A1%AE%E5%AE%9A%E6%80%A7%E7%BB%93%E8%AE%BA%E4%B8%8E%E5%BE%85%E8%AE%A8%E8%AE%BA.md) | 3.2 KB | 本阶段能确定的 4 条、不能确定的 8 条、待讨论的 4 个问题、当前处置。想知道“这一步到底证明了什么”，读这份，不要读 QC 报告。 |

## 五、未入库的产物

| 文件 | 大小 | 角色与作用 |
|---|---:|---|
| `soap_vectors_stable.npz` | 39 KB | `soap_test.py` 存下的 5 个原子 SOAP 向量，附带 rcut / nmax / lmax / sigma。 |
| `soap_vectors_metastable.npz` | 39 KB | 同上，亚稳态版本。 |

这两个 npz 只有 39 KB，但当时按“中间产物”处理没进版本控制，属 `DATA_MANIFEST.md` 的登记范围。
重算代价很小，`soap_test.py` 跑一次即可重建。

## 六、本阶段的关键数字

| 数字 | 值 | 出处 | 含义 |
|---|---|---|---|
| 原子数 | 47,120（两个文件都一致） | `inspection_summary.txt` | 头部声明与实际解析相等 |
| 字段数 | 11 = 3 坐标 + 8 辅助 | 同上 | 辅助字段为 id / type / PTM / CNA / ENT1-3 / PEATOM |
| SOAP 维度 | 1014 | `soap_test_summary.txt` | 单元素、`n_max=12, l_max=12` 下 12×13/2×13 |
| SOAP 参数 | r_cut 3.25 Å / n_max 12 / l_max 12 / sigma 0.5 Å | 同上 | 只是初测值，本阶段不冻结 |
| 稳定性 | 5 原子矩阵全部有限、可复现 | 同上 | 同一输入两次结果一致 |
| 总判定 | PASS × 2（stable + metastable） | `stage1_qc_report.txt` 末尾 | 拿到进入 Stage 02 的资格 |

## 七、本阶段交给下游的东西

- **可信输入**：两个 CFG 可以被同一套解析器读出，坐标约定统一为 scaled。
- **参数初值**：`r_cut=3.25 / n_max=12 / l_max=12 / sigma=0.5`，Stage 02–06 一直沿用，但始终标着“未冻结”。
- **一条已知缺陷**：12/12 触发 DScribe 病态矩阵警告，Stage 02 用参数敏感性（6/6、8/8、10/10）回应了它。

本阶段**没有**产出 4,659 这个人口数。那是 Stage 02 沿 z 轴定界面、再抽样之后才有的事。
