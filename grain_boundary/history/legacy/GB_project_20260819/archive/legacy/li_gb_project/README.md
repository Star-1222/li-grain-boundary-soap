# Li_GB_Project — 锂金属晶界模拟项目

> 大创项目：固态电池锂金属晶界亚稳态模拟  
> BJTU · 计算材料科学方向 · 2026.06 – 2027.06

---

## 项目简介

本项目利用 **SOAP 描述符 + 机器学习势函数 + LAMMPS MD 模拟**，研究 BCC 锂金属对称倾斜晶界（STGB）的亚稳态原子结构与力学行为，建立 **构型 → 机制 → 强度** 的映射关系。

核心技术栈：
- **Li_Bicrystal_Generator** — Python 包，基于 ASE 生成 Σ3/Σ5/Σ7 等 CSL 双晶模型
- **LAMMPS + OpenKIM (Li-SNAP)** — 分子动力学模拟引擎与势函数
- **SOAP** — 局域原子环境描述符，用于结构无序度表征

---

## 目录结构

```
Li_GB_Project/
│
├── Li_Bicrystal_Generator/       # 核心 Python 包：双晶模型生成器
│   ├── __init__.py               #   包初始化与公共 API
│   ├── __main__.py               #   CLI 入口：python -m Li_Bicrystal_Generator
│   ├── main.py                   #   Argparse 命令行逻辑
│   ├── config.py                 #   物理参数 + 晶界数据库（Σ5/Σ13/Σ17/Σ29）
│   ├── builder.py                #   BiCrystalBuilder 核心类
│   ├── geometry.py               #   旋转矩阵、原子去重、居中
│   ├── exporter.py               #   LAMMPS/XYZ/CIF/grainID 导出
│   ├── utils.py                  #   校验与晶界查找工具
│   ├── README.md                 #   包详细文档
│   ├── examples/example.py       #   示例脚本
│   ├── output/                   #   生成的双晶文件（.lmp/.xyz/.cif/_grainID.xyz）
│   └── lammps_workflow/          #   LAMMPS 模拟工作流
│       ├── in.bicrystal_relax    #     晶界弛豫输入脚本
│       ├── in.quench             #     淬火输入脚本
│       ├── in.resume / continue.in  #  续算脚本
│       ├── anneal_save.in        #     退火保存脚本
│       ├── dump2data.in          #     dump 转换脚本
│       ├── run_relax.sh          #     弛豫运行脚本
│       ├── Sigma5_310.lmp        #     Σ5(310) 初始构型
│       ├── Sigma5_310_minimized.data / .restart  #  最小化后构型
│       ├── annealed.data         #     退火后构型
│       ├── final_gb.data / .restart  #  最终构型
│       ├── traj_anneal.dump      #     退火轨迹
│       └── traj_quench.dump      #     淬火轨迹
│
├── lammps_test/                  # Li-SNAP + OpenKIM 环境验证
│   ├── in.test_li                #   BCC Li 2×2×2 + pair_style kim 验证脚本
│   ├── run.sh                    #   便捷运行脚本
│   ├── README.md                 #   环境配置与验证文档
│   ├── final.data                #   验证输出数据
│   └── log.lammps                #   LAMMPS 运行日志
│
├── docs/                         # 开题报告与演示
│   ├── Li_metal_metastable_atomic_structure_mechanical_analysis.pdf
│   └── Li_metal_metastable_GB_atomic_mechanical_analysis.pptx
│
├── literature/                   # 参考文献全文
│   ├── 2024_Metastable grain boundaries...pdf
│   ├── Atomistic modeling of the mechanical.pdf
│   ├── Diffusion and Deformation Mechanism Maps in Li.pdf
│   ├── Elastic and plastic mechanical properties of lithium.pdf
│   ├── Grain Boundary Softening - Lithium Metal.pdf
│   ├── mechanistic-insights-into-li-transport...pdf
│   ├── The growth mechanism of lithium dendrites.pdf
│   └── soap_method/              #   SOAP 方法专题
│       ├── SOAP_structural_disorder_characterization.pdf
│       └── SOAP_structural_disorder_sources.pdf
│
├── notes/                        # 工作笔记（待补充）
├── README.md                     # 本文件
└── .gitignore
```

---

## 快速开始

### 1. 安装依赖

```bash
pip install ase numpy
```

### 2. 生成双晶模型

```bash
cd Li_Bicrystal_Generator

# 列出可用晶界
python -m Li_Bicrystal_Generator --list

# 生成 Σ5(310) 双晶
python -m Li_Bicrystal_Generator --name Sigma5_310
```

输出文件出现在 `Li_Bicrystal_Generator/output/`。

### 3. LAMMPS 模拟

```bash
# 环境验证（确保 Li-SNAP 正确加载）
cd ../lammps_test
bash run.sh

# 晶界弛豫
cd ../Li_Bicrystal_Generator/lammps_workflow
bash run_relax.sh
```

### 4. OVITO 可视化

- 打开 `.xyz` 或 `.lmp` 文件
- 用 `_grainID.xyz` 文件按晶粒着色（grain_id: 1=晶粒1, 2=晶粒2）
- GB 面垂直于 x 轴，位于 cell 中心

---

## 晶界数据库

| 名称 | Σ | 旋转轴 | 角度 (°) | 晶界面 |
|------|---|--------|----------|--------|
| Sigma5_310 | 5 | [001] | 36.87 | (310) |
| Sigma5_210 | 5 | [001] | 53.13 | (210) |
| Sigma13_510 | 13 | [001] | 22.62 | (510) |
| Sigma17_410 | 17 | [001] | 28.07 | (410) |
| Sigma17_530 | 17 | [001] | 61.93 | (530) |
| Sigma29_520 | 29 | [001] | 43.60 | (520) |

---

## 数据来源

本项目于 2026-07-31 由以下三个目录整合而成：

| 原目录 | 整合去向 |
|--------|---------|
| `Documents/Li_Project/` | Li_Bicrystal_Generator/ + lammps_test/ |
| `Desktop/SOAP_innovation_project/` | docs/ + literature/ |
| `Documents/Codex/` | 重复副本已丢弃，空目录已清除 |

原始备份位于 `Documents/_backup_Li_GB_20260731/`。
