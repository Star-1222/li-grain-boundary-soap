# BCC Li Sigma5 [001](310) Grain-Boundary Workflow (Bilingual / 双语版)

本文件为双语版：中文版在前，英文原版在后。This file contains the Chinese version first, followed by the full English original.

# 中文版

# BCC Li Σ5 [001](310) 晶界工作流（中文版）

面向周期 BCC 锂对称倾转晶界的已验证构建、退火、淬火与弛豫工作流。

## 项目状态

本仓库是在验证生成器、文件导出器、LAMMPS/OpenKIM 工作流、断点恢复和最终
弛豫结构之后产生的权威清洁包。

- 自动化回归测试：**63 项通过**
- 初始双晶验证：**通过**
- LAMMPS data、扩展 XYZ 与 CIF 重导入：**通过**
- dump 到 data 重建：**通过**
- Stage 0-4 生产工作流：**完成**
- 最终原子力判据：**通过**
- 最终法向压力判据：**通过**
- 危险近邻列表构建：**0**

最终结构是单条 anneal-quench-minimize 轨迹收敛得到的一个局域极小候选。
它本身不应被解释为亚稳晶界态的统计代表系综。

## 物理模型

| 参数 | 值 |
|---|---:|
| 材料 | BCC Li |
| 晶界类型 | 对称倾转晶界 |
| CSL 名称 | Σ5 |
| 旋转轴 | [001] |
| 晶界平面 | (310) |
| 取向差 | 36.86989765° |
| 初始晶格常数 | 3.506186475312 Å |
| 边界条件 | x、y、z 均周期 |
| 模拟盒子 | 正交 |
| 最终原子数 | 3600 |
| 每晶粒原子数 | 1800 |
| 原子类型 | `atomic` |
| Li 原子质量 | 6.941 |
| 相互作用模型 | `SNAP_ZuoChenLi_2019_Li__MO_732106099012_000` |
| 生产并行 | 12 MPI 任务 × 1 OpenMP 线程 |

由于沿晶界法向的周期性，盒子包含两条晶体学等价的晶界。

### 双晶构建

生成器构建晶粒 1，并通过晶界平面镜像得到晶粒 2。在两个周期晶界平面上独立
检测重合位点，每个重合对恰好删除一个原子，并在两个晶粒间平衡删除。

| 数量 | 值 |
|---|---:|
| 删除前每晶粒位点数 | 1860 |
| 删除前总原子数 | 3720 |
| GB1 重合对 | 60 |
| GB2 重合对 | 60 |
| 晶粒 1 删除 | 60 |
| 晶粒 2 删除 | 60 |
| 最终原子数 | 3600 |

初始验证盒子约为：

```text
66.5252 × 33.2626 × 35.0619 Å
```

验证的最小原子间距为 2.217507 Å。

## 相互作用势与体相参考

所有生产计算使用同一 OpenKIM SNAP 模型：

```lammps
pair_style kim SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
pair_coeff * * Li
```

使用同一模型对周期 128 原子 BCC 体相参考进行零压弛豫。

| 体相参考结果 | 值 |
|---|---:|
| 弛豫晶格常数 | 3.506186475312 Å |
| 每原子势能 | -1.89998153759969 eV/atom |

最终双晶构建使用体相弛豫晶格常数，使面内 CSL 周期不继承可避免的体相应变。

## 模拟工作流

| Stage | 操作 | 主要设置 | 主要输出 |
|---|---|---|---|
| 0 | 势/结构检查 | `run 0` | 成功的力与近邻构建 |
| 1 | 初始极小化 | CG；固定 y/z；弛豫 x；分块 checkpoint | `stage1_minimized.data` |
| 2 | 退火 | NVT，300 K，20 ps，20000 步 | `stage2_annealed.data` |
| 3 | 淬火 | NVT，300 -> 1 K，50 ps，50000 步 | `stage3_quenched.data` |
| 4 | 最终极小化 | 0 K CG；固定 y/z；弛豫 x；严格收敛 | `stage4_final.data` |

Stage 2 和 Stage 3 的时间步为 0.001 ps。

### 盒子弛豫策略

晶界平面为 yz 平面，x 为晶界法向。

- **y 和 z 固定**在零压体相/CSL 周期。
- **x 向零法向压力弛豫**。

这样既保持面内晶体学周期性，又允许双晶容纳晶界过剩体积。报告中的最终收敛
判据适用于 `Pxx`；由于 y 和 z 被刻意约束，非零 `Pyy` 和 `Pzz` 不独立极小化。

### 断点与恢复策略

- Stage 1 和 Stage 4 以可恢复的 500 迭代分块运行极小化。
- Stage 2 和 Stage 3 每 1000 步（1 ps）写 restart。
- Stage 2 和 Stage 3 有专用 resume 输入，保留原目标时间步和温度计划。
- 通用旧版 `in.resume` 已退役，因为单个通用分支无法安全保留各阶段不同的
  物理协议。

## 最终收敛结构

权威文件：

```text
outputs/stage4_final.data
```

| 最终结果 | 值 |
|---|---:|
| 原子数 | 3600 |
| 原子类型 | 1 |
| Lx | 66.01192624065516 Å |
| Ly | 33.2626054898 Å |
| Lz | 35.0618647531 Å |
| 最大原子力 | 3.55281121421774 × 10⁻⁷ eV/Å |
| 力范数 | 1.95219377952745 × 10⁻⁶ |
| Pxx | 3.63132085694783 × 10⁻⁴ bar |
| 最大绝对速度 | 0 |
| 危险构建 | 0 |

最终 data 文件包含有限坐标、严格零速度、3600 个原子和一种原子类型。

## 仓库布局

```text
v2_final/
├── generator/                 Python 双晶生成与验证代码
├── tests/                     自动化回归测试
├── lammps/                    有效 LAMMPS 脚本
│   ├── in.0_check
│   ├── in.1_minimize
│   ├── in.2_anneal
│   ├── in.2_anneal_resume
│   ├── in.3_quench
│   ├── in.3_quench_resume
│   ├── in.4_final
│   ├── in.bulk_energy
│   ├── in.dump2data
│   ├── run_all.sh
│   └── resume.sh
├── outputs/
│   ├── Sigma5_310.*           已验证的初始双晶导出
│   ├── bulk/                  体相参考结构与 restart
│   ├── stage1_minimized.*
│   ├── stage2_annealed.*
│   ├── stage2_anneal.dump
│   ├── stage3_quenched.*
│   ├── stage3_quench.dump
│   └── stage4_final.*         权威最终结构
├── logs/                      成功生产日志与完成标记
├── validation/                dump 到 data 重建证据
├── FINAL_STATUS.md            简明验证摘要
├── ENVIRONMENT.txt            软件/环境记录
├── SHA256SUMS.txt             文件完整性哈希
├── requirements.txt
└── pytest.ini
```

最终包有意排除旧的 a = 3.51 Å 结构、中断运行、过时脚本、语法测试输出、
已完成的 checkpoint 文件、PID 文件、Python 缓存和 `.bak` 文件。

## 验证随附包

激活 Python 环境后运行：

```bash
python -m pytest tests -q

python -m generator.verify \
  outputs/Sigma5_310.data \
  --map outputs/Sigma5_310.grain_map

python -m generator.reimport_check \
  --dir outputs \
  --prefix Sigma5_310

sha256sum -c SHA256SUMS.txt
```

预期高层结果：

```text
all tests passed
[verify] PASSED: 3600 atoms (g1=1800, g2=1800)
[reimport] ALL OK
all SHA-256 checks passed
```

## 在新目录中复现工作流

不要为新运行复用随附的已完成 `outputs/` 和 `logs/`。先生成新输入：

```bash
python -m generator.generate \
  --a 3.506186475312 \
  --outdir outputs_new

python -m generator.bulk \
  --a 3.506186475312 \
  --outdir outputs_new/bulk
```

如需要，运行体相参考：

```bash
OMP_NUM_THREADS=1 mpirun -np 12 lmp \
  -in lammps/in.bulk_energy \
  -var outdir "$PWD/outputs_new" \
  -var kim_model SNAP_ZuoChenLi_2019_Li__MO_732106099012_000 \
  -log logs_new/bulk_energy.lammps
```

运行 Stage 0-4：

```bash
OUTDIR="$PWD/outputs_new" \
LOGDIR="$PWD/logs_new" \
MPI_NP=12 \
OMP_NUM_THREADS=1 \
bash lammps/run_all.sh
```

提供的 runner 每个 stage 写一个日志，只有成功完成后才写 `.done` 标记。

## 恢复中断 stage

使用同一 `OUTDIR` 和 `LOGDIR` 的专用 launcher：

```bash
OUTDIR="$PWD/outputs_new" LOGDIR="$PWD/logs_new" \
  MPI_NP=12 bash lammps/resume.sh stage1

OUTDIR="$PWD/outputs_new" LOGDIR="$PWD/logs_new" \
  MPI_NP=12 bash lammps/resume.sh stage2

OUTDIR="$PWD/outputs_new" LOGDIR="$PWD/logs_new" \
  MPI_NP=12 bash lammps/resume.sh stage3

OUTDIR="$PWD/outputs_new" LOGDIR="$PWD/logs_new" \
  MPI_NP=12 bash lammps/resume.sh stage4
```

resume 输入选择如下：

| Stage | Resume 来源 |
|---|---|
| 1 | `stage1_checkpoint.data` |
| 2 | 最新 `stage2_anneal.checkpoint.*` |
| 3 | 最新 `stage3_quench.checkpoint.*` |
| 4 | `stage4_checkpoint.restart` |

## 从 dump 重建 data 文件

`in.dump2data` 用现有 data 文件作为盒子、质量和原子类型模板，再用 LAMMPS
原生 dump 帧替换坐标和速度：

```bash
lmp -in lammps/in.dump2data \
  -var template outputs/stage1_minimized.data \
  -var dump outputs/stage2_anneal.dump \
  -var timestep 20000 \
  -var out rebuilt_stage2.data
```

该路径已通过重建 timestep 20000 并与 `stage2_annealed.data` 比较原子数、
原子类型和盒子尺寸得到验证。

## 科学范围与局限性

本包为单一弛豫晶界态建立了可复现路径，但以下结论仍需额外工作：

1. **亚稳态采样：** 一条热轨迹不能映射可达局域极小态的总体。需要独立随机
   扰动、不同热历史或 activation-relaxation 搜索。
2. **不同能量盆地验证：** 应使用能量、SOAP 等局部原子描述符以及结构距离/
   去重判据比较结构。
3. **晶界能量：** 已包含匹配体相参考，但此处尚未报告最终 GB 能量分析与
   不确定度处理。
4. **势函数迁移性：** SNAP 模型在本工作流内部自洽，但在强定量结论前，其
   对 Li 晶界构型的精度应使用更高级参考数据评估。
5. **有限尺寸效应：** 周期盒含两条 GB；定量发表前应检查 GB 间距和面内面积
   的收敛性。

## 软件环境

已验证生产运行使用：

- WSL/Linux
- Python 3.14.4
- pytest 9.1.1
- LAMMPS 10 Dec 2025
- OpenKIM 支持
- OpenMPI，12 MPI ranks
- 每 MPI rank 一个 OpenMP 线程

精确记录的环境与运行引用见 `ENVIRONMENT.txt` 和 LAMMPS 日志。

## 引用

使用这些结果时，按 LAMMPS 日志要求引用 OpenKIM 框架和所选模型：

- OpenKIM Project: <https://doi.org/10.1007/s11837-011-0102-6>
- Model citation page: <https://openkim.org/cite/MO_732106099012_000#item-citation>

## 完整性与溯源

`SHA256SUMS.txt` 记录所有随附文件的完整性。成功的 Stage 0-4 日志保留在
`logs/`。更早的工作目录和过时 Windows 副本不属于权威科学包。


# English Original

# BCC Li Σ5 [001](310) Grain-Boundary Workflow

Verified construction, annealing, quenching and relaxation workflow for a periodic BCC lithium symmetric tilt grain boundary.

## Project status

This repository is the authoritative clean package produced after validating the generator, file exporters, LAMMPS/OpenKIM workflow, checkpoint recovery and final relaxed structure.

- Automated regression tests: **63 passed**
- Initial bicrystal verification: **passed**
- LAMMPS-data, extended-XYZ and CIF re-import: **passed**
- Dump-to-data reconstruction: **passed**
- Stage 0–4 production workflow: **completed**
- Final atomic-force criterion: **passed**
- Final normal-pressure criterion: **passed**
- Dangerous neighbor-list builds: **0**

The final structure is one converged local-minimum candidate obtained from one anneal–quench–minimize trajectory. It should not, by itself, be interpreted as a statistically representative ensemble of metastable grain-boundary states.

## Physical model

| Parameter | Value |
|---|---:|
| Material | BCC Li |
| Grain-boundary type | Symmetric tilt boundary |
| CSL designation | Σ5 |
| Rotation axis | [001] |
| Boundary plane | (310) |
| Misorientation | 36.86989765° |
| Initial lattice constant | 3.506186475312 Å |
| Boundary conditions | periodic in x, y and z |
| Simulation cell | orthogonal |
| Final atom count | 3600 |
| Atoms per grain | 1800 |
| Atom style | `atomic` |
| Li atomic mass | 6.941 |
| Interatomic model | `SNAP_ZuoChenLi_2019_Li__MO_732106099012_000` |
| Production parallelism | 12 MPI tasks × 1 OpenMP thread |

The cell contains two crystallographically equivalent grain boundaries because of periodicity along the GB-normal direction.

### Bicrystal construction

The generator builds grain 1 and obtains grain 2 by mirror symmetry across the GB plane. Coincident sites are detected independently at the two periodic GB planes. Exactly one atom is removed from each coincident pair with balanced removal between the two grains.

| Quantity | Value |
|---|---:|
| Sites per grain before removal | 1860 |
| Total atoms before removal | 3720 |
| Coincident pairs at GB1 | 60 |
| Coincident pairs at GB2 | 60 |
| Removed from grain 1 | 60 |
| Removed from grain 2 | 60 |
| Final atoms | 3600 |

The initial verified box is approximately:

```text
66.5252 × 33.2626 × 35.0619 Å
```

Its verified minimum atomic separation is 2.217507 Å.

## Interatomic potential and bulk reference

All production calculations use the same OpenKIM SNAP model:

```lammps
pair_style kim SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
pair_coeff * * Li
```

A periodic 128-atom BCC bulk reference was relaxed at zero pressure using the same model.

| Bulk-reference result | Value |
|---|---:|
| Relaxed lattice constant | 3.506186475312 Å |
| Potential energy per atom | -1.89998153759969 eV/atom |

The bulk-relaxed lattice constant is used for the final bicrystal construction so that the in-plane CSL periods do not inherit avoidable bulk strain.

## Simulation workflow

| Stage | Operation | Main settings | Main output |
|---|---|---|---|
| 0 | Potential/structure check | `run 0` | successful force and neighbor build |
| 1 | Initial minimization | CG; fixed y/z; relaxed x; chunked checkpoints | `stage1_minimized.data` |
| 2 | Annealing | NVT, 300 K, 20 ps, 20,000 steps | `stage2_annealed.data` |
| 3 | Quenching | NVT, 300 → 1 K, 50 ps, 50,000 steps | `stage3_quenched.data` |
| 4 | Final minimization | 0 K CG; fixed y/z; relaxed x; tight convergence | `stage4_final.data` |

The timestep for Stages 2 and 3 is 0.001 ps.

### Cell-relaxation policy

The grain-boundary plane is the yz plane and x is the GB-normal direction.

- **y and z are fixed** at the zero-pressure bulk/CSL periods.
- **x is relaxed** toward zero normal pressure.

This keeps the intended in-plane crystallographic periodicity while allowing the bicrystal to accommodate GB excess volume. The reported final convergence criterion applies to `Pxx`; nonzero `Pyy` and `Pzz` are not independently minimized because y and z are deliberately constrained.

### Checkpoint and recovery policy

- Stage 1 and Stage 4 run minimization in restartable 500-iteration chunks.
- Stage 2 and Stage 3 write restart files every 1000 steps (1 ps).
- Stage 2 and Stage 3 have dedicated resume inputs that preserve the original target timestep and temperature schedule.
- The generic legacy `in.resume` script was retired because one generic branch could not safely preserve every stage's distinct physical protocol.

## Final converged structure

Authoritative file:

```text
outputs/stage4_final.data
```

| Final result | Value |
|---|---:|
| Atoms | 3600 |
| Atom types | 1 |
| Lx | 66.01192624065516 Å |
| Ly | 33.2626054898 Å |
| Lz | 35.0618647531 Å |
| Maximum atomic force | 3.55281121421774 × 10⁻⁷ eV/Å |
| Force norm | 1.95219377952745 × 10⁻⁶ |
| Pxx | 3.63132085694783 × 10⁻⁴ bar |
| Maximum absolute velocity | 0 |
| Dangerous builds | 0 |

The final data file contains finite coordinates, exactly zero velocities, 3600 atoms and one atom type.

## Repository layout

```text
v2_final/
├── generator/                 Python bicrystal generator and verification code
├── tests/                     Automated regression tests
├── lammps/                    Active verified LAMMPS scripts
│   ├── in.0_check
│   ├── in.1_minimize
│   ├── in.2_anneal
│   ├── in.2_anneal_resume
│   ├── in.3_quench
│   ├── in.3_quench_resume
│   ├── in.4_final
│   ├── in.bulk_energy
│   ├── in.dump2data
│   ├── run_all.sh
│   └── resume.sh
├── outputs/
│   ├── Sigma5_310.*           Verified initial bicrystal exports
│   ├── bulk/                  Bulk reference structures and restart
│   ├── stage1_minimized.*
│   ├── stage2_annealed.*
│   ├── stage2_anneal.dump
│   ├── stage3_quenched.*
│   ├── stage3_quench.dump
│   └── stage4_final.*         Authoritative final structure
├── logs/                      Successful production logs and completion markers
├── validation/                Dump-to-data reconstruction evidence
├── FINAL_STATUS.md            Compact validation summary
├── ENVIRONMENT.txt            Software/environment record
├── SHA256SUMS.txt             File-integrity hashes
├── requirements.txt
└── pytest.ini
```

The final package intentionally excludes old a = 3.51 Å structures, interrupted runs, obsolete scripts, syntax-test outputs, completed checkpoint files, PID files, Python caches and `.bak` files.

## Verify the included package

Activate the Python environment and run:

```bash
python -m pytest tests -q

python -m generator.verify \
  outputs/Sigma5_310.data \
  --map outputs/Sigma5_310.grain_map

python -m generator.reimport_check \
  --dir outputs \
  --prefix Sigma5_310

sha256sum -c SHA256SUMS.txt
```

Expected high-level results:

```text
all tests passed
[verify] PASSED: 3600 atoms (g1=1800, g2=1800)
[reimport] ALL OK
all SHA-256 checks passed
```

## Reproduce the workflow in new directories

Do not reuse the included completed `outputs/` and `logs/` directories for a fresh run. Generate a new input set first:

```bash
python -m generator.generate \
  --a 3.506186475312 \
  --outdir outputs_new

python -m generator.bulk \
  --a 3.506186475312 \
  --outdir outputs_new/bulk
```

Run the bulk reference if required:

```bash
OMP_NUM_THREADS=1 mpirun -np 12 lmp \
  -in lammps/in.bulk_energy \
  -var outdir "$PWD/outputs_new" \
  -var kim_model SNAP_ZuoChenLi_2019_Li__MO_732106099012_000 \
  -log logs_new/bulk_energy.lammps
```

Run Stages 0–4:

```bash
OUTDIR="$PWD/outputs_new" \
LOGDIR="$PWD/logs_new" \
MPI_NP=12 \
OMP_NUM_THREADS=1 \
bash lammps/run_all.sh
```

The supplied runner creates one log per stage and writes a `.done` marker only after successful stage completion.

## Resume an interrupted stage

Use the stage-specific launcher with the same `OUTDIR` and `LOGDIR`:

```bash
OUTDIR="$PWD/outputs_new" LOGDIR="$PWD/logs_new" \
  MPI_NP=12 bash lammps/resume.sh stage1

OUTDIR="$PWD/outputs_new" LOGDIR="$PWD/logs_new" \
  MPI_NP=12 bash lammps/resume.sh stage2

OUTDIR="$PWD/outputs_new" LOGDIR="$PWD/logs_new" \
  MPI_NP=12 bash lammps/resume.sh stage3

OUTDIR="$PWD/outputs_new" LOGDIR="$PWD/logs_new" \
  MPI_NP=12 bash lammps/resume.sh stage4
```

Resume inputs are selected as follows:

| Stage | Resume source |
|---|---|
| 1 | `stage1_checkpoint.data` |
| 2 | latest `stage2_anneal.checkpoint.*` |
| 3 | latest `stage3_quench.checkpoint.*` |
| 4 | `stage4_checkpoint.restart` |

## Reconstruct a data file from a dump

`in.dump2data` uses an existing data file as the template for the box, masses and atom types, then replaces coordinates and velocities from a native LAMMPS dump frame:

```bash
lmp -in lammps/in.dump2data \
  -var template outputs/stage1_minimized.data \
  -var dump outputs/stage2_anneal.dump \
  -var timestep 20000 \
  -var out rebuilt_stage2.data
```

This path was validated by reconstructing timestep 20000 and comparing its atom count, atom types and box dimensions with `stage2_annealed.data`.

## Scientific scope and limitations

The package establishes a reproducible route to one relaxed GB state, but several conclusions require additional work:

1. **Metastable-state sampling:** one thermal trajectory does not map the population of accessible local minima. Independent random perturbations, thermal histories or activation-relaxation searches are required.
2. **Distinct-basin verification:** structures should be compared using energy, local atomic descriptors such as SOAP and a structural-distance/deduplication criterion.
3. **Grain-boundary energy:** the matching bulk reference is included, but a final GB-energy analysis and uncertainty treatment are not yet reported here.
4. **Potential transferability:** the SNAP model is internally consistent for this workflow, but its accuracy for Li GB configurations should be assessed against higher-level reference data before strong quantitative claims.
5. **Finite-size effects:** the periodic cell contains two GBs; convergence with respect to GB separation and in-plane area should be checked for quantitative publication results.

## Software environment

The validated production run used:

- WSL/Linux
- Python 3.14.4
- pytest 9.1.1
- LAMMPS 10 Dec 2025
- OpenKIM support
- OpenMPI, 12 MPI ranks
- one OpenMP thread per MPI rank

See `ENVIRONMENT.txt` and the LAMMPS logs for the exact recorded environment and runtime citations.

## Citation

When using these results, cite the OpenKIM framework and the selected model as requested in the LAMMPS logs:

- OpenKIM Project: <https://doi.org/10.1007/s11837-011-0102-6>
- Model citation page: <https://openkim.org/cite/MO_732106099012_000#item-citation>

## Integrity and provenance

`SHA256SUMS.txt` records the integrity of every included file. The successful Stage 0–4 logs are retained in `logs/`. The earlier working directory and obsolete Windows copy are not part of the authoritative scientific package.

