# Provenance

## 旧项目 / Legacy source

```text
authoritative old project:
D:\Research\Lithium_GrainBoundary_Project\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production

old generator version: 2.0.0
old lattice constant:  3.506186475312 A
old supercell:         K=6, J=3, L=10, 3600 atoms
old box:               66.525211 x 33.262605 x 35.061865 A
potential:             SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
```

旧项目本身是权威参考源，只读使用，未做任何修改。

## 新实验 / New experiment

```text
new supercell: K=6, J=6, L=20
actual atoms:  14400
box:           66.52521097958189 x 66.52521097958189 x 70.12372950624 A
```

新生成器源码与旧版完全一致（`generator/` 是 v2.0.0 的副本），没有修改
任何源码；2×2 通过 CLI 参数 `--k 6 --j 6 --l 20` 实现。所谓“new
generator changes”因此为空；参数层面的改动全部记录在
`PARAMETERS.json`。

## 坐标约定差异记录

导师口头说“模型可在 x 与 y 方向上加长，至少 2 倍”。旧版权威约定为
x=GB normal、yz=GB plane。若机械地把“x、y”当成盒子方向加倍，会把
GB normal 厚度也翻倍，导致 4×4 面内、2× GB 间距的模型（也就是队友的
K=12/J=6/L=10 方案），这与任务书“只加倍 GB 面内两个方向”不符。

本实验按权威旧版约定执行：`J=3→6`（y，in-plane）、`L=10→20`（z，
in-plane）、`K=6` 不变（x，GB normal）。

## 环境 / Environment

```text
generation host: Windows 11 (this machine)
Python:          3.12.13 (bundled Codex runtime, numpy 2.3.5)
generator:       v2.0.0 (unchanged)

MD host:         not available on this machine
required:        LAMMPS (validated old environment used 10 Dec 2025 build)
                 OpenKIM support
                 OpenMPI, 12 ranks recommended
```

旧版生产环境参考（`legacy/reports_v2/ENVIRONMENT.txt`）：WSL/Linux、
Python 3.14.4、pytest 9.1.1、LAMMPS 10 Dec 2025、OpenKIM、OpenMPI 12 ranks。

## 时间线 / Timeline

```text
2026-08-30   Stage 0 completed; legacy files copied; protocol docs written;
             MD scripts prepared; LAMMPS not available on this machine.
```
