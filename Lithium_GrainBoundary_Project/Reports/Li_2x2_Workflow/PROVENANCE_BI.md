# Provenance (Bilingual / 双语版)

# 中文版

# 溯源（中文版）

## 旧项目

```text
权威旧项目：
D:\Research\Lithium_GrainBoundary_Project\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production

旧生成器版本：2.0.0
旧晶格常数：  3.506186475312 A
旧超胞：      K=6, J=3, L=10, 3600 atoms
旧盒子：      66.525211 x 33.262605 x 35.061865 A
势函数：      SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
```

旧项目本身是权威参考源，只读使用，未做任何修改。

## 新实验

```text
新超胞：K=6, J=6, L=20
实际原子数：14400
盒子：66.52521097958189 x 66.52521097958189 x 70.12372950624 A
```

新生成器源码与旧版完全一致（`generator/` 是 v2.0.0 的副本），没有修改
任何源码；2×2 通过 CLI 参数 `--k 6 --j 6 --l 20` 实现。所谓 "new
generator changes" 为空；参数层面的改动全部记录在 `PARAMETERS.json`。

## 坐标约定差异记录

导师口头说“模型可在 x 与 y 方向上加长，至少 2 倍”。旧版权威约定为
x=GB normal、yz=GB plane。若机械地把“x、y”当成盒子方向加倍，会把 GB
normal 厚度也翻倍，导致 4×4 面内、2× GB 间距的模型，与任务书“只加倍
GB 面内两个方向”不符。

本实验按权威旧版约定执行：`J=3->6`（y，in-plane）、`L=10->20`（z，
in-plane）、`K=6` 不变（x，GB normal）。

## 环境

```text
生成主机：Windows 11（本机）
Python：  3.12.13（Codex bundled runtime，numpy 2.3.5）
生成器：  v2.0.0（未修改）

MD 主机： 本机不可用
需要：    LAMMPS（已验证旧环境使用 10 Dec 2025 构建）
          OpenKIM 支持
          OpenMPI，推荐 12 ranks
```

旧版生产环境参考：WSL/Linux、Python 3.14.4、pytest 9.1.1、LAMMPS 10 Dec
2025、OpenKIM、OpenMPI 12 ranks。

## 时间线

```text
2026-08-30   Stage 0 完成；旧文件复制完成；协议文档完成；
             MD 脚本已备好；本机无 LAMMPS。
```


# English Version

# Provenance (English)

## Legacy source

```text
authoritative old project:
D:\Research\Lithium_GrainBoundary_Project\01_GB_Construction\Sigma5_310_Metastable_GB\v2_production

old generator version: 2.0.0
old lattice constant:  3.506186475312 A
old supercell:         K=6, J=3, L=10, 3600 atoms
old box:               66.525211 x 33.262605 x 35.061865 A
potential:             SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
```

The old project is the authoritative read-only reference source and was not
modified.

## New experiment

```text
new supercell: K=6, J=6, L=20
actual atoms:  14400
box:           66.52521097958189 x 66.52521097958189 x 70.12372950624 A
```

The new generator source is identical to the legacy version (`generator/` is
a v2.0.0 copy); no source code was modified. The 2x2 extension is achieved
with CLI arguments `--k 6 --j 6 --l 20`. There are therefore no "new
generator changes"; all parameter-level changes are recorded in
`PARAMETERS.json`.

## Coordinate-convention note

The advisor said the model could be elongated in the x and y directions by at
least 2x. The authoritative legacy convention is x = GB normal and yz = GB
plane. Mechanically doubling box directions named "x/y" would also double the
GB-normal thickness, producing a 4x4 in-plane model with twice the GB spacing
(the teammate's K=12/J=6/L=10 scheme), which conflicts with the task of
doubling only the two in-plane GB directions.

This experiment follows the authoritative legacy convention: `J=3 to 6` (y,
in-plane), `L=10 to 20` (z, in-plane), and `K=6` unchanged (x, GB normal).

## Environment

```text
generation host: Windows 11 (this machine)
Python:          3.12.13 (bundled Codex runtime, numpy 2.3.5)
generator:       v2.0.0 (unchanged)

MD host:         not available on this machine
required:        LAMMPS (validated old environment used 10 Dec 2025 build)
                 OpenKIM support
                 OpenMPI, 12 ranks recommended
```

Reference legacy environment: WSL/Linux, Python 3.14.4, pytest 9.1.1, LAMMPS
10 Dec 2025, OpenKIM, OpenMPI 12 ranks.

## Timeline

```text
2026-08-30   Stage 0 completed; legacy files copied; protocol docs written;
             MD scripts prepared; LAMMPS not available on this machine.
```

