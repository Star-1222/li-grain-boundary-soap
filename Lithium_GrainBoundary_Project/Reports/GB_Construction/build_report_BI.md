# Bilingual Report / 双语报告

本文件为双语版：中文版在前，英文原版在后。This file contains the Chinese version first, followed by the full English original.

# 中文版

# 构建报告：v2 Sigma5 (310)[001] BCC Li 双晶（Sigma5_310）（中文版）

## 构建方法

修正后的生成器使用**镜像（CSL）构建**：

1. 晶粒 1 由其自身旋转的 BCC 晶格生成，取向矩阵 M1 = Rz(+theta/2)，
   theta = 2*atan(1/3) = 36.86989765 deg。使用常规立方 BCC 基矢
   (0,0,0)+(1/2,1/2,1/2)；不再重复旧版生成器的原胞错误。
2. 晶粒 2 是晶粒 1 关于 x = Lx/2 的严格镜像。绝不在固定周期盒子内做额外
   旋转。
3. 重叠移除只删除位于两个晶界平面（x=0 和 x=Lx/2）上的跨晶粒重合对，
   不删除其他原子。
4. 移除：120 个原子（3.2258%），g1=60、g2=60，GB1=60、GB2=60。

## 参数

- 晶格常数 a = 3.506186475312 A
- 取向差 theta = 36.86989765 deg，绕 [001]（Sigma5）
- 晶界平面 (310)；盒子周期 K=6、J=3、L=10
- 盒子：66.525211 x 33.262605 x 35.061865 A（正交，p p p）
- 面内周期：y = sqrt(10)*a = 11.087535 A，z = a = 3.506186475312 A
- 原子数：3720 -> 3600（每晶粒 1800/1800）
- 密度：0.046400746239 atoms/A^3 = 2/a^3

## 取向矩阵（行为实验坐标系中的晶轴）

M1（晶粒 1）= Rz(+theta/2)：

    [[0.9486832981, -0.3162277660, 0],
     [0.3162277660,  0.9486832981, 0],
     [0, 0, 1]]

M2（晶粒 2）= Rz(-theta/2)：

    [[0.9486832981,  0.3162277660, 0],
     [-0.3162277660, 0.9486832981, 0],
     [0, 0, 1]]

M2 . M1^T = Rz(-theta)，绕 [001] 角度 36.86989765 deg。

## 输出文件

- `outputs/Sigma5_310.data`       LAMMPS data（正交，Velocities=0）
- `outputs/Sigma5_310.xyz`        含 grain_id 的扩展 XYZ
- `outputs/Sigma5_310.cif`        CIF（P1）
- `outputs/Sigma5_310_grain_id.xyz`  含 grain_id 列的普通 XYZ
- `outputs/Sigma5_310.grain_map`  atom_id -> grain_id 映射
- `outputs/removal_log.csv`       被删原子的完整记录
- `outputs/params.json`           构建参数
- `outputs/verification_report.md`  全部几何检查

## 下一步（LAMMPS）

运行 v2/lammps/ 中的分阶段工作流（见 v2/README.md）。
势函数：OpenKIM SNAP_ZuoChenLi_2019_Li__MO_732106099012_000（暂定）。


# English Original

# Build report: v2 Sigma5 (310)[001] BCC Li bicrystal (Sigma5_310)

## Construction method

The corrected generator uses the **mirror (CSL) construction**:

1. Grain 1 is generated from its own rotated BCC lattice with orientation
   matrix M1 = Rz(+theta/2), theta = 2*atan(1/3) = 36.86989765 deg.
   Conventional cubic BCC basis (0,0,0)+(1/2,1/2,1/2) is used; the
   primitive-cell mistake of the legacy generator is not repeated.
2. Grain 2 is the exact mirror image of grain 1 across x = Lx/2.
   No rotation inside a fixed periodic box is ever performed.
3. Overlap removal deletes ONLY cross-grain coincident pairs sitting on
   the two GB planes (x=0 and x=Lx/2). Nothing else is deleted.
4. Removal: 120 atoms (3.2258%),
   g1=60, g2=60,
   GB1=60,
   GB2=60.

## Parameters

- lattice constant a = 3.506186475312 A
- misorientation theta = 36.86989765 deg about [001] (Sigma5)
- GB plane (310); box periods K=6, J=3, L=10
- box: 66.525211 x 33.262605 x 35.061865 A (orthogonal, p p p)
- in-plane periods: y = sqrt(10)*a = 11.087535 A, z = a = 3.506186475312 A
- atoms: 3720 -> 3600 (per grain 1800/1800)
- density: 0.046400746239 atoms/A^3 = 2/a^3

## Orientation matrices (rows = crystal axes in the lab frame)

M1 (grain 1) = Rz(+theta/2):
    [[0.9486832981, -0.3162277660, 0],
     [0.3162277660,  0.9486832981, 0],
     [0, 0, 1]]
M2 (grain 2) = Rz(-theta/2):
    [[0.9486832981,  0.3162277660, 0],
     [-0.3162277660, 0.9486832981, 0],
     [0, 0, 1]]
M2 . M1^T = Rz(-theta), angle 36.86989765 deg about [001].

## Files written

- `outputs/Sigma5_310.data`       LAMMPS data (orthogonal, Velocities=0)
- `outputs/Sigma5_310.xyz`        extended XYZ with grain_id
- `outputs/Sigma5_310.cif`        CIF (P1)
- `outputs/Sigma5_310_grain_id.xyz`  plain XYZ with grain_id column
- `outputs/Sigma5_310.grain_map`  atom_id -> grain_id map
- `outputs/removal_log.csv`     full record of deleted atoms
- `outputs/params.json`         construction parameters
- `outputs/verification_report.md`  all geometric checks

## Next steps (LAMMPS)

Run the staged workflow in v2/lammps/ (see v2/README.md).
Potential: OpenKIM SNAP_ZuoChenLi_2019_Li__MO_732106099012_000 (tentative).

