# Bilingual Report / 双语报告

本文件为双语版：中文版在前，英文原版在后。This file contains the Chinese version first, followed by the full English original.

# 中文版

# 几何验证报告：v2 Sigma5 (310)[001] BCC Li（中文版）

- a = 3.506186475312 A，theta = 36.86989765 deg，轴 [001]，平面 310
- 盒子 = 66.525211 x 33.262605 x 35.061865 A（K=6 J=3 L=10）
- 面内周期 y = 11.087535 A，z = 3.506186 A

## 1. 取向差与旋转轴

- 角度 = 36.86989765 deg（目标 36.86989765）
- 轴 = [0.0, 0.0, -1.0]（目标 [001]）
- 立方对称最小值 = 36.86989765 deg

## 2. 晶界平面

- 晶粒 1 中的法向 = [0.9486832980505138, -0.31622776601683794, 0.0]
  ~ [3,-1,0]
- 晶粒 2 中的法向 = [0.9486832980505138, 0.31622776601683794, 0.0]
  ~ [3,1,0]

## 3. 两条明确且等价的 GB

- GB1 在 x = 0.0，GB2 在 x = 33.26260548979094 A
- 跨 GB 距离 < 4.5 A 的配对：GB1 = 900，GB2 = 900
- 谱完全相同：True
- 最小跨 GB 距离：GB1 = 2.217507 A，GB2 = 2.217507 A
- 距离 2.2199 A 处的配对：GB1 = 0，GB2 = 0

## 4. 远离 GB 的体相 BCC

- 边界 = 5.0 A；采样 = 2520，完美 = 2520，失败 = 0

## 5. 距离

- 全局最小 = 2.217507 A（解析 2d = 2.217507 A）
- BCC 最近邻 = 3.036447 A，次近邻 = 3.506186 A

## 6. 空洞

- 占据的 x 平面：60；预期间隙 = 1.1087535163263649 A；间隙范围 =
  [1.1087535163263595, 1.1087535163263738] A；缺失平面 = []

## 7. 移除记录

- 总移除 = 120（g1 = 60，g2 = 60；GB1 = 60，GB2 = 60）
- 所有被删原子都在 GB 平面上：True

## 8. 原子数与密度

- 原子 = 3600（g1 = 1800，g2 = 1800）
- 密度 = 0.046400746239 atoms/A^3（2/a^3 = 0.046400746239）

## 9. 面内周期性（y, z）

- 晶格索引未命中：晶粒1 = 0，晶粒2 = 0
- y 周期 = 11.087535 A = sqrt(10)*a；z 周期 = 3.506186 A = a

## 10. 旧版回归

- 本模型：3720 -> 3600 原子，移除 3.2258%
- 旧模型：3600 -> 3224，10.44%（不得复现）
- 旧版错误场景存在：False

**总体：PASSED = True**


# English Original

# Geometric verification report (v2 Sigma5 (310)[001] BCC Li)

- a = 3.506186475312 A, theta = 36.86989765 deg, axis [001], plane 310
- box = 66.525211 x 33.262605 x 35.061865 A (K=6 J=3 L=10)
- in-plane periods y = 11.087535 A, z = 3.506186 A

## 1. Misorientation and rotation axis
- angle = 36.86989765 deg (target 36.86989765)
- axis = [0.0, 0.0, -1.0] (target [001])
- cubic-symmetry minimum = 36.86989765 deg

## 2. GB plane
- normal in crystal 1 = [0.9486832980505138, -0.31622776601683794, 0.0] ~ [3,-1,0]
- normal in crystal 2 = [0.9486832980505138, 0.31622776601683794, 0.0] ~ [3,1,0]

## 3. Two GBs, explicit and equivalent
- GB1 at x = 0.0, GB2 at x = 33.26260548979094 A
- cross-GB pairs < 4.5 A: GB1 = 900, GB2 = 900
- spectra identical: True
- min cross-GB distance: GB1 = 2.217507 A, GB2 = 2.217507 A
- pairs at 2.2199 A: GB1 = 0, GB2 = 0

## 4. Bulk BCC away from the GBs
- margin = 5.0 A; sampled = 2520, perfect = 2520, failed = 0

## 5. Distances
- global min = 2.217507 A (analytic 2d = 2.217507 A)
- BCC NN = 3.036447 A, NNN = 3.506186 A

## 6. Voids
- occupied x planes: 60; expected gap = 1.1087535163263649 A; gap range = [1.1087535163263595, 1.1087535163263738] A; missing planes = []

## 7. Removal log
- total removed = 120 (g1 = 60, g2 = 60; GB1 = 60, GB2 = 60)
- all deleted atoms on GB planes: True

## 8. Counts and density
- atoms = 3600 (g1 = 1800, g2 = 1800)
- density = 0.046400746239 atoms/A^3 (2/a^3 = 0.046400746239)

## 9. In-plane periodicity (y, z)
- lattice-index misses: grain1 = 0, grain2 = 0
- y period = 11.087535 A = sqrt(10)*a; z period = 3.506186 A = a

## 10. Legacy regression
- this model: 3720 -> 3600 atoms, removal 3.2258%
- legacy model: 3600 -> 3224, 10.44% (must NOT be reproduced)
- legacy-bad scenario present: False

**Overall: PASSED = True**

