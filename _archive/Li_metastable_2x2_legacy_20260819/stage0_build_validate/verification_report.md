# Geometric verification report (v2 Sigma5 (310)[001] BCC Li)

- a = 3.506186475312 A, theta = 36.86989765 deg, axis [001], plane 310
- box = 66.525211 x 66.525211 x 70.123730 A (K=6 J=6 L=20)
- in-plane periods y = 11.087535 A, z = 3.506186 A

## 1. Misorientation and rotation axis
- angle = 36.86989765 deg (target 36.86989765)
- axis = [0.0, 0.0, -1.0] (target [001])
- cubic-symmetry minimum = 36.86989765 deg

## 2. GB plane
- normal in crystal 1 = [0.9486832980505139, -0.316227766016838, 0.0] ~ [3,-1,0]
- normal in crystal 2 = [0.9486832980505139, 0.316227766016838, 0.0] ~ [3,1,0]

## 3. Two GBs, explicit and equivalent
- GB1 at x = 0.0, GB2 at x = 33.26260548979094 A
- cross-GB pairs < 4.5 A: GB1 = 3600, GB2 = 3600
- spectra identical: True
- min cross-GB distance: GB1 = 2.217507 A, GB2 = 2.217507 A
- pairs at 2.2199 A: GB1 = 0, GB2 = 0

## 4. Bulk BCC away from the GBs
- margin = 5.0 A; sampled = 10080, perfect = 10080, failed = 0

## 5. Distances
- global min = 2.217507 A (analytic 2d = 2.217507 A)
- BCC NN = 3.036447 A, NNN = 3.506186 A

## 6. Voids
- occupied x planes: 60; expected gap = 1.1087535163263649 A; gap range = [1.1087535163263595, 1.1087535163263738] A; missing planes = []

## 7. Removal log
- total removed = 480 (g1 = 240, g2 = 240; GB1 = 240, GB2 = 240)
- all deleted atoms on GB planes: True

## 8. Counts and density
- atoms = 14400 (g1 = 7200, g2 = 7200)
- density = 0.046400746239 atoms/A^3 (2/a^3 = 0.046400746239)

## 9. In-plane periodicity (y, z)
- lattice-index misses: grain1 = 0, grain2 = 0
- y period = 11.087535 A = sqrt(10)*a; z period = 3.506186 A = a

## 10. Legacy regression
- this model: 14880 -> 14400 atoms, removal 3.2258%
- legacy model: 3600 -> 3224, 10.44% (must NOT be reproduced)
- legacy-bad scenario present: False

**Overall: PASSED = True**
