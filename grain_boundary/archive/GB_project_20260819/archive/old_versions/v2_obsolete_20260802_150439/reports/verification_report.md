# Verification report — v2 Sigma5 (310)[001] BCC Li bicrystal

Reference values from the independent Node implementation
(`work/cs5_310_ref_result.json`); the Python pipeline (`python -m
generator.generate` + `python -m generator.verify`) re-runs every check.

## 1. Misorientation and rotation axis

| quantity | value |
|---|---|
| misorientation angle | 36.8698976458 deg (target 2*atan(1/3)) |
| rotation axis | [0, 0, -1] (∥ [001]) |
| cubic-symmetry minimum angle | 36.8698976458 deg (24 rotations scanned) |

## 2. GB plane

| frame | normal (numeric) | Miller |
|---|---|---|
| crystal 1 | (0.9486832981, -0.3162277660, 0) | [3,-1,0] |
| crystal 2 | (0.9486832981, 0.3162277660, 0) | [3,1,0] |

## 3. Two GBs: positions and equivalence

- GB1 at x = 0 (periodic), GB2 at x = Lx/2 = 33.2988 A
- cross-GB pairs < 4.5 A: **900** at GB1, **900** at GB2
- distance spectra **identical** (pair-by-pair, atol 1e-9)
- min cross-GB distance: 2.219919 A at both GBs; 60 pairs at 2.2199 A per GB

## 4. Bulk BCC away from the GBs

- margin = 5 A; sampled = **2520** atoms, perfect BCC (8 NN + 6 NNN) = **2520**, failed = 0

## 5. Distances

- global min = 2.219918917438207 A = analytic 2d = 2*3.51/sqrt(10) = 2.219918917438202 A
- BCC NN = 3.039749167283379 A, NNN = 3.51 A
- neighbour histogram (A): {2.2: 120, 3.0: 14160, 3.5: 10560, 4.2: 240, 4.4: 600}

## 6. Voids

- 1-A bins along x: 67 bins, **no empty bin**, max empty gap = 0.5 A
- density profile (1-A bins) starts: [30, 60, 60, 60, 60, ...] atoms/bin

## 7. Removal log

- total removed = 120 (g1 = 60, g2 = 60; GB1 = 60, GB2 = 60)
- all deleted atoms lie on a GB plane; all are cross-grain coincidences
- full record: `outputs/removal_log.csv`

## 8. Counts and density

- atoms: 3600 (g1 = 1800, g2 = 1800)
- density = 0.0462496708468738 A^-3 = 2/a^3 (exact)

## 9. In-plane periodicity

- lattice-index misses: grain1 = 0, grain2 = 0
- y period = sqrt(10)*a = 11.099594587191012 A; z period = a = 3.51 A

## 10. Legacy regression

- this model: 3720 → 3600, removal 3.2258 %
- legacy model: 3600 → 3224, removal 10.44 % — **not reproduced**

## 11. Old vs new model

| quantity | legacy | v2 |
|---|---|---|
| atoms before removal | 3600 (claimed) | 3720 |
| atoms after removal | 3224 | 3600 |
| removal fraction | 10.44 % | 3.2258 % |
| box | 85.977 x 60.795 x 14.892 A triclinic (xy=-17.195, xz=-8.598, yz=-6.079) | 66.5976 x 33.2988 x 35.1 A orthogonal |
| GBs | unclear (tilt factors) | two explicit equivalent GBs at x=0 and x=Lx/2 |
| boundary | ASE pbc (F,T,T) vs LAMMPS p p p | p p p everywhere |

**Overall: PASSED (all 11 checks).**
