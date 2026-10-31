# Build report — v2 Sigma5 (310)[001] BCC Li symmetric tilt grain boundary

Status: **files generated; geometry verified against the independent Node reference
results** (`work/cs5_310_ref_result.json`).  Python/ASE/LAMMPS runs are pending
(execute the WSL command sequence in `v2/README.md`).

## 1. Why the legacy generator was wrong

| # | Legacy bug | Consequence |
|---|---|---|
| 1 | `ase.build.bulk("Li", "bcc")` returns the **primitive** BCC cell, but the code treated it as the conventional cubic cell | wrong Miller-index arithmetic, wrong lattice geometry |
| 2 | the two half-crystals were **rotated inside an already-built periodic box** | box vectors incompatible with the rotated grains; periodicity violated |
| 3 | overlap removal deleted **every atom closer than a cutoff anywhere in the system** | 3600 → 3224 atoms, 10.44 % removed (376 atoms), including interior atoms |
| 4 | ASE cell used `pbc=(False, True, True)` while LAMMPS ran with `boundary p p p` | inconsistent periodicity between the generator and the MD engine |
| 5 | source code said `SUPERCELL=(20,12,12)` but the actual output used `(10,6,6)` | untraceable model |
| 6 | `in.resume` was identical to the main input; `dump2data.in` used an empty box with default `add no` | resume and dump→data workflows broken |

## 2. Corrected crystallographic construction (mirror / CSL method)

- BCC Li, **conventional cubic cell** with basis `(0,0,0)` and `(1/2,1/2,1/2)`.
- Lattice constant `a = 3.51 A` (same as legacy; consistent with the ZuoChenLi SNAP family).
- Misorientation `theta = 2*atan(1/3) = 36.86989764584402 deg` about `[001]`, Sigma 5.
- GB plane `(310)`; the GB normal in crystal-1 frame is `[3,-1,0]`, in crystal-2 frame `[3,1,0]`.
- Orientation matrices (rows = crystal axes in lab frame):

```
M1 = Rz(+theta/2) =
  [ 0.9486832981  -0.3162277660  0 ]
  [ 0.3162277660   0.9486832981  0 ]
  [ 0             0             1 ]

M2 = Rz(-theta/2) =
  [ 0.9486832981   0.3162277660  0 ]
  [-0.3162277660   0.9486832981  0 ]
  [ 0              0             1 ]

M2 . M1^T = Rz(-theta):  36.86989765 deg about [001]
```

- Box (orthogonal, periodic in all three directions):

```
Lx = 6*sqrt(10)*a = 66.59756752314607 A      (GB normal, 2 GBs at x=0 and x=Lx/2)
Ly = 3*sqrt(10)*a = 33.298783761573034 A     (in-plane, J=3 CSL periods)
Lz = 10*a         = 35.1 A                   (in-plane, L=10 periods)
```

- In-plane periods: `y = sqrt(10)*a = 11.099594587191012 A`, `z = a = 3.51 A`.

### Construction steps

1. **Grain 1** is generated from its own rotated lattice, one atom at a time
   (no rotation of any pre-built box):
   ```
   x = d*(3p + q + 2b),   y = d*(-p + 3q + b),   z = a*(r + b/2)
   p, q, r in Z, b in {0,1}, d = a/sqrt(10)
   ```
   domain: `0 <= x <= Lx/2`, `0 <= y < Ly`, `0 <= z < Lz`.
2. **Grain 2** is the exact mirror image of grain 1 across `x = Lx/2`:
   `(x, y, z) -> (Lx - x, y, z)`, with `x = Lx` wrapped to `x = 0`.
   This guarantees the two GBs are strictly equivalent and the in-plane
   periods match exactly.  No global deletion is ever performed.
3. **Overlap removal** only considers cross-grain coincident pairs sitting
   on the GB planes (`|x| < 0.1 A` or `|x - Lx/2| < 0.1 A`, pair distance
   `< 0.05 A`).  One atom of each pair is removed, alternating the grain so
   the mirror symmetry is preserved.  Every deletion is logged to
   `outputs/removal_log.csv` (atom_id, grain_id, x, y, z, nearest_distance,
   GB, partner_atom_id).

### Counts

- sites per grain: 1860 (31 (310) planes x 60 atoms/plane) → 3720 total
- coincident pairs: 60 on GB1 + 60 on GB2 = 120
- removed: 120 (g1=60, g2=60; GB1=60, GB2=60) — 3.2258 % (not the legacy 10.44 %)
- final: **3600 atoms** (g1=1800, g2=1800)
- density: 3600 / V = **0.0462496708468738 A^-3 = 2/a^3** (exact)

## 3. Exported files (written by `python -m generator.generate`)

- `Sigma5_310.data` — LAMMPS data (orthogonal box, Velocities = 0)
- `Sigma5_310.xyz` — extended XYZ with grain_id per atom
- `Sigma5_310.cif` — CIF (P1)
- `Sigma5_310_grain_id.xyz` — plain XYZ with grain_id column
- `Sigma5_310.grain_map` — atom_id → grain_id map
- `removal_log.csv` — full deletion record (120 rows)
- `params.json` — construction parameters
- `verification_report.md` — all geometric checks

## 4. Report provenance

The numbers in this report and in `verification_report.md` were produced by an
independent implementation (Node/JavaScript) that mirrors the closed-form
geometry of this package, and are stored in `work/cs5_310_ref_result.json`.
`python -m generator.generate` re-runs the same checks in Python and writes the
reports into the output directory.
