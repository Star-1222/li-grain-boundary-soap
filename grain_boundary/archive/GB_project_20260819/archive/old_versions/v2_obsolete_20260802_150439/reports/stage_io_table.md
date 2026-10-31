# Stage input/output table (v2 LAMMPS workflow)

All stages use `boundary p p p`, `units metal`, `atom_style atomic`, and the
same OpenKIM SNAP potential `SNAP_ZuoChenLi_2019_Li__MO_732106099012_000`.

## Cell relaxation policy

- **x (GB normal): fixed.** The two GBs are separated by Lx/2; relaxing x
  would arbitrarily change the GB spacing, so x is never relaxed.
- **y, z (in-plane): relaxed at P = 0** in stage 1 and stage 4. The in-plane
  CSL periods sqrt(10)*a and a are lattice degrees of freedom; relaxing them
  removes artificial in-plane stress.
- Full pressure tensor (pxx pyy pzz pxy pxz pyz) is printed in every thermo block.

| stage | input | output | log |
|---|---|---|---|
| 0 check | `Sigma5_310.data` | run 0 sanity (forces OK) | `logs/stage0_check.lammps` |
| 1 minimize | `Sigma5_310.data` | `stage1_minimized.data`, `stage1_minimized.restart` | `logs/stage1_minimize.lammps` |
| 2 anneal 300 K 20 ps | `stage1_minimized.restart` | `stage2_annealed.data`, `stage2_annealed.restart`, `stage2_anneal.dump`, periodic `stage2_anneal.checkpoint` | `logs/stage2_anneal.lammps` |
| 3 quench 300→1 K 50 ps | `stage2_annealed.restart` | `stage3_quenched.data`, `stage3_quenched.restart`, `stage3_quench.dump`, periodic `stage3_quench.checkpoint` | `logs/stage3_quench.lammps` |
| 4 final 0 K minimize | `stage3_quenched.restart` | `stage4_final.data`, `stage4_final.restart`, `stage4_final.dump` | `logs/stage4_final.lammps` |
| resume | any restart/checkpoint | `resume_<stage>.data`, `resume_<stage>.restart`, `resume_<stage>.dump` | user-chosen `-log` |
| dump2data | lammpstrj dump + frame | rebuilt data file | user-chosen `-log` |
| bulk energy | `bulk_Li.data` | `bulk_Li_minimized.data`, `.restart` | user-chosen `-log` |

## Safety / resume

- `run_all.sh` writes `logs/<stage>.done` only after a stage finishes; re-running
  skips finished stages and continues from the first unfinished one.
- stages 2 and 3 write periodic checkpoints (`restart 10000 <file>.checkpoint`);
  `in.resume` can continue from any checkpoint with a chosen number of extra steps.
- every stage uses its own log file; anneal/quench/final outputs never overwrite.
