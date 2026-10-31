# v2 — corrected Sigma5 (310)[001] BCC Li symmetric tilt grain boundary

Re-implementation of the bicrystal generator and the full LAMMPS workflow.
Legacy results are untouched in `../legacy_results/`; all new results go to
`./outputs/`.

## Layout

```
v2/
  generator/       builder, geometry core, exporters, verification, CLI
  lammps/          staged LAMMPS inputs + run/resume shell scripts
  tests/           pytest suite (no external tools needed except numpy)
  reports/         build/verification reports, stage I/O table
  outputs/         generated structures + reports (empty until you run)
  logs/            LAMMPS logs + checkpoint markers (empty until you run)
```

## Dependencies (WSL Ubuntu)

```bash
pip install --upgrade numpy ase pytest

# OpenKIM model for the SNAP potential
kim-api-collections-management install SNAP_ZuoChenLi_2019_Li__MO_732106099012_000

# LAMMPS must be built with the KIM package
which lmp && lmp -h | grep -i KIM
```

## Quick start (in WSL, from the v2 directory)

```bash
cd /mnt/c/Users/1/Documents/Codex/2026-08-01/bcc-li-5-001-sigma5-310/v2

# 1) tests
python -m pytest tests -v

# 2) small model + geometry verification
python -m generator.generate --model-size small --outdir outputs/small
python -m generator.verify outputs/small/Sigma5_310_small.data \
    --map outputs/small/Sigma5_310_small.grain_map

# 3) full model + geometry verification
python -m generator.generate --outdir outputs
python -m generator.verify outputs/Sigma5_310.data \
    --map outputs/Sigma5_310.grain_map

# 4) re-import the exported files with ASE
python -m generator.reimport_check --dir outputs --prefix Sigma5_310

# 5) LAMMPS stage 0 (run 0: load potential + forces)
cd lammps
lmp -in in.0_check -var outdir ../outputs \
    -var kim_model SNAP_ZuoChenLi_2019_Li__MO_732106099012_000 \
    -log ../logs/stage0_check.lammps

# 6) short minimization (stage 1 only; full anneal/quench later)
lmp -in in.1_minimize -var outdir ../outputs \
    -var kim_model SNAP_ZuoChenLi_2019_Li__MO_732106099012_000 \
    -log ../logs/stage1_minimize.lammps
```

## Full staged workflow (after geometry validation)

```bash
cd /mnt/c/Users/1/Documents/Codex/2026-08-01/bcc-li-5-001-sigma5-310/v2/lammps
bash run_all.sh                 # stage0 -> stage4, checkpoint-safe
bash resume.sh                  # re-run: skips finished stages
```

## Mid-stage resume example (quench interrupted at a checkpoint)

```bash
lmp -in in.resume \
    -var restart ../outputs/stage3_quench.checkpoint -var stage quench \
    -var nsteps 30000 -var tstart 300.0 -var tstop 1.0 \
    -var outdir ../outputs \
    -var kim_model SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
```

## Forbidden / mandated workflow rules

- NVT target temperature is never 0 K (quench goes to 1 K, then velocities
  are zeroed and a 0 K minimization is performed).
- final minimization uses `minimize 0.0 1.0e-6 10000 100000` (etol = 0 so the
  minimizer cannot stop early on the energy criterion).
- x (GB normal) is fixed; y and z are relaxed at P = 0; the full pressure
  tensor is printed at every stage.
- the final data file velocities are all zero (checked by tests).

## Final review fixes (2026-08-01)

During the final static review of this v2 package, two bugs were found and
fixed before any run was attempted:

1. `generator/geometry.cubic_rotation_group()` used only 4 sign patterns, so
   it returned the 12-element tetrahedral subgroup instead of the full
   24-element proper cubic (octahedral) group.  It now scans all 8 sign
   patterns and keeps the 24 matrices with det > 0.  This affects the
   cubic-symmetry misorientation scan (verified to be exactly 36.86989765 deg
   for Rz(-36.87 deg) under the full group).
2. `generator/lammps_io.read_lammps_data()` matched the `Atoms # atomic`
   section header with the atom-count pattern and would raise
   `int("Atoms")`.  The header pattern is now
   `^d+ atoms# v2 — corrected Sigma5 (310)[001] BCC Li symmetric tilt grain boundary

Re-implementation of the bicrystal generator and the full LAMMPS workflow.
Legacy results are untouched in `../legacy_results/`; all new results go to
`./outputs/`.

## Layout

```
v2/
  generator/       builder, geometry core, exporters, verification, CLI
  lammps/          staged LAMMPS inputs + run/resume shell scripts
  tests/           pytest suite (no external tools needed except numpy)
  reports/         build/verification reports, stage I/O table
  outputs/         generated structures + reports (empty until you run)
  logs/            LAMMPS logs + checkpoint markers (empty until you run)
```

## Dependencies (WSL Ubuntu)

```bash
pip install --upgrade numpy ase pytest

# OpenKIM model for the SNAP potential
kim-api-collections-management install SNAP_ZuoChenLi_2019_Li__MO_732106099012_000

# LAMMPS must be built with the KIM package
which lmp && lmp -h | grep -i KIM
```

## Quick start (in WSL, from the v2 directory)

```bash
cd /mnt/c/Users/1/Documents/Codex/2026-08-01/bcc-li-5-001-sigma5-310/v2

# 1) tests
python -m pytest tests -v

# 2) small model + geometry verification
python -m generator.generate --model-size small --outdir outputs/small
python -m generator.verify outputs/small/Sigma5_310_small.data \
    --map outputs/small/Sigma5_310_small.grain_map

# 3) full model + geometry verification
python -m generator.generate --outdir outputs
python -m generator.verify outputs/Sigma5_310.data \
    --map outputs/Sigma5_310.grain_map

# 4) re-import the exported files with ASE
python -m generator.reimport_check --dir outputs --prefix Sigma5_310

# 5) LAMMPS stage 0 (run 0: load potential + forces)
cd lammps
lmp -in in.0_check -var outdir ../outputs \
    -var kim_model SNAP_ZuoChenLi_2019_Li__MO_732106099012_000 \
    -log ../logs/stage0_check.lammps

# 6) short minimization (stage 1 only; full anneal/quench later)
lmp -in in.1_minimize -var outdir ../outputs \
    -var kim_model SNAP_ZuoChenLi_2019_Li__MO_732106099012_000 \
    -log ../logs/stage1_minimize.lammps
```

## Full staged workflow (after geometry validation)

```bash
cd /mnt/c/Users/1/Documents/Codex/2026-08-01/bcc-li-5-001-sigma5-310/v2/lammps
bash run_all.sh                 # stage0 -> stage4, checkpoint-safe
bash resume.sh                  # re-run: skips finished stages
```

## Mid-stage resume example (quench interrupted at a checkpoint)

```bash
lmp -in in.resume \
    -var restart ../outputs/stage3_quench.checkpoint -var stage quench \
    -var nsteps 30000 -var tstart 300.0 -var tstop 1.0 \
    -var outdir ../outputs \
    -var kim_model SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
```

## Forbidden / mandated workflow rules

- NVT target temperature is never 0 K (quench goes to 1 K, then velocities
  are zeroed and a 0 K minimization is performed).
- final minimization uses `minimize 0.0 1.0e-6 10000 100000` (etol = 0 so the
  minimizer cannot stop early on the energy criterion).
- x (GB normal) is fixed; y and z are relaxed at P = 0; the full pressure
  tensor is printed at every stage.
- the final data file velocities are all zero (checked by tests).

-style (line ends with "atoms", exactly two tokens, first
   token numeric).  Re-import verification (`verify_data_file`,
   `reimport_check`) depends on this parser.

## Reports

- `reports/build_report.md` — construction method, parameters, counts
- `reports/verification_report.md` — all geometric checks with numbers
- `reports/stage_io_table.md` — per-stage inputs/outputs, resume policy
- `reports/reference_energy_bulk.md` — bulk reference energy protocol

The numeric values in the reports were cross-checked with an independent
implementation stored in `../work/cs5_310_ref_result.json`.
