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
