# Li-SNAP + OpenKIM + LAMMPS Environment Verification

## Project Structure

```
project/
├── in.test_li      # LAMMPS input: BCC Li 2x2x2, pair_style kim Li-SNAP
├── run.sh          # Convenience run script
└── README.md       # This file
```

## Prerequisites

### 1. LAMMPS with KIM Package

Make sure LAMMPS was compiled with the KIM package:

```bash
lmp -h | grep -i kim
```

Expected keywords in output:
- KIM
- kim_init
- kim_query
- kim_property
- kim_interactions

### 2. OpenKIM API

Check if the KIM API collections tool is installed:

```bash
which kim-api-collections-management
```

If not found, install on Ubuntu:

```bash
sudo apt update
sudo apt install kim-api
```

### 3. Install Li-SNAP Model

List installed KIM models:

```bash
kim-api-collections-management list
```

Search for `SNAP_ZuoChenLi_2019_Li__MO_732106099012_000`.

If not found, install it:

```bash
kim-api-collections-management install SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
```

## Run the Test

```bash
cd project/
lmp -in in.test_li
```

Or use the run script:

```bash
cd project/
bash run.sh
```

## How to Confirm Li-SNAP is Loaded Correctly

### 1. KIM Model Loading

Look for a message near the top of LAMMPS stdout:

```
KIM: Loading model ...
KIM: Model name: SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
```

### 2. Minimization Convergence

The minimization should converge cleanly. Sample output:

```
Step         Temp          PotEng         Press         Lx            Ly            Lz            Volume
0            ...           ...            ...           ...           ...           ...           ...
...
Minimization stats:
  E_init, E_final, ...
```

### 3. Output Quantities

- **Total Energy (pe)**: For 16 Li atoms (2x2x2 BCC), expect a negative value around -20 to -30 eV
- **Pressure**: Should be near 0 bar after minimization (indicating equilibrium)
- **Lattice constant**: Should converge to approximately 3.5 Angstroms per cell

### 4. No Error Messages

Watch for these errors:
- `ERROR: Unknown pair style kim` -- KIM package not compiled into LAMMPS
- `ERROR: Could not find KIM model` -- model not installed in KIM database
- `ERROR: pair_style kim ...` -- syntax error

## Troubleshooting

### "Unknown pair style kim"

Rebuild LAMMPS with KIM package enabled:

```bash
make yes-KIM
make mpi   # or your build target
```

### "Could not find KIM model"

Install the model:

```bash
kim-api-collections-management install SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
```

### pair_style kim syntax error

The current `in.test_li` uses the LAMMPS 10 Dec 2025 documented syntax:

```
pair_style kim SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
pair_coeff * * Li
```

If this does not work for your LAMMPS version, try the alternative form (model name in pair_coeff instead):

```
pair_style kim
pair_coeff * * SNAP_ZuoChenLi_2019_Li__MO_732106099012_000 Li
```

### kim-api-collections-management not found

Install OpenKIM API:

```bash
sudo apt update
sudo apt install kim-api
```

## References

- OpenKIM Li-SNAP: `SNAP_ZuoChenLi_2019_Li__MO_732106099012_000` (https://openkim.org)
- LAMMPS KIM Package: https://docs.lammps.org/KIM.html
- LAMMPS pair_style kim: https://docs.lammps.org/pair_kim.html
- Zuo et al. (2019), SNAP Li interatomic potential