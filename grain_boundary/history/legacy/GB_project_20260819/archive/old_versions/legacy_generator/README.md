# Li Bicrystal Generator

A reusable, modular Python framework for generating **BCC lithium bicrystal models** for LAMMPS + OpenKIM simulations, built on [ASE](https://wiki.fysik.dtu.dk/ase/) (Atomic Simulation Environment).

Designed for computational materials science research on lithium metal grain boundaries. The architecture supports future extension to other crystal structures (FCC, HCP), elements, and grain boundary types without modifying existing code.

---

## Project Structure

```
Li_Bicrystal_Generator/
│
├── __init__.py          Package init, version, public API
├── __main__.py          CLI:  python -m Li_Bicrystal_Generator
├── main.py              Argparse CLI entry point
├── config.py            All physical parameters + GB database
├── builder.py           Core BiCrystalBuilder class
├── geometry.py          Rotation matrices, overlap removal, centering
├── exporter.py          LAMMPS data / XYZ / CIF export
├── utils.py             Validation, GB lookup helpers
│
├── output/              Generated bicrystal files go here
├── examples/            Ready-to-run example scripts
│   └── example.py
│
└── README.md            This file
```

---

## Installation

### 1. Install ASE

```bash
pip install ase numpy
```

ASE is the only hard dependency. NumPy is used internally for matrix operations.

### 2. Verify

```bash
python -c "import ase; print(ase.__version__)"
```

---

## Usage

### Command line

```bash
# List available grain boundaries
python -m Li_Bicrystal_Generator --list

# Generate Sigma5 (310) bicrystal
python -m Li_Bicrystal_Generator --name Sigma5_310

# Override parameters
python -m Li_Bicrystal_Generator --name Sigma13_510 --lattice 3.52 --supercell 30 10 10

# Custom output directory and formats
python -m Li_Bicrystal_Generator --name Sigma5_210 --output-dir my_gbs --formats lammps xyz
```

### As a Python library

```python
from Li_Bicrystal_Generator.builder import BiCrystalBuilder, GrainBoundaryDescriptor
from Li_Bicrystal_Generator.config import get_default_config
from Li_Bicrystal_Generator.exporter import export

config = get_default_config()
config["lattice_constant"] = 3.51
config["supercell"] = (20, 12, 12)
config["overlap_cutoff"] = 2.3

gb = GrainBoundaryDescriptor(
    name="Sigma5_310",
    sigma=5, axis=(0, 0, 1), angle=36.87, plane=(3, 1, 0),
)

builder = BiCrystalBuilder(config)
atoms = builder.build(gb)

export(atoms, gb.name, output_dir="output", formats=("lammps", "xyz", "cif"))
```

### Running the example

```bash
python examples/example.py
```

Output files appear in `output/`.

---

## Grain Boundary Database

All available grain boundaries are defined in `config.py` in the `GRAIN_BOUNDARIES` list.  Currently implemented (all [001] tilt GBs in BCC):

| Name          | Sigma | Axis   | Angle (°) | Plane  |
|---------------|-------|--------|-----------|--------|
| Sigma5_310    | 5     | [001]  | 36.87     | (310)  |
| Sigma5_210    | 5     | [001]  | 53.13     | (210)  |
| Sigma13_510   | 13    | [001]  | 22.62     | (510)  |
| Sigma17_410   | 17    | [001]  | 28.07     | (410)  |
| Sigma17_530   | 17    | [001]  | 61.93     | (530)  |
| Sigma29_520   | 29    | [001]  | 43.60     | (520)  |

Only physically meaningful CSL grain boundaries are included. Random generation is not supported.

---

## Export Formats

| Format     | Extension | Description                           |
|------------|-----------|---------------------------------------|
| LAMMPS     | `.lmp`    | `read_data` format, ready for LAMMPS  |
| XYZ        | `.xyz`    | Human-readable, OVITO/VMD compatible  |
| CIF        | `.cif`    | Crystallographic Information File     |

Adding a new format requires:
1. Write an export function in `exporter.py` (e.g. `export_vasp`)
2. Register it in the `FORMAT_MAP` dictionary

No other code changes are needed.

---

## Visualisation with OVITO

1. Open the `.xyz` or `.lmp` file in OVITO.
2. Color by particle type (should be all Li).
3. Use the Slice or Identify Diamond tools to highlight the GB region.
4. The GB plane lies perpendicular to the x-axis (centered in the cell).

---

## Preparing for LAMMPS + OpenKIM (Li-SNAP)

The exported `.lmp` file can be used directly in a LAMMPS input script:

```lammps
units           metal
atom_style      atomic
boundary        p p p

read_data       Sigma5_310.lmp

pair_style      kim SNAP_ZuoChenLi_2019_Li__MO_732106099012_000
pair_coeff      * * Li

minimize        1.0e-8 1.0e-8 10000 100000
```

Note: the exported bicrystal uses `atom_style atomic` by default (element type 1 = Li). Ensure this matches your LAMMPS pair style setup.

---

## Architecture & Extensibility

### Adding a new grain boundary

Edit `config.py` and append to `GRAIN_BOUNDARIES`:

```python
{
    "name": "MyGB",
    "sigma": 37,
    "axis": (0, 0, 1),
    "angle": 18.92,
    "plane": (6, 1, 0),
}
```

No other code changes are required.

### Adding a new crystal structure

In `builder.py`, the `_create_oriented_cell` method uses `ase.build.bulk()`.  Changing the `crystal` config parameter (e.g. to `"fcc"`) will use the corresponding ASE lattice builder.  The orientation logic in `geometry.py` may need adjustment for non-cubic systems.

### Adding a new export format

Write a function in `exporter.py` and register it in `FORMAT_MAP`.

### Batch / high-throughput generation

The `BiCrystalBuilder` class is stateless between `build()` calls.  Loop over grain boundaries in a script:

```python
for gb_data in GRAIN_BOUNDARIES:
    gb = GrainBoundaryDescriptor.from_dict(gb_data)
    atoms = builder.build(gb)
    export(atoms, gb.name, ...)
```

---

## Dependencies

- Python >= 3.10
- ASE >= 3.22
- NumPy

No GUI libraries, visualization tools, or external potential files are required.

---

## License

Academic use.  No warranty.

---

## References

- Sutton & Balluffi, *Interfaces in Crystalline Materials* (1995)
- OpenKIM: https://openkim.org
- ASE: https://wiki.fysik.dtu.dk/ase/
- LAMMPS KIM package: https://docs.lammps.org/KIM.html
- R. Zuo et al., *SNAP Li interatomic potential* (2019)