"""Export utilities for bicrystal structures.

Each format gets a dedicated function.  Adding a new format requires
only writing one function and registering it in the FORMAT_MAP below.
"""

from __future__ import annotations

import os
from typing import Callable

from ase import Atoms
from ase.io import write as ase_write


# ── dedicated exporters ─────────────────────────────────────────────

def export_lammps(atoms: Atoms, filename: str) -> None:
    """Export to LAMMPS data format.

    The file can be read directly with ``read_data`` in a LAMMPS input
    script, or visualised with OVITO / VMD.

    Parameters
    ----------
    atoms : Atoms
        Structure to export.
    filename : str
        Output path (e.g. ``"Sigma5_310.lmp"``).
    """
    ase_write(filename, atoms, format="lammps-data")
    print(f"  [export] LAMMPS data → {filename}")


def export_xyz(atoms: Atoms, filename: str) -> None:
    """Export to XYZ format (human-readable, compatible with OVITO)."""
    ase_write(filename, atoms, format="xyz")
    print(f"  [export] XYZ         → {filename}")


def export_cif(atoms: Atoms, filename: str) -> None:
    """Export to CIF format (crystallographic information file)."""
    ase_write(filename, atoms, format="cif")
    print(f"  [export] CIF         → {filename}")


# ── grain_id XYZ for OVITO ──────────────────────────────────────────

def export_grainid_xyz(atoms: Atoms, filename: str) -> None:
    """Export an extended XYZ file that includes the ``grain_id`` array.

    The file uses the extended XYZ format recognised by OVITO:
    the ``grain_id`` column appears as a property so that OVITO's
    Color Coding modifier can colour atoms by grain membership.

    Parameters
    ----------
    atoms : Atoms
        Structure with a ``grain_id`` integer array attached.
    filename : str
        Output path (e.g. ``"Sigma5_310_grainID.xyz"``).
    """
    if "grain_id" not in atoms.arrays:
        print("  [export] WARNING: no 'grain_id' array found — writing plain XYZ")
        ase_write(filename, atoms, format="xyz")
        print(f"  [export] XYZ         → {filename}")
        return

    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    grain_ids = atoms.arrays["grain_id"]
    cell = atoms.cell
    n = len(atoms)

    # Build the OVITO-compatible comment line
    lattice_str = (
        f"{cell[0,0]:.10f} {cell[0,1]:.10f} {cell[0,2]:.10f} "
        f"{cell[1,0]:.10f} {cell[1,1]:.10f} {cell[1,2]:.10f} "
        f"{cell[2,0]:.10f} {cell[2,1]:.10f} {cell[2,2]:.10f}"
    )
    comment = (
        f'Lattice="{lattice_str}" '
        f'Properties=species:S:1:pos:R:3:grain_id:R:1'
    )

    with open(filename, "w") as f:
        f.write(f"{n}\n")
        f.write(f"{comment}\n")
        for sym, (x, y, z), gid in zip(symbols, positions, grain_ids):
            f.write(f"{sym} {x:.10f} {y:.10f} {z:.10f} {gid}\n")

    print(f"  [export] XYZ (grain_id) → {filename}")


# ── dispatch map (add new formats here) ─────────────────────────────

FORMAT_MAP: dict[str, Callable[[Atoms, str], None]] = {
    "lammps": export_lammps,
    "xyz": export_xyz,
    "cif": export_cif,
    "grainid_xyz": export_grainid_xyz,
}

SUPPORTED_FORMATS: tuple[str, ...] = tuple(FORMAT_MAP.keys())


# ── generic entry point ─────────────────────────────────────────────

# ── special extension mapping ────────────────────────────────────────
# Formats whose file name differs from simply ``{base_name}.{ext}``.
_SPECIAL_EXT: dict[str, str] = {
    "grainid_xyz": "_grainID.xyz",
}


def export(
    atoms: Atoms,
    base_name: str,
    output_dir: str = "output",
    formats: tuple[str, ...] = ("lammps", "xyz", "cif", "grainid_xyz"),
) -> None:
    """Export *atoms* to one or more file formats.

    Parameters
    ----------
    atoms : Atoms
        Structure to export.
    base_name : str
        Base file name (without extension), e.g. ``"Sigma5_310"``.
    output_dir : str
        Directory for output files (created if it does not exist).
    formats : tuple of str
        Which formats to write.  Supported: ``"lammps"``, ``"xyz"``,
        ``"cif"``, ``"grainid_xyz"``.
    """
    os.makedirs(output_dir, exist_ok=True)

    for fmt in formats:
        if fmt not in FORMAT_MAP:
            print(f"  [export] WARNING: unknown format '{fmt}', skipping")
            continue

        if fmt in _SPECIAL_EXT:
            path = os.path.join(output_dir, f"{base_name}{_SPECIAL_EXT[fmt]}")
        else:
            ext = {"lammps": "lmp"}.get(fmt, fmt)
            path = os.path.join(output_dir, f"{base_name}.{ext}")
        FORMAT_MAP[fmt](atoms, path)