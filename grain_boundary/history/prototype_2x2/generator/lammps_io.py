"""Minimal readers for LAMMPS data files and grain-id maps.

Used by the re-import verification (verify_data_file) so the exported
files are checked independently of the in-memory builder object.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

from .builder import Atom


def read_lammps_data(path: str | Path) -> Tuple[List[Atom], Tuple[float, float, float], int]:
    """Parse an orthogonal LAMMPS data file (atom_style atomic).

    Returns (atoms, box, ntypes).  Atom velocities are read if present.
    """
    path = Path(path)
    lines = path.read_text(encoding="utf-8").splitlines()
    natoms = None
    ntypes = 1
    box = None
    section = None
    atoms: List[Atom] = []
    velocities: Dict[int, Tuple[float, float, float]] = {}
    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        low = line.lower()
        if low.endswith("atoms") and len(line.split()) == 2 and line.split()[0].isdigit():
            natoms = int(line.split()[0])
            continue
        if low.startswith("atom types"):
            ntypes = int(line.split()[0])
            continue
        if low.endswith("xlo xhi") or (len(line.split()) == 3 and low.endswith("xhi")):
            if box is None:
                parts = line.split()
                xlo, xhi = float(parts[0]), float(parts[1])
                continue
        if low.endswith("ylo yhi"):
            parts = line.split()
            ylo, yhi = float(parts[0]), float(parts[1])
            continue
        if low.endswith("zlo zhi"):
            parts = line.split()
            zlo, zhi = float(parts[0]), float(parts[1])
            box = (xhi - xlo, yhi - ylo, zhi - zlo)
            continue
        if low in ("atoms", "atoms # atomic", "atoms # full"):
            section = "atoms"
            continue
        if low.startswith("velocities"):
            section = "velocities"
            continue
        if low in ("masses", "bonds", "angles", "dihedrals", "impropers"):
            section = None
            continue
        if section == "atoms":
            parts = line.split()
            if len(parts) >= 5:
                aid = int(parts[0])
                atype = int(parts[1])
                atoms.append(Atom(id=aid, type=atype, x=float(parts[2]),
                                 y=float(parts[3]), z=float(parts[4]), grain_id=0))
            continue
        if section == "velocities":
            parts = line.split()
            if len(parts) >= 4:
                velocities[int(parts[0])] = (float(parts[1]), float(parts[2]), float(parts[3]))
    if box is None:
        raise ValueError(f"no box bounds found in {path}")
    for at in atoms:
        v = velocities.get(at.id)
        if v is not None:
            at.vx, at.vy, at.vz = v
    if natoms is not None and natoms != len(atoms):
        raise ValueError(f"{path}: header says {natoms} atoms, parsed {len(atoms)}")
    return atoms, box, ntypes


def read_grain_map(path: str | Path) -> Dict[int, int]:
    """Read an atom_id -> grain_id map (whitespace separated)."""
    path = Path(path)
    mapping: Dict[int, int] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            mapping[int(parts[0])] = int(parts[1])
    return mapping
