"""Readers for LAMMPS output consumed by the SOAP analysis.

[NEW 2026-08-21]

Supported files:

1. LAMMPS data files as written by the v2 generator / write_data
   (orthogonal box, "Atoms # atomic", optional Masses / Velocities).
   -> read_structure()

2. Custom atom dumps as written by the modified stage scripts
   (dump custom ... id type x y z vx vy vz c_peatom c_keatom)
   -> read_peratom_dump()

The readers are deliberately independent of the generator package so the
SOAP environment only needs numpy.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, IO, List, Optional, Tuple

import numpy as np


# ----------------------------------------------------------------------------
# LAMMPS data file
# ----------------------------------------------------------------------------

def read_structure(path, element: str = "Li") -> Dict:
    """Read an orthogonal LAMMPS data file into a structure dict.

    Returns
    -------
    {"ids": (N,) int, "positions": (N, 3) float, "types": (N,) int,
     "cell": (3,) float, "element": str, "natoms": int}
    """
    path = Path(path)
    with path.open() as fh:
        lines = fh.readlines()

    box: List[float] = []
    types: List[int] = []
    rows: List[Tuple[int, float, float, float]] = []
    section: Optional[str] = None
    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        low = line.lower()
        # box bound lines look like "<lo> <hi> xlo xhi" (numbers first)
        if low.endswith("xlo xhi") or low.endswith("ylo yhi") or low.endswith("zlo zhi"):
            box.extend(float(v) for v in line.split()[:2])
            continue
        # section headers: "Atoms # atomic", "Velocities", "Masses", ...
        first_word = low.split()[0]
        if first_word in ("atoms", "velocities", "masses", "bonds",
                          "angles", "dihedrals", "impropers"):
            section = first_word if first_word == "atoms" else None
            continue
        if section == "atoms":
            parts = line.split()
            if len(parts) < 5:
                continue
            rows.append((int(parts[0]), float(parts[2]),
                         float(parts[3]), float(parts[4])))
            types.append(int(parts[1]))

    if len(box) != 6:
        raise ValueError(f"could not parse box bounds from {path}")
    if not rows:
        raise ValueError(f"no atoms parsed from {path}")

    rows.sort(key=lambda t: t[0])
    ids = np.array([r[0] for r in rows], dtype=int)
    pos = np.array([[r[1], r[2], r[3]] for r in rows], dtype=float)
    cell = np.array([box[1] - box[0], box[3] - box[2], box[5] - box[4]],
                    dtype=float)
    # wrap positions into the primary cell (data files can carry slightly
    # out-of-box coordinates after box/relax in x)
    for k in range(3):
        pos[:, k] -= np.floor(pos[:, k] / cell[k]) * cell[k]
    return {"ids": ids, "positions": pos,
            "types": np.array(types, dtype=int),
            "cell": cell, "element": element, "natoms": len(ids)}


# ----------------------------------------------------------------------------
# custom dump with per-atom energies
# ----------------------------------------------------------------------------

def _parse_dump_frames(fh: IO[str]) -> List[Dict]:
    """Parse all frames of a LAMMPS custom dump into dicts."""
    frames: List[Dict] = []
    while True:
        line = fh.readline()
        if not line:
            break
        if not line.startswith("ITEM: TIMESTEP"):
            continue
        timestep = int(fh.readline().split()[0])
        # ITEM: NUMBER OF ATOMS
        natoms = int(fh.readline().split()[0])
        # ITEM: BOX BOUNDS ... followed by 2 (orthogonal) or 3 (triclinic)
        # bound lines, then ITEM: ATOMS <columns>.
        bound_lines: List[str] = []
        atoms_header = ""
        while True:
            bline = fh.readline()
            if not bline:
                raise ValueError("truncated dump: no ITEM: ATOMS header")
            if bline.startswith("ITEM:"):
                atoms_header = bline
                break
            bound_lines.append(bline)
        box = [float(v) for v in bound_lines[0].split()[:2]]
        cols = atoms_header.split()[2:]
        rows = [fh.readline().split() for _ in range(natoms)]
        data = {c: np.array([r[i] for r in rows])
                for i, c in enumerate(cols)}
        frames.append({"timestep": timestep, "columns": cols, "data": data,
                       "box": box})
    return frames


def read_peratom_dump(path, timestep: Optional[int] = None) -> Dict:
    """Read a per-atom custom dump (optionally a specific frame).

    Returns the frame dict {"timestep", "columns", "data"} where "data"
    maps column names (id, type, x, y, z, c_peatom, c_keatom, ...) to
    arrays.  With timestep=None the LAST frame is returned.
    """
    path = Path(path)
    with path.open() as fh:
        frames = _parse_dump_frames(fh)
    if not frames:
        raise ValueError(f"no frames parsed from {path}")
    if timestep is None:
        return frames[-1]
    for fr in frames:
        if fr["timestep"] == timestep:
            return fr
    raise ValueError(f"timestep {timestep} not found in {path}")


def peratom_energies(frame: Dict) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract (ids, per-atom PE [eV], per-atom KE [eV]) from a dump frame.

    Returns columns c_peatom / c_keatom when present; KE falls back to 0.
    """
    data = frame["data"]
    ids = np.array(data["id"], dtype=int)
    pe = np.array(data.get("c_peatom", np.zeros(len(ids))), dtype=float)
    ke = np.array(data.get("c_keatom", np.zeros(len(ids))), dtype=float)
    order = np.argsort(ids)
    return ids[order], pe[order], ke[order]
