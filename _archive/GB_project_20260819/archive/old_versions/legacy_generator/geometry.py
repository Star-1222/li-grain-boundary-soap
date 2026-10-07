"""Geometric utilities for bicrystal construction.

Provides rotation matrices, orientation helpers, and overlap removal.
All operations use NumPy only (no external geometry libraries).
"""

from __future__ import annotations

import numpy as np
from ase import Atoms
from ase.neighborlist import primitive_neighbor_list


def rotation_matrix_around_axis(
    axis: np.ndarray, angle: float, deg: bool = False
) -> np.ndarray:
    """Build a 3x3 rotation matrix for a rotation around *axis* by *angle*.

    Uses Rodrigues' rotation formula.

    Parameters
    ----------
    axis : ndarray, shape (3,)
        Rotation axis (does not need to be unit length).
    angle : float
        Rotation angle.
    deg : bool, default False
        If True, *angle* is in degrees; otherwise radians.

    Returns
    -------
    R : ndarray, shape (3, 3)
        Rotation matrix such that ``v_rot = R @ v``.
    """
    if deg:
        angle = np.radians(angle)

    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)

    c = np.cos(angle)
    s = np.sin(angle)
    t = 1.0 - c

    x, y, z = axis

    R = np.array([
        [t*x*x + c,    t*x*y - s*z,  t*x*z + s*y],
        [t*x*y + s*z,  t*y*y + c,    t*y*z - s*x],
        [t*x*z - s*y,  t*y*z + s*x,  t*z*z + c  ],
    ])
    return R


def orientation_matrix(h: int, k: int, l: int) -> np.ndarray:
    """Build the 3x3 rotation matrix that orients a cubic cell for GB construction.

    The resulting frame has:
        x  —— GB plane normal (h, k, l) direction
        y  —— in-plane direction perpendicular to the tilt axis
        z  —— tilt axis direction [001]

    This function is designed for [001] tilt grain boundaries in cubic crystals.

    Parameters
    ----------
    h, k, l : int
        Miller indices of the GB plane.  For [001] tilt, l = 0.

    Returns
    -------
    R : ndarray, shape (3, 3)
        Orthonormal rotation matrix (rows: new basis vectors in old coordinates).
    """
    n = np.array([float(h), float(k), float(l)], dtype=float)
    n_norm = n / np.linalg.norm(n)

    t = np.array([0.0, 0.0, 1.0])

    u = np.cross(t, n)
    u_norm = u / np.linalg.norm(u)

    R = np.vstack([n_norm, u_norm, t])
    return R


def remove_overlapping_atoms(
    atoms: Atoms, cutoff: float
) -> Atoms:
    """Remove atoms that are closer than *cutoff* to another atom.

    Uses ASE's cell-list-based neighbor search (O(N log N)) instead of
    a brute-force O(N^2) loop, making it suitable for large systems.

    Parameters
    ----------
    atoms : Atoms
        Input atomic structure.
    cutoff : float
        Minimum allowed interatomic distance in Angstrom.

    Returns
    -------
    atoms : Atoms
        Structure with overlapping atoms removed.
    """
    positions = atoms.get_positions()
    cell = atoms.get_cell()
    pbc = atoms.pbc

    # Use ASE's primitive neighbor list (cell-list based, O(N log N))
    i_idx, j_idx, *_ = primitive_neighbor_list(
        "ijd", pbc, cell, positions, cutoff=cutoff,
    )

    # Mark atoms that have a neighbor within cutoff
    # Keep the lower index, remove the higher index
    remove = set()
    for i, j in zip(i_idx, j_idx):
        if i != j:
            remove.add(max(i, j))

    mask = np.ones(len(atoms), dtype=bool)
    for idx in remove:
        mask[idx] = False

    return atoms[mask]


def center_gb(atoms: Atoms, shift: float = 0.0) -> Atoms:
    """Translate atoms so the GB plane is centred in the cell along x.

    Parameters
    ----------
    atoms : Atoms
        Input structure with GB region somewhere in the cell.
    shift : float, default 0.0
        Additional fractional shift in x (between -0.5 and 0.5).

    Returns
    -------
    atoms : Atoms
        Structure with positions adjusted.
    """
    frac = atoms.get_scaled_positions()
    x_frac = frac[:, 0]
    x_mid = 0.5 * (x_frac.max() + x_frac.min())
    frac[:, 0] -= x_mid - 0.5 + shift
    frac[:, 0] = np.mod(frac[:, 0], 1.0)
    atoms.set_scaled_positions(frac)
    atoms.wrap()
    return atoms