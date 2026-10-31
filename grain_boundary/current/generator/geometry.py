"""Crystallographic core of the current Sigma5 (310)[001] bicrystal generator.

Everything is expressed with explicit Miller indices and Cartesian vectors.
No ASE is required: the conventional BCC cell, the orientation matrices and
the mirror construction are written out in closed form. ASE is only used
(optionally) in the tests for cross-validation of the exported files.
"""
from __future__ import annotations

import math
from typing import Iterable, Iterator, Sequence

import numpy as np


def rotz(angle: float) -> np.ndarray:
    """Rotation matrix Rz(angle); rows = new axes expressed in old basis."""
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def conventional_bcc_basis() -> np.ndarray:
    """BCC conventional cubic cell: two atoms at (0,0,0) and (1/2,1/2,1/2).

    Returns fractional coordinates (2, 3). This is the correct cell for the
    Miller arithmetic below; the legacy generator mistakenly used the
    primitive cell of ase.build.bulk("Li", "bcc") as if it were conventional.
    """
    return np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])


def orientation_matrices(theta: float) -> tuple[np.ndarray, np.ndarray]:
    """Orientation matrices of the two grains (rows = crystal -> lab).

    M1 = Rz(+theta/2)  (grain 1),  M2 = Rz(-theta/2)  (grain 2).
    The two grains are mirror images across the GB plane; their relative
    rotation is M2 . M1^T = Rz(-theta), i.e. 36.87 deg about [001].
    """
    return rotz(theta / 2.0), rotz(-theta / 2.0)


def box_dims(a: float, k: int, j: int, l: int) -> tuple[float, float, float]:
    """Orthogonal box (Lx, Ly, Lz). Lx = k*sqrt(10)*a, Ly = j*sqrt(10)*a, Lz = l*a."""
    d = a / math.sqrt(10.0)
    return k * math.sqrt(10.0) * a, j * math.sqrt(10.0) * a, l * a


def plane_spacing(a: float) -> float:
    """d = a/sqrt(10): spacing between (310) planes along the GB normal."""
    return a / math.sqrt(10.0)


def grain1_position(p: int, q: int, r: int, b: int, a: float, d: float) -> tuple[float, float, float]:
    """Lab position of grain-1 lattice point (p, q, r, b).

    x = d(3p + q + 2b),  y = d(-p + 3q + b),  z = a(r + b/2).
    """
    return (d * (3 * p + q + 2 * b), d * (-p + 3 * q + b), a * (r + b / 2.0))


def lattice_indices_grain1(x: float, y: float, z: float, a: float, d: float) -> tuple[int, int, int, int] | None:
    """Inverse of grain1_position: return integer (p, q, r, b) or None."""
    xd, yd = x / d, y / d
    for b in (0, 1):
        p = (3 * xd - yd - 5 * b) / 10.0
        q = (xd + 3 * yd - 5 * b) / 10.0
        r = z / a - b / 2.0
        if abs(p - round(p)) < 1e-8 and abs(q - round(q)) < 1e-8 and abs(r - round(r)) < 1e-8:
            return round(p), round(q), round(r), b
    return None


def grain2_position(p: int, q: int, r: int, b: int, a: float, d: float) -> tuple[float, float, float]:
    """Lab position of grain-2 lattice point (M2 = Rz(-theta/2)).

    x = d(3p - q + b),  y = d(p + 3q + 2b),  z = a(r + b/2).
    """
    return (d * (3 * p - q + b), d * (p + 3 * q + 2 * b), a * (r + b / 2.0))


def lattice_indices_grain2(x: float, y: float, z: float, a: float, d: float) -> tuple[int, int, int, int] | None:
    """Inverse of grain2_position: return integer (p, q, r, b) or None."""
    xd, yd = x / d, y / d
    for b in (0, 1):
        p = (3 * xd + yd - 5 * b) / 10.0
        q = (3 * yd - xd - 5 * b) / 10.0
        r = z / a - b / 2.0
        if abs(p - round(p)) < 1e-8 and abs(q - round(q)) < 1e-8 and abs(r - round(r)) < 1e-8:
            return round(p), round(q), round(r), b
    return None


def wrap(x: float, lo: float, hi: float) -> float:
    """Wrap x into [lo, hi)."""
    period = hi - lo
    return lo + (x - lo) % period


def pbc_vector(ax: float, ay: float, az: float, bx: float, by: float, bz: float,
               lx: float, ly: float, lz: float) -> tuple[float, float, float]:
    """Minimum-image vector b - a under 3D periodicity."""
    dx = bx - ax
    dy = by - ay
    dz = bz - az
    dx -= math.floor(dx / lx + 0.5) * lx
    dy -= math.floor(dy / ly + 0.5) * ly
    dz -= math.floor(dz / lz + 0.5) * lz
    return dx, dy, dz


def pbc_dist(ax: float, ay: float, az: float, bx: float, by: float, bz: float,
             lx: float, ly: float, lz: float) -> float:
    dx, dy, dz = pbc_vector(ax, ay, az, bx, by, bz, lx, ly, lz)
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def rotation_angle_axis(M: np.ndarray) -> tuple[float, np.ndarray]:
    """Angle (deg) and unit axis of rotation matrix M."""
    trace = float(np.trace(M))
    cosang = max(-1.0, min(1.0, (trace - 1.0) / 2.0))
    angle = math.degrees(math.acos(cosang))
    skew = np.array([M[2, 1] - M[1, 2], M[0, 2] - M[2, 0], M[1, 0] - M[0, 1]])
    norm = float(np.linalg.norm(skew))
    axis = skew / norm if norm > 1e-12 else np.array([0.0, 0.0, 1.0])
    return angle, axis


def cubic_rotation_group() -> list[np.ndarray]:
    """The 24 proper rotations of the cubic point group (m3m)."""
    perms = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
    signs = [(1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1),
             (1, -1, -1), (-1, 1, -1), (-1, -1, 1), (-1, -1, -1)]
    group = []
    for p in perms:
        for s in signs:
            m = np.zeros((3, 3))
            for i in range(3):
                m[i, p[i]] = s[i]
            if np.linalg.det(m) > 0.5:
                group.append(m)
    return group


def disorientation_angle_deg(M: np.ndarray, group: Sequence[np.ndarray] | None = None) -> float:
    """Minimum rotation angle (deg) of M*g over the cubic group."""
    group = cubic_rotation_group() if group is None else group
    best = 180.0
    for g in group:
        angle, _ = rotation_angle_axis(M @ g)
        best = min(best, angle)
    return best


def fit_orientation(center: tuple[float, float, float], positions: Sequence[tuple[float, float, float]],
                    ref: np.ndarray, lx: float, ly: float, lz: float, a: float) -> np.ndarray:
    """Fit the cubic axes of a perfect BCC neighbourhood around `center`.

    The six neighbours at distance a form three +/- pairs; the canonical
    signed permutation that best matches the reference matrix `ref` is kept.
    Returns the fitted orientation matrix (rows = crystal -> lab).
    """
    cx, cy, cz = center
    nbrs = []
    for p in positions:
        dist = pbc_dist(cx, cy, cz, p[0], p[1], p[2], lx, ly, lz)
        if 0.98 * a < dist < 1.02 * a:
            nbrs.append(pbc_vector(cx, cy, cz, p[0], p[1], p[2], lx, ly, lz))
    axes: list[np.ndarray] = []
    used = [False] * len(nbrs)
    for i, vi in enumerate(nbrs):
        if used[i]:
            continue
        for j in range(i + 1, len(nbrs)):
            if used[j]:
                continue
            if abs(sum(vi[k] * nbrs[j][k] for k in range(3)) + a * a) < 0.01:
                used[i] = used[j] = True
                axes.append(np.array(vi) / a)
                break
    if len(axes) != 3:
        raise ValueError(f"expected 3 cubic axes, found {len(axes)}")
    perms = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
    signs = [(1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1), (-1, -1, -1)]
    best: np.ndarray | None = None
    best_tr = -math.inf
    for p in perms:
        for s in signs:
            rows = np.array([axes[p[0]] * s[0], axes[p[1]] * s[1], axes[p[2]] * s[2]])
            if np.linalg.det(rows) < 0:
                continue
            tr = float(np.sum(ref * rows))
            if tr > best_tr:
                best_tr = tr
                best = rows
    return best


class CellList:
    """Simple linked-cell neighbour search with minimum-image periodicity."""

    def __init__(self, positions: Sequence[tuple[float, float, float]],
                 lx: float, ly: float, lz: float, rcut: float):
        self.lx, self.ly, self.lz = lx, ly, lz
        self.rcut = rcut
        self.nx = max(1, int(lx // rcut))
        self.ny = max(1, int(ly // rcut))
        self.nz = max(1, int(lz // rcut))
        self.cells: dict[tuple[int, int, int], list[int]] = {}
        self.positions = list(positions)
        for idx, pos in enumerate(self.positions):
            key = self._cell_key(pos)
            self.cells.setdefault(key, []).append(idx)

    def _cell_key(self, pos: tuple[float, float, float]) -> tuple[int, int, int]:
        return (int(pos[0] // self.rcut) % self.nx,
                int(pos[1] // self.rcut) % self.ny,
                int(pos[2] // self.rcut) % self.nz)

    def neighbors(self, idx: int, rmax: float) -> Iterator[tuple[int, float]]:
        """Yield (neighbour_index, distance) within rmax of atom idx."""
        pos = self.positions[idx]
        cx, cy, cz = self._cell_key(pos)
        r2 = rmax * rmax
        for ix in (-1, 0, 1):
            for iy in (-1, 0, 1):
                for iz in (-1, 0, 1):
                    key = ((cx + ix) % self.nx, (cy + iy) % self.ny, (cz + iz) % self.nz)
                    for j in self.cells.get(key, []):
                        if j == idx:
                            continue
                        d2 = 0.0
                        for comp in range(3):
                            delta = self.positions[j][comp] - pos[comp]
                            box = (self.lx, self.ly, self.lz)[comp]
                            delta -= math.floor(delta / box + 0.5) * box
                            d2 += delta * delta
                        if d2 < r2:
                            yield j, math.sqrt(d2)

    def min_distance(self) -> float:
        """Global minimum pair distance (pbc)."""
        best = math.inf
        for i in range(len(self.positions)):
            for _, dist in self.neighbors(i, best):
                if dist < best:
                    best = dist
        return best
