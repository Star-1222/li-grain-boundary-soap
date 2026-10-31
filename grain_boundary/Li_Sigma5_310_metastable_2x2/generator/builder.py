"""Bicrystal builder: strict mirror construction of the Sigma5 (310)[001] STGB.

Corrected geometry core
-----------------------
1. Grain 1 is generated from its own rotated BCC lattice
       x = d*(3p + q + 2b),  y = d*(-p + 3q + b),  z = a*(r + b/2)
   inside the half-box  0 <= x <= Lx/2  (p, q, r in Z, b in {0, 1},
   d = a/sqrt(10)).  No ASE, no rotation of a pre-built box.
2. Grain 2 is the exact mirror image of grain 1 across the mid-plane
   x = Lx/2:  (x, y, z) -> (Lx - x, y, z),  with x = Lx wrapped to x = 0.
   The two grains therefore share the same in-plane periods and the two GBs
   (x = 0 and x = Lx/2) are strictly equivalent.
3. Overlap removal only deletes cross-grain coincident pairs that lie on a
   GB plane (|x| < GB_TOL or |x - Lx/2| < GB_TOL) and are closer than
   COINC_TOL.  For the default box there are 60 coincident pairs per GB
   (120 total); one atom per pair is removed, alternating the grain so the
   mirror symmetry is preserved.  Every deletion is logged.

The old generator instead (a) used the BCC primitive cell of
ase.build.bulk() as if it were the conventional cell, (b) rotated the two
half-crystals inside an already built periodic box, (c) deleted every atom
closer than a cutoff anywhere in the system (3600 -> 3224, 10.44 %), and
(d) built the ASE cell with pbc = (False, True, True) while LAMMPS ran with
boundary p p p.  None of that is repeated here.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterator, List, Tuple

import numpy as np

from .config import DEFAULT_A, DEFAULT_J, DEFAULT_K, DEFAULT_L, THETA
from .geometry import box_dims, orientation_matrices, pbc_dist, plane_spacing

GB_TOL = 0.10      # Angstrom: atoms within this distance of a GB plane are "on the GB"
COINC_TOL = 0.05   # Angstrom: max separation of a coincident (overlapping) pair


@dataclass
class Atom:
    """One atom of the final bicrystal (ids renumbered to 1..N)."""

    id: int
    type: int
    x: float
    y: float
    z: float
    grain_id: int
    vx: float = 0.0
    vy: float = 0.0
    vz: float = 0.0


@dataclass
class Bicrystal:
    a: float
    K: int
    J: int
    L: int
    box: Tuple[float, float, float]   # (Lx, Ly, Lz) in Angstrom
    atoms: List[Atom]                  # final atoms, ids 1..N
    removal_log: List[dict]            # one dict per deleted atom
    M1: np.ndarray
    M2: np.ndarray
    theta: float                       # radians
    meta: dict = field(default_factory=dict)


def _mirror_x(x: float, lx: float) -> float:
    """Mirror across x = Lx/2 and wrap x = Lx back to x = 0 (periodic GB1)."""
    xm = lx - x
    if xm >= lx - 1e-9:
        xm -= lx
    return xm


def _iter_grain1_sites(a: float, K: int, J: int, L: int) -> Iterator[Tuple[float, float, float]]:
    """Yield (x, y, z) of every grain-1 BCC site inside the half-box.

    A generous window of lattice indices is swept and the periodic domain
    [0, Lx/2] x [0, Ly) x [0, Lz) is applied.  The same window is used for
    both grains, so the mirrored grain-2 set is an exact mirror image.
    """
    d = plane_spacing(a)
    Lx, Ly, Lz = box_dims(a, K, J, L)
    # Derive a window that provably contains every in-domain site.
    # Reduced units: X = x/d in [0, 5K], Y = y/d in [0, 10J], Z = z/a in [0, L].
    #   p = (3X - Y - 5b)/10 in [-(J+1), 1.5K],  q = (X + 3Y - 5b)/10 in [-1, K/2 + 3J]
    p_lo, p_hi = -(J + 2), math.floor(1.5 * K) + 2
    q_lo, q_hi = -2, math.floor(K / 2.0 + 3 * J) + 2
    r_lo, r_hi = -2, L + 2
    eps = 1e-9
    for p in range(p_lo, p_hi + 1):
        for q in range(q_lo, q_hi + 1):
            for r in range(r_lo, r_hi + 1):
                for b in (0, 1):
                    x = d * (3 * p + q + 2 * b)
                    y = d * (-p + 3 * q + b)
                    z = a * (r + b / 2.0)
                    if -eps <= x <= Lx / 2.0 + eps and -eps <= y < Ly - eps and -eps <= z < Lz - eps:
                        yield (x, y, z)


def _coincident_pairs(raw_g1: List[dict], raw_g2: List[dict],
                      lx: float, ly: float, lz: float) -> dict:
    """Cross-grain coincident pairs sitting on the two GB planes.

    Returns {gb_index: [(g1_atom, g2_atom, distance), ...]} where gb_index
    is 0 for x = 0 (GB1) and 1 for x = Lx/2 (GB2).  Only pairs closer than
    COINC_TOL are returned; exact coincidences have distance 0.
    """
    def plane_of(a: dict) -> int | None:
        if abs(a["x"]) < GB_TOL:
            return 0
        if abs(a["x"] - lx / 2.0) < GB_TOL:
            return 1
        return None

    g2_by: dict = {}
    for b in raw_g2:
        pl = plane_of(b)
        if pl is None:
            continue
        g2_by.setdefault((pl, round(b["y"], 6), round(b["z"], 6)), []).append(b)

    pairs: dict = {0: [], 1: []}
    for a in raw_g1:
        pl = plane_of(a)
        if pl is None:
            continue
        for b in g2_by.get((pl, round(a["y"], 6), round(a["z"], 6)), []):
            dist = pbc_dist(a["x"], a["y"], a["z"], b["x"], b["y"], b["z"], lx, ly, lz)
            if dist < COINC_TOL:
                pairs[pl].append((a, b, dist))
    return pairs


def build_bicrystal(a: float = DEFAULT_A, K: int = DEFAULT_K, J: int = DEFAULT_J,
                    L: int = DEFAULT_L) -> Bicrystal:
    """Build the Sigma5 (310)[001] bicrystal (mirror method)."""
    if a <= 0 or K < 1 or J < 1 or L < 1:
        raise ValueError(f"invalid parameters a={a}, K={K}, J={J}, L={L}")
    Lx, Ly, Lz = box_dims(a, K, J, L)
    M1, M2 = orientation_matrices(THETA)

    sites1 = list(_iter_grain1_sites(a, K, J, L))
    sites2 = [(_mirror_x(x, Lx), y, z) for (x, y, z) in sites1]

    raw_g1 = [{"id": i, "grain_id": 1, "x": x, "y": y, "z": z}
              for i, (x, y, z) in enumerate(sites1, start=1)]
    n1 = len(raw_g1)
    raw_g2 = [{"id": i + n1, "grain_id": 2, "x": x, "y": y, "z": z}
              for i, (x, y, z) in enumerate(sites2, start=1)]

    pairs = _coincident_pairs(raw_g1, raw_g2, Lx, Ly, Lz)

    removal_log: List[dict] = []
    removed_ids: set = set()
    per_plane_removed = {0: {"g1": 0, "g2": 0}, 1: {"g1": 0, "g2": 0}}
    for pl, label in ((0, "GB1"), (1, "GB2")):
        plist = sorted(pairs[pl], key=lambda t: (t[0]["id"], t[1]["id"]))
        if len(plist) % 2 != 0:
            raise RuntimeError(f"{label}: odd number of coincident pairs ({len(plist)})")
        for k, (atom_a, atom_b, dist) in enumerate(plist):
            if dist >= COINC_TOL:
                raise RuntimeError(f"{label}: pair distance {dist} not a coincidence")
            victim, partner = (
                (atom_a, atom_b) if k % 2 == 0 else (atom_b, atom_a)
            )
            removed_ids.add(victim["id"])
            per_plane_removed[pl][f"g{victim['grain_id']}"] += 1
            removal_log.append({
                "atom_id": victim["id"],
                "grain_id": victim["grain_id"],
                "x": victim["x"],
                "y": victim["y"],
                "z": victim["z"],
                "nearest_distance_A": dist,
                "gb": label,
                "partner_atom_id": partner["id"],
            })

    final_atoms: List[Atom] = []
    for raw_atom in list(raw_g1) + list(raw_g2):
        if raw_atom["id"] in removed_ids:
            continue
        final_atoms.append(
            Atom(
                id=0,
                type=1,
                x=raw_atom["x"],
                y=raw_atom["y"],
                z=raw_atom["z"],
                grain_id=raw_atom["grain_id"],
            )
        )
    # deterministic order: all grain-1 atoms first, then all grain-2 atoms
    final_atoms.sort(key=lambda at: (at.grain_id, at.x, at.y, at.z))
    for idx, at in enumerate(final_atoms, start=1):
        at.id = idx

    n_removed_g1 = sum(per_plane_removed[p]["g1"] for p in (0, 1))
    n_removed_g2 = sum(per_plane_removed[p]["g2"] for p in (0, 1))
    if n_removed_g1 != n_removed_g2:
        raise RuntimeError(f"asymmetric removal: g1={n_removed_g1} vs g2={n_removed_g2}")

    meta = {
        "a_A": a,
        "box_A": [Lx, Ly, Lz],
        "in_plane_periods_A": [math.sqrt(10.0) * a, a],
        "K": K, "J": J, "L": L,
        "theta_deg": math.degrees(THETA),
        "sites_per_grain": n1,
        "total_before_removal": len(raw_g1) + len(raw_g2),
        "coincident_pairs_gb1": len(pairs[0]),
        "coincident_pairs_gb2": len(pairs[1]),
        "removed_total": len(removal_log),
        "removed_grain1": n_removed_g1,
        "removed_grain2": n_removed_g2,
        "final_total": len(final_atoms),
        "final_grain1": sum(1 for at in final_atoms if at.grain_id == 1),
        "final_grain2": sum(1 for at in final_atoms if at.grain_id == 2),
        "removal_fraction_pct": 100.0 * len(removal_log) / (len(raw_g1) + len(raw_g2)),
    }

    return Bicrystal(a=a, K=K, J=J, L=L, box=(Lx, Ly, Lz), atoms=final_atoms,
                     removal_log=removal_log, M1=M1, M2=M2, theta=THETA, meta=meta)
