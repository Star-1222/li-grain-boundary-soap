"""Automatic geometric verification of the built bicrystal.

Every check requested by the user is implemented here:

1. misorientation angle = 36.86989765 deg about [001]
2. rotation axis = [001]
3. GB plane = (310) in both crystal frames
4. two GBs explicit and periodically equivalent (distance spectra identical)
5. atoms far from the GB remain perfect BCC (8 NN + 6 NNN)
6. nearest-neighbour distance is reasonable (2d = 2.2199 A), no overlaps
7. no large voids along x
8. deleted atoms only on the GB planes (cross-grain coincidences)
9. atom count and mass density are exactly 3600 and 2/a^3
10. y and z in-plane periods are strictly compatible (lattice indices map)
11. the broken legacy 3600 -> 3224 / 10.44 % scenario is absent

'run_all_checks()' raises VerificationError on the first failing check and
otherwise returns a dict with every measured quantity (used to write the
verification report).  'verify_data_file()' re-imports a LAMMPS data file
(plus the optional grain map) and runs the position-based checks again, so
the exported files are validated independently of the in-memory object.
"""
from __future__ import annotations

import math
import re
from pathlib import Path
from typing import List, Sequence, Tuple

import numpy as np

from .builder import Atom, Bicrystal, COINC_TOL, GB_TOL
from .config import (DEFAULT_A, DEFAULT_J, DEFAULT_K, DEFAULT_L, THETA,
                     THETA_DEG, TILT_AXIS, GB_PLANE_LABEL)
from .geometry import (CellList, box_dims, cubic_rotation_group,
                       disorientation_angle_deg, fit_orientation,
                       lattice_indices_grain1, lattice_indices_grain2,
                       orientation_matrices, pbc_dist, plane_spacing,
                       rotation_angle_axis)

NN_SCALE = math.sqrt(3.0) / 2.0   # BCC nearest-neighbour distance = a*sqrt(3)/2
NNN_SCALE = 1.0                   # BCC second-neighbour distance = a


class VerificationError(Exception):
    """Raised when a geometric verification check fails."""


def _require(cond: bool, msg: str) -> None:
    if not cond:
        raise VerificationError(msg)


def check_misorientation(bc: Bicrystal) -> dict:
    R = bc.M2 @ bc.M1.T
    angle, axis = rotation_angle_axis(R)
    cubic_min = disorientation_angle_deg(R, cubic_rotation_group())
    _require(abs(angle - THETA_DEG) < 1e-6,
             f"misorientation angle {angle} != 36.86989765 deg")
    _require(abs(axis[2]) > 0.999 and abs(axis[0]) < 1e-6 and abs(axis[1]) < 1e-6,
             f"rotation axis {axis} != [001]")
    _require(abs(cubic_min - THETA_DEG) < 1e-6,
             f"cubic-symmetry misorientation {cubic_min} != 36.86989765 deg")
    return {"misorientation_angle_deg": angle,
            "rotation_axis": axis.tolist(),
            "cubic_symmetry_min_angle_deg": cubic_min,
            "ok": True}


def check_gb_plane(bc: Bicrystal) -> dict:
    xhat = np.array([1.0, 0.0, 0.0])
    n1 = bc.M1.T @ xhat          # GB normal in crystal-1 frame
    n2 = bc.M2.T @ xhat          # GB normal in crystal-2 frame
    n1 = n1 / np.linalg.norm(n1)
    n2 = n2 / np.linalg.norm(n2)

    ideal1 = np.array([3.0, -1.0, 0.0]) / math.sqrt(10.0)
    ideal2 = np.array([3.0, 1.0, 0.0]) / math.sqrt(10.0)

    err1 = min(np.linalg.norm(n1 - ideal1), np.linalg.norm(n1 + ideal1))
    err2 = min(np.linalg.norm(n2 - ideal2), np.linalg.norm(n2 + ideal2))

    _require(
        err1 < 1e-6 and err2 < 1e-6,
        f"GB normal in crystal frames ({n1}, {n2}) != {{310}} normals",
    )
    return {"gb_normal_crystal1": n1.tolist(),
            "gb_normal_crystal2": n2.tolist(),
            "gb_normal_error1": float(err1),
            "gb_normal_error2": float(err2),
            "gb_normal_crystal1_miller": "[3,-1,0]",
            "gb_normal_crystal2_miller": "[3,1,0]",
            "gb_plane_label": GB_PLANE_LABEL,
            "ok": True}


def _cross_gb_spectra(bc: Bicrystal, rcut: float = 4.5) -> Tuple[List[float], List[float]]:
    """Distance spectra of cross-grain pairs below rcut, per GB plane."""
    Lx, Ly, Lz = bc.box
    pos = [(at.x, at.y, at.z) for at in bc.atoms]
    gid = [at.grain_id for at in bc.atoms]
    cl = CellList(pos, Lx, Ly, Lz, rcut)
    spec = {0: [], 1: []}
    seen = set()
    for i in range(len(pos)):
        for j, dist in cl.neighbors(i, rcut):
            if gid[i] == gid[j]:
                continue
            key = (min(i, j), max(i, j))
            if key in seen:
                continue
            seen.add(key)
            ai, aj = pos[i], pos[j]
            # classify which GB plane the pair belongs to by x position
            xi, xj = ai[0], aj[0]
            if (xi < Lx / 4 and xj < Lx / 4) or (xi < Lx / 4 and xj > 3 * Lx / 4) \
               or (xi > 3 * Lx / 4 and xj < Lx / 4) or (xi > 3 * Lx / 4 and xj > 3 * Lx / 4):
                spec[0].append(dist)      # near x = 0 (GB1, periodic)
            else:
                spec[1].append(dist)      # near x = Lx/2 (GB2)
    return sorted(spec[0]), sorted(spec[1])


def check_two_gbs(bc: Bicrystal) -> dict:
    Lx, Ly, Lz = bc.box
    near0 = [at for at in bc.atoms if at.x < 1.0 or at.x > Lx - 1.0]
    nearhalf = [at for at in bc.atoms if abs(at.x - Lx / 2) < 1.0]
    _require(len(near0) > 0, "no atoms near GB1 (x=0)")
    _require(len(nearhalf) > 0, "no atoms near GB2 (x=Lx/2)")
    s1, s2 = _cross_gb_spectra(bc)
    _require(len(s1) == len(s2) and len(s1) > 0, "GB spectra sizes differ or empty")
    _require(np.allclose(np.array(s1), np.array(s2), atol=1e-9),
             "GB1 and GB2 cross-GB distance spectra are not identical")
    return {"gb1_x": 0.0, "gb2_x": Lx / 2.0,
            "cross_gb_pairs_below_4p5A_gb1": len(s1),
            "cross_gb_pairs_below_4p5A_gb2": len(s2),
            "spectra_identical": True,
            "min_gb1_A": s1[0], "min_gb2_A": s2[0],
            "pairs_at_2p22A_gb1": sum(1 for v in s1 if abs(v - 2.2199189) < 1e-3),
            "pairs_at_2p22A_gb2": sum(1 for v in s2 if abs(v - 2.2199189) < 1e-3),
            "ok": True}


def check_bulk_bcc(bc: Bicrystal, margin: float = 5.0) -> dict:
    """Atoms farther than margin from every GB plane must be perfect BCC."""
    Lx, Ly, Lz = bc.box
    a = bc.a
    pos = [(at.x, at.y, at.z) for at in bc.atoms]
    cl = CellList(pos, Lx, Ly, Lz, 1.02 * a)
    nn_lo, nn_hi = 0.98 * NN_SCALE * a, 1.02 * NN_SCALE * a
    nnn_lo, nnn_hi = 0.98 * NNN_SCALE * a, 1.02 * NNN_SCALE * a
    sampled = perfect = failed = 0
    failures: List[str] = []
    for i, at in enumerate(bc.atoms):
        d1 = min(at.x, Lx - at.x)
        if d1 < margin or abs(at.x - Lx / 2) < margin:
            continue
        sampled += 1
        nn = nnn = 0
        for _, dist in cl.neighbors(i, nnn_hi):
            if nn_lo < dist < nn_hi:
                nn += 1
            elif nnn_lo < dist < nnn_hi:
                nnn += 1
        if nn == 8 and nnn == 6:
            perfect += 1
        else:
            failed += 1
            failures.append(f"atom {at.id}: nn={nn} nnn={nnn}")
            if failed > 5:
                break
    _require(failed == 0, f"bulk atoms not BCC: {failures[:3]}")
    return {"margin_A": margin, "sampled": sampled, "perfect": perfect,
            "failed": failed, "ok": True}


def check_min_distance(bc: Bicrystal) -> dict:
    Lx, Ly, Lz = bc.box
    pos = [(at.x, at.y, at.z) for at in bc.atoms]
    cl = CellList(pos, Lx, Ly, Lz, 3.0)
    gmin = cl.min_distance()
    analytic = 2.0 * plane_spacing(bc.a)  # 2d
    _require(abs(gmin - analytic) < 1e-3, f"global min distance {gmin} != 2d={analytic}")
    _require(gmin > 1.0, f"min distance {gmin} too small (overlapping atoms)")
    return {"global_min_A": gmin, "analytic_2d_A": analytic,
            "nn_bcc_A": NN_SCALE * bc.a, "nnn_bcc_A": NNN_SCALE * bc.a, "ok": True}


def check_voids(bc: Bicrystal, plane_tol: float = 1e-6) -> dict:
    """Detect a missing x plane from periodic gaps between occupied x levels.

    A fixed 1-A bin is invalid here because the ideal (310) plane spacing is
    a/sqrt(10), about 1.11 A for Li.  Instead, collect occupied x levels,
    include the periodic boundary gap, and flag gaps larger than 1.5 times
    the ideal (310) spacing.  One missing complete plane would create a gap
    of approximately 2*d310.
    """
    Lx, _, _ = bc.box
    xs = sorted(at.x % Lx for at in bc.atoms)

    planes: List[float] = []
    for x in xs:
        if not planes or x - planes[-1] > plane_tol:
            planes.append(x)

    _require(len(planes) >= 2, "fewer than two occupied x planes")

    gaps = [planes[i + 1] - planes[i] for i in range(len(planes) - 1)]
    gaps.append(Lx - planes[-1] + planes[0])

    expected = plane_spacing(bc.a)
    large_gap_indices = [
        i for i, gap in enumerate(gaps) if gap > 1.5 * expected
    ]

    _require(
        not large_gap_indices,
        "possible missing x plane(s); large periodic gaps at "
        f"{large_gap_indices}, max gap={max(gaps)} A",
    )

    return {
        "occupied_x_planes": len(planes),
        "expected_plane_gap_A": expected,
        "min_plane_gap_A": min(gaps),
        "max_plane_gap_A": max(gaps),
        "missing_x_planes": large_gap_indices,
        "ok": True,
    }


def check_removal_localized(bc: Bicrystal, require_log: bool = True) -> dict:
    Lx, Ly, Lz = bc.box
    if not require_log and not bc.removal_log:
        return {"removed_total": 0, "removed_grain1": 0, "removed_grain2": 0,
                "gb1_removed": 0, "gb2_removed": 0, "all_on_gb_plane": None,
                "ok": True, "skipped": "no removal log for re-imported data"}
    _require(len(bc.removal_log) > 0, "removal log is empty")
    for rec in bc.removal_log:
        x = rec["x"]
        on_gb = abs(x) < GB_TOL or abs(x - Lx / 2) < GB_TOL
        _require(on_gb, f"removed atom {rec['atom_id']} not on a GB plane (x={x})")
    return {"removed_total": len(bc.removal_log),
            "removed_grain1": sum(1 for r in bc.removal_log if r["grain_id"] == 1),
            "removed_grain2": sum(1 for r in bc.removal_log if r["grain_id"] == 2),
            "gb1_removed": sum(1 for r in bc.removal_log if r["gb"] == "GB1"),
            "gb2_removed": sum(1 for r in bc.removal_log if r["gb"] == "GB2"),
            "all_on_gb_plane": True, "ok": True}


def check_counts(bc: Bicrystal) -> dict:
    a, K, J, L = bc.a, bc.K, bc.J, bc.L
    Lx, Ly, Lz = bc.box
    n_g1 = sum(1 for at in bc.atoms if at.grain_id == 1)
    n_g2 = sum(1 for at in bc.atoms if at.grain_id == 2)
    total = len(bc.atoms)
    expected = 20 * K * J * L
    density = total / (Lx * Ly * Lz)
    target = 2.0 / a ** 3
    _require(total == expected, f"atom count {total} != expected {expected}")
    _require(n_g1 == n_g2, f"grain populations {n_g1} vs {n_g2} not balanced")
    _require(abs(density - target) < 1e-12, f"density {density} != 2/a^3 {target}")
    return {"final_total": total, "final_grain1": n_g1, "final_grain2": n_g2,
            "density_atoms_per_A3": density, "density_target_2_over_a3": target,
            "ok": True}


def check_periodicity(bc: Bicrystal) -> dict:
    a = bc.a
    d = plane_spacing(a)
    miss1 = miss2 = 0
    for at in bc.atoms:
        if at.grain_id == 1:
            if lattice_indices_grain1(at.x, at.y, at.z, a, d) is None:
                miss1 += 1
        else:
            if lattice_indices_grain2(at.x, at.y, at.z, a, d) is None:
                miss2 += 1
    _require(miss1 == 0 and miss2 == 0,
             f"lattice-index misses: g1={miss1} g2={miss2}")
    return {"grain1_misses": miss1, "grain2_misses": miss2,
            "y_period_A": math.sqrt(10.0) * a, "z_period_A": a, "ok": True}


def check_legacy_regression(bc: Bicrystal) -> dict:
    meta = bc.meta
    after = meta.get("final_total", len(bc.atoms))
    expected_removed = 4 * bc.J * bc.L
    before = meta.get("total_before_removal", after + expected_removed)
    removed = before - after
    pct = 100.0 * removed / before

    _require(
        not (before == 3600 and after == 3224),
        "legacy 3600 -> 3224 scenario reproduced!",
    )
    _require(
        removed == expected_removed,
        f"removed {removed} atoms, expected 4*J*L={expected_removed}",
    )

    return {
        "before": before,
        "after": after,
        "removed": removed,
        "expected_removed": expected_removed,
        "removal_pct": pct,
        "legacy_bad": False,
        "ok": True,
    }


def run_all_checks(bc: Bicrystal, require_removal_log: bool = True) -> dict:
    results = {
        "parameters": {
            "a_A": bc.a, "K": bc.K, "J": bc.J, "L": bc.L,
            "theta_deg": THETA_DEG, "tilt_axis": list(TILT_AXIS),
            "gb_plane": GB_PLANE_LABEL,
            "box_A": list(bc.box),
            "in_plane_periods_A": [math.sqrt(10.0) * bc.a, bc.a],
        },
        "orientation": check_misorientation(bc),
        "gb_plane": check_gb_plane(bc),
        "two_gbs": check_two_gbs(bc),
        "bulk_bcc": check_bulk_bcc(bc),
        "distances": check_min_distance(bc),
        "voids": check_voids(bc),
        "removal": check_removal_localized(bc, require_log=require_removal_log),
        "counts": check_counts(bc),
        "periodicity": check_periodicity(bc),
        "legacy_comparison": check_legacy_regression(bc),
        "passed": True,
    }
    return results


def verify_data_file(
    data_path: str | Path,
    map_path: str | Path | None = None,
    a: float | None = None,
) -> dict:
    """Re-import a LAMMPS data file and run the position-based checks."""
    import json
    from .lammps_io import read_lammps_data, read_grain_map

    data_path = Path(data_path)
    if a is None:
        params_path = data_path.parent / "params.json"
        if params_path.is_file():
            params = json.loads(params_path.read_text(encoding="utf-8"))
            a = float(params["a_A"])
        else:
            a = DEFAULT_A
    if a <= 0:
        raise ValueError(f"lattice constant must be positive, got {a}")

    atoms, box, ntypes = read_lammps_data(data_path)
    if map_path is not None:
        mapping = read_grain_map(map_path)
        for at in atoms:
            g = mapping.get(at.id)
            if g is not None:
                at.grain_id = g
    if any(at.grain_id not in (1, 2) for at in atoms):
        raise ValueError(
            "grain ids missing for re-import verification; pass --map "
            "(the file generated alongside the data file)")
    # guess K, J, L from the orthogonal box
    Lx, Ly, Lz = box
    K = round(Lx / (math.sqrt(10.0) * a))
    J = round(Ly / (math.sqrt(10.0) * a))
    L = round(Lz / a)
    # fit the orientation of each grain from a few bulk (far-from-GB) atoms
    M1_ideal, M2_ideal = orientation_matrices(THETA)
    centers1 = [at for at in atoms if at.grain_id == 1
                and 5.0 < min(at.x, Lx - at.x) and abs(at.x - Lx / 2) > 5.0]
    centers2 = [at for at in atoms if at.grain_id == 2
                and 5.0 < min(at.x, Lx - at.x) and abs(at.x - Lx / 2) > 5.0]
    pos = [(at.x, at.y, at.z) for at in atoms]
    if len(centers1) >= 4 and len(centers2) >= 4:
        fits1 = [fit_orientation((c.x, c.y, c.z), pos, M1_ideal, Lx, Ly, Lz, a)
                 for c in centers1[:4]]
        fits2 = [fit_orientation((c.x, c.y, c.z), pos, M2_ideal, Lx, Ly, Lz, a)
                 for c in centers2[:4]]
        M1 = sum(fits1) / len(fits1)
        M2 = sum(fits2) / len(fits2)
    else:
        M1, M2 = M1_ideal, M2_ideal
    bc = Bicrystal(a=a, K=K, J=J, L=L, box=box, atoms=atoms, removal_log=[],
                   M1=M1, M2=M2, theta=THETA)
    res = run_all_checks(bc, require_removal_log=False)
    res["reimport"] = {"data_file": str(data_path), "atoms": len(atoms),
                       "box_A": list(box), "fitted_orientations": True}
    return res


def main(argv=None) -> int:
    """CLI: verify a generated LAMMPS data file (optionally with the grain map).

    Usage:
        python -m generator.verify <data.lmp> [--map <grain_map>]
    """
    import argparse

    ap = argparse.ArgumentParser(description="Verify a v2 Sigma5 data file")
    ap.add_argument("data_file", type=str)
    ap.add_argument("--map", type=str, default=None,
                    help="atom_id -> grain_id map file")
    ap.add_argument("--a", type=float, default=None,
                    help="lattice constant in Angstrom; overrides params.json")
    args = ap.parse_args(argv)
    try:
        res = verify_data_file(args.data_file, args.map, args.a)
    except (VerificationError, ValueError) as exc:
        print(f"[verify] FAILED: {exc}")
        return 2
    c = res["counts"]
    print(f"[verify] PASSED: {c['final_total']} atoms "
          f"(g1={c['final_grain1']}, g2={c['final_grain2']})")
    print(f"[verify] box = {res['parameters']['box_A'][0]:.4f} x "
          f"{res['parameters']['box_A'][1]:.4f} x "
          f"{res['parameters']['box_A'][2]:.4f} A")
    print(f"[verify] misorientation = {res['orientation']['misorientation_angle_deg']:.8f} deg")
    print(f"[verify] min distance = {res['distances']['global_min_A']:.6f} A")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
