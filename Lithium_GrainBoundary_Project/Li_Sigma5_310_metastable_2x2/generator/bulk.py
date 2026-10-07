"""Bulk BCC Li reference structure for the same SNAP potential.

Builds a perfect BCC single crystal (conventional cubic cell repeated
n x n x n) with the same lattice constant a as the bicrystal, and
exports it as LAMMPS data / extended XYZ / CIF.  The matching LAMMPS
input in.bulk_energy.lmp computes the cohesive reference energy.

Usage:
    python -m generator.bulk [--a 3.51] [--n 4] [--outdir outputs/bulk]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from .builder import Atom
from .config import DEFAULT_A, ELEMENT, MASS_LI
from .exporters import write_cif, write_extxyz, write_lammps_data


def build_bulk(a: float = DEFAULT_A, n: int = 4):
    """Return (atoms, box) of a perfect BCC crystal n^3 conventional cells."""
    box = (n * a, n * a, n * a)
    basis = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)]
    atoms = []
    i = 1
    for ix in range(n):
        for iy in range(n):
            for iz in range(n):
                for (fx, fy, fz) in basis:
                    atoms.append(Atom(id=i, type=1,
                                     x=(ix + fx) * a, y=(iy + fy) * a,
                                     z=(iz + fz) * a, grain_id=1))
                    i += 1
    return atoms, box


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="BCC Li bulk reference cell")
    ap.add_argument("--a", type=float, default=DEFAULT_A)
    ap.add_argument("--n", type=int, default=4)
    ap.add_argument("--outdir", type=str, default="outputs/bulk")
    args = ap.parse_args(argv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    atoms, box = build_bulk(a=args.a, n=args.n)

    class _Bulk:
        pass
    b = _Bulk()
    b.atoms = atoms
    b.box = box
    b.removal_log = []
    b.meta = {"a_A": args.a, "n": args.n, "atoms": len(atoms)}
    b.M1 = np.eye(3)
    b.M2 = np.eye(3)
    b.theta = 0.0

    write_lammps_data(b, outdir / "bulk_Li.data")
    write_extxyz(b, outdir / "bulk_Li.xyz")
    write_cif(b, outdir / "bulk_Li.cif")
    print(f"[ok] bulk BCC Li: {len(atoms)} atoms, box {box[0]:.4f} A^3")
    print("      expected density 2/a^3 =", 2.0 / args.a ** 3)
    print("      run lmp -in lammps/in.bulk_energy.lmp for the reference energy")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
