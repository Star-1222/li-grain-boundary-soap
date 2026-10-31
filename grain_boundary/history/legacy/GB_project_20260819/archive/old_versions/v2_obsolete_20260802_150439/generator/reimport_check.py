"""Re-import the exported structure files with ASE and cross-check them.

Reads the LAMMPS data, extended XYZ and CIF produced by generate.py back
with ASE and verifies that (a) the atom count matches, (b) the box matches,
(c) the cell is periodic in all three directions.  This validates the
exporters independently of the builder.

Usage:
    python -m generator.reimport_check [--dir outputs] [--prefix Sigma5_310]

Requires ASE (pip install ase).  Exits non-zero if ASE is missing or a
check fails.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

from .config import DEFAULT_A, KIM_MODEL


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="ASE re-import cross-check")
    ap.add_argument("--dir", type=str, default="outputs")
    ap.add_argument("--prefix", type=str, default="Sigma5_310")
    args = ap.parse_args(argv)

    try:
        from ase.io import read
    except ImportError:
        print("[reimport] ERROR: ASE is not installed.  Run: pip install ase")
        return 2

    d = Path(args.dir)
    data = d / f"{args.prefix}.data"
    xyz = d / f"{args.prefix}.xyz"
    cif = d / f"{args.prefix}.cif"
    for f in (data, xyz, cif):
        if not f.exists():
            print(f"[reimport] ERROR: missing {f}")
            return 2

    # reference values from the generated params.json
    params = d / "params.json"
    if params.exists():
        import json
        meta = json.loads(params.read_text(encoding="utf-8"))
        n_ref = meta["final_total"]
        box_ref = [meta.get("box_A") or [0, 0, 0]]
    else:
        # fallback: parse the data header
        from .lammps_io import read_lammps_data
        atoms0, box0, _ = read_lammps_data(data)
        n_ref, box_ref = len(atoms0), box0

    ok = True
    results = {}
    for name, f, fmt in (("lammps-data", data, None),
                         ("extxyz", xyz, "extxyz"),
                         ("cif", cif, "cif")):
        try:
            ats = read(f, format=fmt) if fmt else read(f, format="lammps-data")
        except Exception as exc:  # noqa: BLE001
            print(f"[reimport] FAIL {name}: {exc}")
            ok = False
            continue
        n = len(ats)
        cell = ats.get_cell()
        lengths = cell.lengths()
        pbc = tuple(bool(p) for p in ats.pbc)
        results[name] = {"atoms": n, "lengths": lengths, "pbc": pbc}
        if n != n_ref:
            print(f"[reimport] FAIL {name}: {n} atoms != expected {n_ref}")
            ok = False
        # pbc is physically encoded in extxyz/cif; a LAMMPS data file carries
        # no boundary info (the in.* scripts declare boundary p p p), so for
        # lammps-data we check the cell and that all atoms lie inside it.
        if name != "lammps-data" and pbc != (True, True, True):
            print(f"[reimport] FAIL {name}: pbc={pbc} != (T,T,T)")
            ok = False
        if name == "lammps-data":
            pos = ats.get_positions()
            outside = sum(1 for pp in pos
                          if any(c < -1e-6 or c > length + 1e-6
                                 for c, length in zip(pp, lengths)))
            if outside:
                print(f"[reimport] FAIL {name}: {outside} atoms outside the cell")
                ok = False
        if box_ref:
            bx = box_ref[0] if isinstance(box_ref[0], (list, tuple)) else box_ref
            for i, lb in enumerate(bx):
                if abs(lengths[i] - float(lb)) > 1e-6:
                    print(f"[reimport] FAIL {name}: length[{i}] {lengths[i]} != {lb}")
                    ok = False
        print(f"[reimport] ok {name}: {n} atoms, cell={lengths[0]:.4f} {lengths[1]:.4f} {lengths[2]:.4f}, pbc={pbc}")

    # grain_id consistency in extxyz vs map
    mapf = d / f"{args.prefix}.grain_map"
    if mapf.exists():
        from .lammps_io import read_grain_map
        mapping = read_grain_map(mapf)
        if len(mapping) != n_ref:
            print(f"[reimport] FAIL grain map: {len(mapping)} entries != {n_ref} atoms")
            ok = False
        else:
            print(f"[reimport] ok grain map: {len(mapping)} entries")

    print("[reimport] ALL OK" if ok else "[reimport] FAILURES DETECTED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
