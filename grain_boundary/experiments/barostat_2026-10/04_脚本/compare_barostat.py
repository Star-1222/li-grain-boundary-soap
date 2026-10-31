#!/usr/bin/env python3
"""Parse the barostat-comparison LAMMPS logs into a clean table.

Reads logs/baro_<tag>.lammps for every tag, re-assembles the thermo tables
(LAMMPS splits long tables into repeated header blocks), computes averages
over the second half of each run, and writes

    <outdir>/thermo_<tag>.csv     tidy per-step data
    <outdir>/summary.json         averaged quantities per branch
    <outdir>/summary.md           the same as a markdown table
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np


def parse_thermo(path: Path):
    """Return (steps, dict_of_arrays) from a LAMMPS log file."""
    rows = []
    fields = None
    last_fields = None
    with Path(path).open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if re.match(r"^\s*Step\s+\w", line):
                fields = line.split()
                last_fields = fields
                continue
            if fields is None:
                continue
            parts = line.split()
            if len(parts) != len(fields):
                fields = None          # left the table block
                continue
            try:
                vals = [float(p) for p in parts]
            except ValueError:
                fields = None
                continue
            rows.append(vals)

    if not rows or last_fields is None:
        return None
    arr = np.array(rows)
    data = {f: arr[:, i] for i, f in enumerate(last_fields)}
    return data


def stats(vals, n_last):
    tail = vals[-n_last:] if n_last > 0 else vals
    return {"mean": float(np.mean(tail)), "std": float(np.std(tail, ddof=1)),
            "last": float(vals[-1]), "first": float(vals[0])}


KEYS = ["Temp", "Press", "Pxx", "Pyy", "Pzz", "v_pMean", "v_devP",
        "Lx", "Ly", "Lz", "Volume", "v_eXX", "v_eYY", "v_eZZ",
        "PotEng", "v_peAtom", "v_fBCC"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--logs", default="logs")
    ap.add_argument("--tags", default="xonly,aniso")
    ap.add_argument("--outdir", default="exp_barostat/analysis")
    ap.add_argument("--frac", type=float, default=0.5,
                    help="fraction of the tail used for averaging")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    summary = {}

    for tag in args.tags.split(","):
        path = Path(args.logs) / f"baro_{tag}.lammps"
        if not path.is_file():
            print(f"[warn] missing {path}")
            continue
        data = parse_thermo(path)
        if data is None:
            print(f"[warn] no thermo data in {path}")
            continue
        n = len(data["Step"])
        n_last = max(2, int(n * args.frac))
        summary[tag] = {
            "steps": n,
            "averaged_over_last_steps": int(
                (data["Step"][-1] - data["Step"][-n_last])),
            **{k: stats(data[k], n_last) for k in KEYS if k in data},
        }
        head = ["Step"] + [k for k in KEYS if k in data]
        cols = [data[k] for k in head]
        np.savetxt(outdir / f"thermo_{tag}.csv", np.column_stack(cols),
                   delimiter=",", header=",".join(head), comments="")
        np.savez(outdir / f"thermo_{tag}.npz", **data)
        print(f"[ok] {tag}: {n} rows, averaged over last "
              f"{summary[tag]['averaged_over_last_steps']:.0f} steps")

    (outdir / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    # ---- markdown table -------------------------------------------------
    tags = list(summary.keys())
    if tags:
        lines = []
        for k in KEYS:
            if k not in summary[tags[0]]:
                continue
            row = [f"`{k}`"]
            for t in tags:
                s = summary[t][k]
                row.append(f"{s['mean']:.1f} ± {s['std']:.1f}")
            lines.append("| " + " | ".join(row) + " |")
        header = "| quantity | " + " | ".join(tags) + " |"
        md = (header + "\n|" + "---|" * (len(tags) + 1) + "\n" + "\n".join(lines))
        (outdir / "summary.md").write_text(md + "\n", encoding="utf-8")
        print("\n" + md)


if __name__ == "__main__":
    main()
