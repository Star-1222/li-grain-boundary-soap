"""Command-line driver for the SOAP grain-boundary characterization.

[NEW 2026-08-21]

Single-structure analysis (research item (e)):

    python -m soap_analysis.run_analysis \
        --structure path/to/stage4_final.data \
        [--peratom path/to/stage4_final.peratom] \
        [--params  path/to/build/params.json] \
        [--outdir  soap] [--backend numpy] \
        [--rcut 4.0 --n-max 6 --l-max 6 --sigma 0.5]

produces
    soap/summary.json          all results (also consumed by --compare)
    soap/disorder_peratom.csv  id, x, y, z, disorder, pe, ke per atom
    soap/disorder_profile.csv  x bin centres and binned mean disorder

Study-wide comparison across GBs:

    python -m soap_analysis.run_analysis --compare path/to/gb_study

scans <study>/<gb>/soap/summary.json files and writes
    <study>/comparison.csv     disorder statistics vs misorientation angle
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from .soap_backend import BACKENDS, SoapHyper, available_backends
from .bulk_reference import bulk_soap_reference
from .disorder_profile import (disorder_profile, gb_width_regions,
                               pearson_correlation, per_atom_disorder)
from .io_lammps import peratom_energies, read_peratom_dump, read_structure
from .compare_structures import disorder_vs_theta

# lattice constant of the relaxed bulk BCC Li reference (v2_production
# in.bulk_energy result); overridden by --params / --a.
DEFAULT_A = 3.506186475312


def _load_a(args) -> float:
    if args.a is not None:
        return float(args.a)
    if args.params is not None:
        with open(args.params) as fh:
            p = json.load(fh)
        if "a_A" in p:
            return float(p["a_A"])
        if "a" in p:
            return float(p["a"])
    return DEFAULT_A


def analyze_structure(structure_path: Path, *, hyper: SoapHyper,
                      backend: str, a: float, element: str,
                      peratom_path: Optional[Path],
                      params_path: Optional[Path],
                      outdir: Path, bin_width: float) -> Dict:
    """Full single-structure SOAP characterization.  Returns the summary."""
    structure = read_structure(structure_path, element=element)
    positions = structure["positions"]
    cell = structure["cell"]
    ids = structure["ids"]
    n = structure["natoms"]
    print(f"[soap] structure: {structure_path} ({n} atoms, "
          f"box {cell[0]:.3f} x {cell[1]:.3f} x {cell[2]:.3f} A)")

    # ---- descriptors -----------------------------------------------------
    X = None
    if backend in ("auto",):
        chosen = available_backends()[0]
    else:
        chosen = backend
    print(f"[soap] backend: {chosen} "
          f"(rcut={hyper.rcut} n_max={hyper.n_max} l_max={hyper.l_max} "
          f"sigma={hyper.sigma})")
    from .soap_backend import soap_descriptors
    X = soap_descriptors(positions, cell, [element] * n,
                         hyper=hyper, backend=backend, normalize=True)

    # ---- bulk reference & disorder ---------------------------------------
    bulk_vec = bulk_soap_reference(a, hyper=hyper, backend=backend,
                                   element=element)
    disorder = per_atom_disorder(X, bulk_vec)
    print(f"[soap] disorder: mean={disorder.mean():.6f} "
          f"max={disorder.max():.6f} (bulk reference a={a:.6f} A)")

    # ---- profile along the GB normal (x) ---------------------------------
    prof = disorder_profile(positions[:, 0], disorder,
                            bin_width=bin_width, box_x=float(cell[0]))
    regions = gb_width_regions(prof)
    print(f"[soap] GB profile peak={prof['profile_max']:.6f}; "
          f"{len(regions)} region(s) above half max: "
          + "; ".join(f"[{r['x_lo']:.1f},{r['x_hi']:.1f}] w={r['width_A']:.1f}A"
                      for r in regions))

    # ---- optional per-atom energies --------------------------------------
    pe = ke = None
    energy_stats: Dict = {}
    if peratom_path is not None and Path(peratom_path).is_file():
        frame = read_peratom_dump(peratom_path)
        eids, pe, ke = peratom_energies(frame)
        if len(eids) == n and np.array_equal(eids, ids):
            r_pe = pearson_correlation(disorder, pe)
            energy_stats = {"peratom_file": str(peratom_path),
                            "pearson_r_disorder_vs_pe": r_pe,
                            "mean_pe_eV": float(np.mean(pe)),
                            "max_pe_eV": float(np.max(pe))}
            print(f"[soap] disorder/energy Pearson r = {r_pe:.4f}")
        else:
            print(f"[soap] WARNING: per-atom ids in {peratom_path} do not "
                  f"match the structure; energy correlation skipped")
            pe = ke = None

    # ---- metadata --------------------------------------------------------
    meta: Dict = {}
    if params_path is not None and Path(params_path).is_file():
        try:
            with open(params_path) as fh:
                p = json.load(fh)
            meta = {"theta_deg": p.get("theta_deg"),
                    "sigma": p.get("sigma"),
                    "plane": p.get("plane", p.get("gb_plane"))}
        except (OSError, json.JSONDecodeError):
            pass

    if meta.get("plane"):
        label = f"Sigma{meta.get('sigma')}_{meta['plane']}"
    else:
        label = structure_path.stem

    summary = {
        "label": label,
        "structure": str(structure_path),
        "peratom": str(peratom_path) if peratom_path else None,
        "params": str(params_path) if params_path else None,
        "backend": chosen,
        "hyperparameters": {"rcut": hyper.rcut, "n_max": hyper.n_max,
                            "l_max": hyper.l_max, "sigma": hyper.sigma},
        "bulk_reference_a_A": a,
        "natoms": n,
        "cell_A": [float(c) for c in cell],
        "theta_deg": meta.get("theta_deg"),
        "sigma": meta.get("sigma"),
        "plane": meta.get("plane"),
        "disorder": {"mean": float(disorder.mean()),
                     "max": float(disorder.max()),
                     "std": float(disorder.std())},
        "gb": {"bin_width_A": bin_width,
               "profile_peak": prof["profile_max"],
               "regions": regions},
        "energy": energy_stats,
    }

    # ---- outputs ----------------------------------------------------------
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    with open(outdir / "disorder_peratom.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        header = ["id", "x", "y", "z", "disorder"]
        if pe is not None:
            header += ["pe_eV", "ke_eV"]
        w.writerow(header)
        for i in range(n):
            row = [int(ids[i]), f"{positions[i, 0]:.6f}",
                   f"{positions[i, 1]:.6f}", f"{positions[i, 2]:.6f}",
                   f"{disorder[i]:.8f}"]
            if pe is not None:
                row += [f"{pe[i]:.8f}", f"{ke[i]:.8f}"]
            w.writerow(row)
    with open(outdir / "disorder_profile.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["x_center_A", "mean_disorder", "std_disorder",
                    "n_atoms", "mean_pe_eV"])
        mpe = None
        if pe is not None:
            x = positions[:, 0]
            n_bins = len(prof["bin_centers"])
            idx = np.clip(np.floor(np.mod(x, cell[0]) / bin_width).astype(int),
                          0, n_bins - 1)
            mpe = np.full(n_bins, np.nan)
            s = np.bincount(idx, weights=pe, minlength=n_bins)
            c = np.bincount(idx, minlength=n_bins)
            occ = c > 0
            mpe[occ] = s[occ] / c[occ]
        for i, xc in enumerate(prof["bin_centers"]):
            w.writerow([f"{xc:.3f}", f"{prof['mean_disorder'][i]:.8f}",
                        f"{prof['std_disorder'][i]:.8f}",
                        int(prof["counts"][i]),
                        "" if mpe is None else f"{mpe[i]:.8f}"])
    print(f"[soap] wrote {outdir}/summary.json, disorder_peratom.csv, "
          f"disorder_profile.csv")
    return summary


def compare_study(study_root: Path) -> None:
    """Aggregate per-GB summaries into a comparison CSV."""
    summaries = []
    for sj in sorted(study_root.glob("*/soap/summary.json")):
        with open(sj) as fh:
            summaries.append(json.load(fh))
    if not summaries:
        print(f"[compare] no */soap/summary.json found under {study_root}")
        return
    rows = disorder_vs_theta(summaries)
    out = study_root / "comparison.csv"
    fields = ["label", "plane", "sigma", "theta_deg", "natoms",
              "disorder_mean", "disorder_max", "disorder_gb_peak",
              "gb_width_A", "energy_disorder_r"]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in fields})
    print(f"[compare] wrote {out} ({len(rows)} grain boundaries)")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="SOAP characterization of equilibrated Li grain boundaries")
    ap.add_argument("--structure", type=str, default=None,
                    help="LAMMPS data file of the (relaxed) GB structure")
    ap.add_argument("--peratom", type=str, default=None,
                    help="per-atom dump with c_peatom/c_keatom columns "
                         "(stage4_final.peratom)")
    ap.add_argument("--params", type=str, default=None,
                    help="generator params.json (a, theta, sigma, plane)")
    ap.add_argument("--outdir", type=str, default="soap")
    ap.add_argument("--backend", type=str, default="auto",
                    choices=list(BACKENDS))
    ap.add_argument("--rcut", type=float, default=4.0)
    ap.add_argument("--n-max", type=int, default=6)
    ap.add_argument("--l-max", type=int, default=6)
    ap.add_argument("--sigma", type=float, default=0.5)
    ap.add_argument("--a", type=float, default=None,
                    help="bulk lattice constant [A] (default: from --params "
                         f"or {DEFAULT_A})")
    ap.add_argument("--element", type=str, default="Li")
    ap.add_argument("--bin-width", type=float, default=1.0)
    ap.add_argument("--compare", type=str, default=None,
                    help="study root with */soap/summary.json files; writes "
                         "comparison.csv")
    args = ap.parse_args(argv)

    if args.compare is not None:
        compare_study(Path(args.compare))
        return 0
    if args.structure is None:
        ap.error("--structure is required (or use --compare STUDY_ROOT)")

    hyper = SoapHyper(rcut=args.rcut, n_max=args.n_max, l_max=args.l_max,
                      sigma=args.sigma)
    a = _load_a(args)
    analyze_structure(Path(args.structure), hyper=hyper, backend=args.backend,
                      a=a, element=args.element,
                      peratom_path=Path(args.peratom) if args.peratom else None,
                      params_path=Path(args.params) if args.params else None,
                      outdir=Path(args.outdir), bin_width=args.bin_width)
    return 0


if __name__ == "__main__":
    sys.exit(main())
