#!/usr/bin/env python3
"""Structural analysis for the Li Sigma5(310)[001] bicrystal.

Four structural descriptors are computed for every snapshot.  They answer
different questions, and the disagreement between two of them is itself the
result -- see the notes below.

  A.  nn8+nnn6 (absolute)  -- the project's own criterion, identical to
      generator/verify.py::check_bulk_bcc : exactly 8 neighbours inside
      [1-tol,1+tol]*(sqrt3/2)*a and exactly 6 inside [1-tol,1+tol]*a, with the
      a used at build time.  Validated against perfect lattices.  Because the
      bands are absolute, this responds to *elastic* strain as well as to
      disorder.  On a thermally populated MD snapshot it is far too strict:

  A'. nn8+nnn6 (rescaled) -- identical rule, but the bands are rescaled to the
      lattice constant the box itself implies, a* = (2V/N)^(1/3).  Removes the
      uniform-strain part of A.

  C.  coordination within the first g(r) minimum == 8 -- a topological rule.
      r_cut is taken at the first minimum of the structure's own radial
      distribution, so a *uniform* strain cannot move atoms across it (the
      1NN/2NN shells sit at 3.03 and 3.50 A, i.e. 15 % apart).  Insensitive to
      elastic strain, sensitive to genuine coordination defects.

  B.  Steinhardt q4/q6 from the 14 nearest neighbours (8x1NN + 6x2NN).  Purely
      angular -> invariant under any affine strain.  Implementation validated
      against perfect lattices: FCC -> (0.1909, 0.5745), BCC -> (0.0364,
      0.5107), both matching the published values.

Usage:
  python3 analyze_structure.py FILE.data [--a 3.506] [--tol 0.02]
                               [--rcut 3.35] [--name tag] [--outdir DIR]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import KDTree

try:                                   # scipy >= 1.15
    from scipy.special import sph_harm_y as _sphharm

    def _Y(l, theta_colat, phi_azim):
        # sph_harm_y(n, m, theta=polar, phi=azimuthal)
        return _sphharm(l, np.arange(-l, l + 1)[None, :],
                        theta_colat.reshape(-1)[:, None],
                        phi_azim.reshape(-1)[:, None])
except ImportError:                    # fallback for older scipy
    from scipy.special import sph_harm as _sphharm

    def _Y(l, theta_colat, phi_azim):
        # deprecated sph_harm(m, n, theta=azimuthal, phi=polar)
        return _sphharm(np.arange(-l, l + 1)[None, :], l,
                        phi_azim.reshape(-1)[:, None],
                        theta_colat.reshape(-1)[:, None])


NN_SCALE = math.sqrt(3.0) / 2.0
NNN_SCALE = 1.0
Q_REF = {"fcc": (0.1909, 0.5745), "bcc": (0.0364, 0.5107)}


# ------------------------------------------------------------------ IO ----
def read_lammps_data(path: str):
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    bounds = {"x": None, "y": None, "z": None}
    natoms = None
    start = None
    for ln, raw in enumerate(lines):
        s = raw.strip()
        head = s.split()
        if len(head) >= 2 and head[-1] == "atoms":
            try:
                natoms = int(head[0])
            except ValueError:
                pass
        if len(head) >= 4 and head[-2] in ("xlo", "ylo", "zlo"):
            try:
                bounds[head[-2][0]] = (float(head[0]), float(head[1]))
            except ValueError:
                pass
        if s.startswith("Atoms"):
            start = ln + 1
            break
    if any(v is None for v in bounds.values()):
        raise SystemExit(f"[error] cannot parse box bounds from {path}")
    if start is None:
        raise SystemExit(f"[error] no Atoms section in {path}")

    L = np.array([bounds[k][1] - bounds[k][0] for k in "xyz"], dtype=float)
    origin = np.array([bounds[k][0] for k in "xyz"])

    rows = []
    for raw in lines[start:]:
        p = raw.strip().split()
        if len(p) < 5:
            if rows:
                break
            continue
        try:
            rows.append((int(p[0]), float(p[2]), float(p[3]), float(p[4])))
        except ValueError:
            if rows:
                break
            continue
    arr = np.array(rows, dtype=float)
    if natoms is not None and len(arr) != natoms:
        print(f"[warn] parsed {len(arr)} atoms, header says {natoms}", file=sys.stderr)
    arr = arr[np.argsort(arr[:, 0])]
    pos = np.mod(arr[:, 1:4] - origin, L)
    return np.ascontiguousarray(pos), L


def min_image(disp, L):
    return disp - L * np.round(disp / L)


# ------------------------------------------- order parameter A -----------
def criterion_A(pos, L, tree, a, tol):
    nn_lo, nn_hi = (1 - tol) * NN_SCALE * a, (1 + tol) * NN_SCALE * a
    nnn_lo, nnn_hi = (1 - tol) * NNN_SCALE * a, (1 + tol) * NNN_SCALE * a
    pairs = tree.query_pairs(nnn_hi, output_type="ndarray")
    n = len(pos)
    n1 = np.zeros(n)
    n2 = np.zeros(n)
    if len(pairs):
        ii, jj = pairs[:, 0], pairs[:, 1]
        dist = np.linalg.norm(min_image(pos[ii] - pos[jj], L), axis=1)
        n1 = np.bincount(ii[(dist >= nn_lo) & (dist < nn_hi)], minlength=n).astype(float)
        n2 = np.bincount(ii[(dist >= nnn_lo) & (dist < nnn_hi)], minlength=n).astype(float)
    return ((n1 == 8) & (n2 == 6)), n1, n2


# ------------------------------------------- order parameter C -----------
def coordination(pos, L, tree, rcut):
    pairs = tree.query_pairs(rcut, output_type="ndarray")
    n = len(pos)
    cn = np.zeros(n)
    if len(pairs):
        ii = pairs[:, 0]
        cn = np.bincount(ii, minlength=n).astype(float)
        # each pair found once by query_pairs -> both partners need it
        jj = pairs[:, 1]
        cn += np.bincount(jj, minlength=n).astype(float)
    return cn


# ------------------------------------------- order parameter B -----------
def steinhardt(pos, L, tree, k=14, lvals=(4, 6)):
    _, idx = tree.query(pos, k=k + 1)
    nb = idx[:, 1:]
    disp = min_image(pos[nb] - pos[:, None, :], L)
    dist = np.linalg.norm(disp, axis=2)
    unit = disp / dist[:, :, None]
    theta_colat = np.arccos(np.clip(unit[:, :, 2], -1.0, 1.0))
    phi_azim = np.arctan2(unit[:, :, 1], unit[:, :, 0])
    n = len(pos)
    out = {}
    for l in lvals:
        Y = _Y(l, theta_colat, phi_azim).reshape(n, k, 2 * l + 1)
        qbar = Y.mean(axis=1)
        out[l] = np.sqrt(4.0 * np.pi / (2 * l + 1) *
                         np.sum(np.abs(qbar) ** 2, axis=1))
    return out


# -------------------------------------------------- radial g(r) ----------
def rdf(pos, L, tree, rmax=7.0, nbins=350):
    pairs = tree.query_pairs(rmax, output_type="ndarray")
    ii, jj = pairs[:, 0], pairs[:, 1]
    dist = np.linalg.norm(min_image(pos[ii] - pos[jj], L), axis=1)
    hist, edges = np.histogram(dist, bins=nbins, range=(0.0, rmax))
    r = 0.5 * (edges[1:] + edges[:-1])
    shell = 4.0 / 3.0 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    rho = len(pos) / np.prod(L)
    return r, hist / (shell * rho * len(pos)) * 2.0


def first_minimum(r, g, rlo=2.5, rhi=4.2, smooth=7):
    """Locate the first *local* minimum of g(r) after the first peak.

    A plain argmin over the whole window picks up the much deeper minimum that
    follows the second shell (~3.9 A for BCC Li); that would put the 2NN shell
    inside the coordination sphere and is not what we want.  The curve is
    lightly smoothed first so that bin noise cannot stop the walk early.
    """
    sel = (r >= rlo) & (r <= rhi)
    rs, gs = r[sel], g[sel]
    kernel = np.ones(smooth) / smooth
    gs_s = np.convolve(gs, kernel, mode="same")
    ipk = int(np.argmax(gs_s))
    i = ipk
    while i + 1 < len(gs_s) and gs_s[i + 1] <= gs_s[i]:
        i += 1
    return float(rs[ipk]), float(rs[i])


# --------------------------------------- GB-normal profile & width -------
def profile_x(pos, L, mask, nbins=120):
    edges = np.linspace(0.0, L[0], nbins + 1)
    idx = np.clip(np.digitize(pos[:, 0], edges) - 1, 0, nbins - 1)
    cnt = np.bincount(idx, minlength=nbins).astype(float)
    hit = np.bincount(idx[mask.astype(bool)], minlength=nbins).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        q = np.where(cnt > 0, hit / np.maximum(cnt, 1), np.nan)
    return 0.5 * (edges[1:] + edges[:-1]), q, cnt


def gb_width(ctr, q, nbins):
    finite = np.isfinite(q)
    qf = np.interp(ctr, ctr[finite], q[finite])
    q_bulk = 0.5 * (qf[int(0.25 * nbins)] + qf[int(0.75 * nbins)])
    thr = 0.5 * (q_bulk + qf.min())
    bw = ctr[1] - ctr[0]
    return {"w_GB_A": float((qf < thr).sum() * bw),
            "bulk_plateau_Q": float(q_bulk),
            "min_Q": float(qf.min()),
            "threshold_Q": float(thr)}


# ------------------------------------------------------------- main ------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("datafile")
    ap.add_argument("--a", type=float, default=3.506)
    ap.add_argument("--tol", type=float, default=0.02)
    ap.add_argument("--rcut", type=float, default=None,
                    help="coordination cutoff; default = first g(r) minimum")
    ap.add_argument("--nbins", type=int, default=120)
    ap.add_argument("--name", default=None)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    name = args.name or Path(args.datafile).stem
    pos, L = read_lammps_data(args.datafile)
    n = len(pos)
    print(f"[info] {name}: N={n} box={L[0]:.5f} x {L[1]:.5f} x {L[2]:.5f} A "
          f"V/atom={np.prod(L)/n:.4f} A^3")

    tree = KDTree(pos, boxsize=L)
    a_star = (2.0 * np.prod(L) / n) ** (1.0 / 3.0)

    r, g = rdf(pos, L, tree)
    rp, rm = first_minimum(r, g)
    rcut = args.rcut if args.rcut else rm
    print(f"[rdf] 1st peak {rp:.3f} A, 1st minimum {rm:.3f} A -> rcut {rcut:.3f} A")

    isA, n1, n2 = criterion_A(pos, L, tree, args.a, args.tol)
    isA2, _, _ = criterion_A(pos, L, tree, a_star, args.tol)
    cn = coordination(pos, L, tree, rcut)
    isC = cn == 8
    qs = steinhardt(pos, L, tree, k=14)
    q4, q6 = qs[4], qs[6]
    # BCC-like: close to the perfect-BCC point in (q4,q6) space
    d_ref = np.sqrt((q4 - Q_REF["bcc"][0]) ** 2 + (q6 - Q_REF["bcc"][1]) ** 2)
    isB = d_ref < 0.12

    cx, qxC, _ = profile_x(pos, L, isC, args.nbins)
    _, qxB, _ = profile_x(pos, L, isB, args.nbins)
    _, qxA, _ = profile_x(pos, L, isA2, args.nbins)
    wC = gb_width(cx, qxC, args.nbins)
    wB = gb_width(cx, qxB, args.nbins)

    print(f"[A ] nn8+nnn6 a={args.a:.4f}      frac={isA.mean():.4f} "
          f"({int(isA.sum())}/{n})")
    print(f"[A'] nn8+nnn6 a*={a_star:.4f}     frac={isA2.mean():.4f} "
          f"({int(isA2.sum())}/{n})")
    print(f"[C ] CN(r<{rcut:.2f})==8          frac={isC.mean():.4f} "
          f"({int(isC.sum())}/{n})  <CN>={cn.mean():.3f}")
    print(f"[B ] near perfect-BCC in (q4,q6) frac={isB.mean():.4f} "
          f"({int(isB.sum())}/{n})  <q4>={q4.mean():.4f} <q6>={q6.mean():.4f}")
    print(f"[wGB] coordination: {wC['w_GB_A']:.2f} A   "
          f"plateau={wC['bulk_plateau_Q']:.3f} min={wC['min_Q']:.3f}")

    res = {
        "name": name, "file": str(args.datafile), "n_atoms": int(n),
        "box_A": [float(v) for v in L],
        "volume_per_atom_A3": float(np.prod(L) / n),
        "a_implied_if_bcc_A": float(a_star),
        "reference_a_A": float(args.a),
        "rdf_first_peak_A": rp, "rdf_first_minimum_A": rm, "coordination_rcut_A": rcut,
        "A_nn8nnn6_absolute": {"fraction": float(isA.mean()), "count": int(isA.sum())},
        "Aprime_nn8nnn6_rescaled": {"fraction": float(isA2.mean()),
                                    "count": int(isA2.sum())},
        "C_coordination8": {"fraction": float(isC.mean()), "count": int(isC.sum()),
                            "mean_CN": float(cn.mean())},
        "B_steinhardt_bcc": {"fraction": float(isB.mean()), "count": int(isB.sum()),
                             "q4_mean": float(q4.mean()), "q6_mean": float(q6.mean()),
                             "q4_std": float(q4.std()), "q6_std": float(q6.std())},
        "w_GB_A": {"coordination": wC["w_GB_A"], "steinhardt": wB["w_GB_A"]},
        "profile_coordination": {k: v for k, v in wC.items() if k != "w_GB_A"},
        "profile_steinhardt": {k: v for k, v in wB.items() if k != "w_GB_A"},
    }
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / f"{name}.json").write_text(
        json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
    np.savez(outdir / f"{name}_profiles.npz",
             x=cx, qA=qxA, qAprime=qxA, qC=qxC, qB=qxB, r=r, g=g,
             q4=q4, q6=q6, cn=cn)
    print(f"[out] {outdir / (name + '.json')}")


if __name__ == "__main__":
    main()
