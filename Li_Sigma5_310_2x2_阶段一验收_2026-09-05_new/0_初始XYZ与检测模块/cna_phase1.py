# -*- coding: utf-8 -*-
"""CNA-lite v2 -- fixes reviewer P1 findings.

Fixes vs v1:
  * GB regions are now SYMMETRIC and equal-width for both equivalent GBs.
    v1 bug: gb2 = |xc - Lx/2| < 8 on the folded coordinate xc in (-Lx/2, Lx/2]
    only captured ONE side of GB2 (8 A, ~1920 atoms) while gb1 = |xc| < 8
    captured both sides (16 A, ~3600 atoms). v2 uses true periodic distance
    to each GB plane: d = |mod(x, L) - x_GB|, capped at L/2.
  * Adds Q_BCC(x) folded profile vs distance-to-GB (1 A bins, both GBs folded)
    and estimates the disordered GB width w_GB (half-depth criterion).

Files analyzed (positional args, see bottom): initial / mincheck / quenched / final.
"""
import numpy as np
from scipy.spatial import cKDTree
import sys, io, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

FILES = [
    ("initial   (0 K built GB reference)",     "initial.data"),
    ("p1_min0   (0K all-dim zero-stress relaxed)", "p1_min0_sel.data"),
    ("casecheck (363K anneal 25ps -> 0K)",        "p1_casecheck_0K.data"),
    ("p1_final   (10K hold -> 0K force-minimized)",       "p1_held_0K.data"),
]
FILES = [(l, p) for l, p in FILES if os.path.exists(p)]
W_GB = 8.0   # half-width of GB sampling region [A]
PROF_BIN = 1.0

def read_data(path):
    with open(path, "r") as f:
        lines = f.readlines()
    box, atoms, in_atoms = {}, [], False
    for ln in lines:
        s = ln.strip()
        if len(s.split()) == 4 and s.split()[2] in ("xlo", "ylo", "zlo"):
            p = s.split()
            box[p[2]] = (float(p[0]), float(p[1]))
        if s.startswith("Atoms"):
            in_atoms = True
            continue
        if in_atoms:
            if s == "" and atoms:
                break
            if s.startswith(("Velocities", "Masses", "Bonds")) or (s and not s[0].isdigit()):
                if atoms:
                    break
                continue
            if s:
                p = s.split()
                atoms.append((int(p[0]), float(p[2]), float(p[3]), float(p[4])))
    a = np.array([(x[1], x[2], x[3]) for x in atoms])
    L = np.array([box["xlo"][1] - box["xlo"][0],
                  box["ylo"][1] - box["ylo"][0],
                  box["zlo"][1] - box["zlo"][0]])
    return a, L, box

def bcc_mask(pos, L, rcut1, rcut2):
    tree = cKDTree(pos, boxsize=L)
    n1 = np.array([len(x) for x in tree.query_ball_point(pos, rcut1)]) - 1
    n2 = np.array([len(x) for x in tree.query_ball_point(pos, rcut2)]) - 1
    return (n1 == 8) & (n2 == 14)

# calibration on initial
a0, L0, _ = read_data("initial.data")
tree = cKDTree(a0, boxsize=L0)
d, _ = tree.query(a0, k=2)
r_nn = np.median(d[:, 1])
a_lat = r_nn * 2.0 / np.sqrt(3.0)
rcut1, rcut2 = r_nn * 1.08, a_lat * 1.10
print("=" * 74)
print("CNA-lite v2 report -- symmetric GB sampling (reviewer P1 fix)")
print("=" * 74)
print(f"calibration : initial.data  r_nn = {r_nn:.4f} A  a_BCC = {a_lat:.4f} A")
print(f"criteria    : CN(r < {rcut1:.2f}) == 8  AND  CN(r < {rcut2:.2f}) == 14")
print(f"GB sampling : periodic distance to GB plane, half-width {W_GB} A, both GBs")
print()

# profiles keyed by label for later comparison
profiles = {}

for label, path in FILES:
    pos, L, box = read_data(path)
    lo = np.array([box["xlo"][0], box["ylo"][0], box["zlo"][0]])
    pos = np.mod(pos - lo, L)
    is_bcc = bcc_mask(pos, L, rcut1, rcut2)
    x = pos[:, 0]
    lx = L[0]
    # true periodic distance to each GB plane (GB1 at x~0, GB2 at x~Lx/2)
    u = np.mod(x, lx)                       # [0, lx)
    d1 = np.minimum(u, lx - u)              # distance to GB1
    du2 = np.abs(u - lx / 2.0)
    d2 = np.minimum(du2, lx - du2)          # distance to GB2
    gb1, gb2 = d1 < W_GB, d2 < W_GB
    interior = ~(gb1 | gb2)

    # folded profile vs distance-to-nearest-GB
    dg = np.minimum(d1, d2)                 # distance to nearest GB, [0, lx/2]
    edges = np.arange(0, lx / 2 + PROF_BIN, PROF_BIN)
    idx = np.digitize(dg, edges) - 1
    nb = np.maximum(idx.max() + 1, 1)
    cnt = np.bincount(idx, minlength=nb)
    bcc = np.bincount(idx, weights=is_bcc.astype(float), minlength=nb)
    prof = np.where(cnt > 0, bcc / np.maximum(cnt, 1), np.nan)
    centers = edges[:nb] + PROF_BIN / 2
    profiles[label] = (centers, prof, cnt)

    print(f"--- {label} ({path}) ---")
    print(f"    N = {len(pos)}   total BCC = {is_bcc.mean()*100:6.2f}%")
    print(f"    GB region  x ~ 0    (|d| < {W_GB:.0f} A): BCC = {is_bcc[gb1].mean()*100:5.1f}%  [{gb1.sum()} atoms]")
    print(f"    GB region  x ~ Lx/2 (|d| < {W_GB:.0f} A): BCC = {is_bcc[gb2].mean()*100:5.1f}%  [{gb2.sum()} atoms]")
    print(f"    grain interior              : BCC = {is_bcc[interior].mean()*100:5.1f}%  [{interior.sum()} atoms]")
    # w_GB: half-depth between GB minimum and interior plateau
    plateau = np.nanmean(prof[centers > 12.0])
    gmin = np.nanmin(prof[centers < 6.0])
    half = 0.5 * (plateau + gmin)
    below = centers[(prof < half) & ~np.isnan(prof)]
    w_gb = below.max() if len(below) else 0.0
    print(f"    Q_BCC plateau (d>12 A) = {plateau*100:.1f}%, GB core min = {gmin*100:.1f}%, "
          f"half-depth threshold = {half*100:.1f}%")
    print(f"    w_GB (half-depth, last bin below threshold) = {w_gb:.1f} A")
    print()

print("=" * 74)
print("Q_BCC profile vs distance-to-GB (folded, 1 A bins) -- BCC fraction [%]")
print("=" * 74)
print(f"{'d_to_GB[A]':>10} " + " ".join(f"{l.split()[0]:>10}" for l, _ in FILES))
for i in range(len(centers)):
    row = " ".join(f"{profiles[l][1][i]*100:10.1f}" if i < len(profiles[l][1]) and not np.isnan(profiles[l][1][i]) else f"{'--':>10}" for l, _ in FILES)
    print(f"{centers[i]:10.1f} {row}")
