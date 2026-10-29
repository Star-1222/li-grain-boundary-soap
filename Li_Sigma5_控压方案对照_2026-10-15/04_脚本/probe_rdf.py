#!/usr/bin/env python3
"""Probe g(r) and coordination number vs cutoff on one snapshot."""
import sys
import numpy as np

from analyze_structure import read_lammps_data, rdf, coordination, criterion_A
from scipy.spatial import KDTree

path = sys.argv[1] if len(sys.argv) > 1 else "stage3b_quench/stage3b_quenched.data"
pos, L = read_lammps_data(path)
tree = KDTree(pos, boxsize=L)
r, g = rdf(pos, L, tree, rmax=6.0, nbins=600)

print("== g(r) from 2.4 to 5.2 A ==")
sel = (r >= 2.4) & (r <= 5.2)
for rr, gg in zip(r[sel][::4], g[sel][::4]):
    bar = "#" * int(gg * 40)
    print(f"{rr:6.3f}  {gg:6.3f}  {bar}")

print()
print("== coordination number vs cutoff ==")
print(" rcut   <CN>    frac(CN==8)   frac(CN>=8)")
for rc in np.arange(2.9, 4.4, 0.1):
    cn = coordination(pos, L, tree, float(rc))
    print(f" {rc:4.2f}  {cn.mean():6.3f}   {np.mean(cn==8):8.4f}    {np.mean(cn>=8):8.4f}")

print()
print("== criterion A vs tolerance (a = box-implied) ==")
a_star = (2.0 * np.prod(L) / len(pos)) ** (1.0 / 3.0)
for tol in (0.02, 0.05, 0.08, 0.10, 0.15):
    isA, n1, n2 = criterion_A(pos, L, tree, a_star, tol)
    print(f" tol={tol:.2f}  frac={isA.mean():.4f}  <n1NN>={n1.mean():.2f}  <n2NN>={n2.mean():.2f}")
