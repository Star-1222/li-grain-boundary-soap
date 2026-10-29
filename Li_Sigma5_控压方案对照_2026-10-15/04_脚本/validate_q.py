#!/usr/bin/env python3
"""Validate the q4/q6 implementation against perfect BCC / FCC / HCP lattices
and locate the first minimum of g(r) in the quenched Li bicrystal.

Reported invariants for the *first shell only* depend on how many neighbours
are included, so we print k = 8, 12, 14 for BCC and k = 12 for FCC/HCP.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import KDTree
from scipy.special import sph_harm_y

NN = np.sqrt(3.0) / 2.0


def make_bcc(a=3.506, ncell=6):
    base = np.array([[0, 0, 0], [0.5, 0.5, 0.5]]) * a
    pos = []
    for i in range(ncell):
        for j in range(ncell):
            for k in range(ncell):
                pos.append(base + np.array([i, j, k]) * a)
    return np.array(pos).reshape(-1, 3), a * ncell


def make_fcc(a=4.0, ncell=4):
    base = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]) * a
    pos = []
    for i in range(ncell):
        for j in range(ncell):
            for k in range(ncell):
                pos.append(base + np.array([i, j, k]) * a)
    return np.array(pos).reshape(-1, 3), a * ncell


def make_hcp(a=3.2, ncell=4):
    # AA stacking repeated -> ideal hcp with c/a = sqrt(8/3)
    c = a * np.sqrt(8.0 / 3.0)
    base = np.array([[0, 0, 0],
                     [1.0 / 3.0, 2.0 / 3.0, 0.5],
                     [2.0 / 3.0, 1.0 / 3.0, 0.5]]) * np.array([a, a, c])
    # orthogonal cell: a x a*sqrt(3) x c  (hexagonal primitive -> orthorhombic 1x1)
    pos = []
    for i in range(ncell):
        for j in range(ncell):
            for k in range(ncell):
                pos.append(base + np.array([i * a, j * a * np.sqrt(3.0), k * c]))
    return np.array(pos).reshape(-1, 3), None


def steinhardt(pos, L, k, lvals=(4, 6)):
    tree = KDTree(pos, boxsize=L)
    _, idx = tree.query(pos, k=k + 1)
    nb = idx[:, 1:]
    disp = pos[nb] - pos[:, None, :]
    disp -= L * np.round(disp / L)
    dist = np.linalg.norm(disp, axis=2)
    unit = disp / dist[:, :, None]
    theta = np.arctan2(unit[:, :, 1], unit[:, :, 0])      # azimuthal
    phi_colat = np.arccos(np.clip(unit[:, :, 2], -1, 1))  # polar
    n = len(pos)
    out = {}
    for l in lvals:
        Y = sph_harm_y(l, np.arange(-l, l + 1)[None, :],
                       phi_colat.reshape(-1)[:, None], theta.reshape(-1)[:, None])
        Y = Y.reshape(n, k, 2 * l + 1)
        qbar = Y.mean(axis=1)
        out[l] = np.sqrt(4.0 * np.pi / (2 * l + 1) * np.sum(np.abs(qbar) ** 2, axis=1))
    return out


def report(tag, pos, Lcell, ks):
    for k in ks:
        if k >= len(pos):
            continue
        q = steinhardt(pos, Lcell if np.ndim(Lcell) else np.array([Lcell] * 3), k)
        print(f"  {tag:12s} k={k:2d}   q4 = {q[4].mean():.4f}   q6 = {q[6].mean():.4f}")


if __name__ == "__main__":
    pos, Lcell = make_bcc()
    print("perfect BCC (a=3.506):")
    report("BCC", pos, np.array([Lcell] * 3), [8, 12, 14])

    pos, Lcell = make_fcc()
    print("perfect FCC (a=4.0):")
    report("FCC", pos, np.array([Lcell] * 3), [12])

    try:
        pos, _ = make_hcp()
        Lcell = np.array([3.2 * 4, 3.2 * np.sqrt(3.0) * 4, 3.2 * np.sqrt(8 / 3) * 4])
        print("perfect HCP (a=3.2):")
        report("HCP", pos, Lcell, [12])
    except Exception as exc:
        print("HCP skipped:", exc)
