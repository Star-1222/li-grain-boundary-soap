#!/usr/bin/env python3
"""Extract box dimensions vs step from a LAMMPS dump file."""
import sys
import numpy as np

path = sys.argv[1]
NAT = 14400
mode = None
steps, boxes = [], []
cur = []

with open(path) as fh:
    for line in fh:
        if line.startswith("ITEM: TIMESTEP"):
            mode = "step"
            continue
        if line.startswith("ITEM: BOX BOUNDS"):
            mode = "box"
            cur = []
            continue
        if line.startswith("ITEM: ATOMS") or line.startswith("ITEM: NUMBER"):
            mode = None
            continue
        if mode == "step":
            steps.append(int(line.split()[0]))
            mode = None
        elif mode == "box":
            p = line.split()
            if len(p) >= 2:
                cur.append(float(p[1]) - float(p[0]))
            if len(cur) == 3:
                boxes.append(cur)
                mode = None

steps = np.array(steps, dtype=int)
boxes = np.array(boxes, dtype=float)
n = min(len(steps), len(boxes))
steps, boxes = steps[:n], boxes[:n]

if len(boxes):
    L0 = boxes[0]
    print(f"{'step':>8} {'Lx':>10} {'Ly':>10} {'Lz':>10} "
          f"{'exx%':>9} {'eyy%':>9} {'ezz%':>9}")
    for s, b in zip(steps, boxes):
        e = (b - L0) / L0 * 100.0
        print(f"{s:8d} {b[0]:10.5f} {b[1]:10.5f} {b[2]:10.5f} "
              f"{e[0]:9.4f} {e[1]:9.4f} {e[2]:9.4f}")
    if len(boxes) >= 2:
        d = (boxes[-1] - boxes[-2]) / boxes[-2] * 100.0
        print(f"\nlast increment (frame-to-frame): "
              f"dLx={d[0]:+.4f}%  dLy={d[1]:+.4f}%  dLz={d[2]:+.4f}%")
