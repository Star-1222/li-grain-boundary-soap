#!/usr/bin/env python3
"""Fixed-threshold GB width: same rule for every structure.

The adaptive (plateau-vs-minimum) threshold used in analyze_structure.py is
not comparable across structures because the minimum Q differs between them.
Here the threshold is pinned at Q < 0.5 for all three structures, which is the
convention proposed in the bi-weekly report.
"""
import json
from pathlib import Path

import numpy as np

ANA = Path("exp_barostat/analysis")
NAMES = ["start_quenched", "final_xonly_fast", "final_aniso_fast"]

out = {}
for n in NAMES:
    z = np.load(ANA / f"{n}_profiles.npz")
    x, qC, qB = z["x"], z["qC"], z["qB"]
    bw = x[1] - x[0]
    finite = np.isfinite(qC)
    qf = np.interp(x, x[finite], qC[finite])
    w = float((qf < 0.5).sum() * bw)
    out[n] = {
        "w_GB_fixed50_A_coordination": w,
        "min_Q_coordination": float(qf.min()),
        "bulk_Q_coordination": float(0.5 * (qf[int(0.25 * len(qf))] +
                                            qf[int(0.75 * len(qf))])),
        "w_GB_fixed50_A_steinhardt": float(
            (np.interp(x, x[np.isfinite(qB)], qB[np.isfinite(qB)]) < 0.5).sum() * bw),
    }
    print(f"{n:20s} w_GB(Q<0.5, coordination) = {w:6.2f} A   "
          f"min Q = {out[n]['min_Q_coordination']:.3f}   "
          f"w_GB(Q<0.5, steinhardt) = {out[n]['w_GB_fixed50_A_steinhardt']:6.2f} A")

(ANA / "gb_width_fixed.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
print("\n[out]", ANA / "gb_width_fixed.json")
