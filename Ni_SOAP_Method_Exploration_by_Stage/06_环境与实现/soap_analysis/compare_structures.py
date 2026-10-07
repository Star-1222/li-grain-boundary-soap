"""Cross-structure comparison helpers (initial vs relaxed, GB vs GB).

[NEW 2026-08-21]
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np

from .similarity import structure_similarity
from .disorder_profile import (disorder_profile, gb_width_regions,
                               pearson_correlation, per_atom_disorder)


def compare_two_structures(Xa: np.ndarray, Xb: np.ndarray,
                           bulk_vector: np.ndarray,
                           x_coords: Optional[np.ndarray] = None,
                           bin_width: float = 1.0,
                           box_x: Optional[float] = None) -> Dict:
    """Compare two states of the SAME model (matched atoms).

    Typical use: initial (stage-1 minimized) vs equilibrated (stage-4)
    GB -- research question: did annealing+quenching relax the GB into a
    distinct metastable basin, and how did the disorder profile change?

    Returns a dict with matched-atom similarity, mean/max disorder of
    each state, GB width regions, and the change of the profile peak.
    """
    sim = structure_similarity(Xa, Xb)
    da = per_atom_disorder(Xa, bulk_vector)
    db = per_atom_disorder(Xb, bulk_vector)
    out: Dict = {
        "matched_atom_similarity": sim,
        "disorder_mean_a": float(da.mean()),
        "disorder_mean_b": float(db.mean()),
        "disorder_max_a": float(da.max()),
        "disorder_max_b": float(db.max()),
        "disorder_mean_change": float(db.mean() - da.mean()),
    }
    if x_coords is not None:
        pa = disorder_profile(x_coords, da, bin_width=bin_width, box_x=box_x)
        pb = disorder_profile(x_coords, db, bin_width=bin_width, box_x=box_x)
        out["profile_peak_a"] = pa["profile_max"]
        out["profile_peak_b"] = pb["profile_max"]
        out["gb_regions_a"] = gb_width_regions(pa)
        out["gb_regions_b"] = gb_width_regions(pb)
        out["profile_peak_change"] = pb["profile_max"] - pa["profile_max"]
    return out


def disorder_vs_theta(summaries: Sequence[Dict]) -> List[Dict]:
    """Build the disorder-vs-misorientation table from per-GB summaries.

    Each input summary (run_analysis output "summary.json") must contain
    "theta_deg" and disorder statistics.  Returns rows sorted by theta.
    """
    rows: List[Dict] = []
    for s in summaries:
        if not isinstance(s, dict):
            continue
        if "theta_deg" not in s or "disorder" not in s:
            continue
        d = s["disorder"]
        row = {
            "label": s.get("label", s.get("structure", "?")),
            "theta_deg": s.get("theta_deg"),
            "sigma": s.get("sigma"),
            "plane": s.get("plane"),
            "natoms": s.get("natoms"),
            "disorder_mean": d.get("mean"),
            "disorder_max": d.get("max"),
            "disorder_gb_peak": s.get("gb", {}).get("profile_peak"),
            "gb_width_A": (s.get("gb", {}).get("regions") or [{}])[0].get("width_A"),
            "energy_disorder_r": s.get("energy", {}).get("pearson_r"),
        }
        rows.append(row)
    rows.sort(key=lambda r: (r["theta_deg"] is None, r["theta_deg"]))
    return rows
