"""Duplicate-structure detection for repeated quenches (basin sampling).

[NEW 2026-08-21]

The metastable-GB workflow (v2 docs, "Scientific scope and limitations")
requires several independent anneal+quench replicas per GB to sample the
metastable-basin landscape.  Two relaxed replicas belong to the same
basin when their matched-atom SOAP similarity is essentially 1.  This
module performs greedy clustering of descriptor sets by that criterion.
"""

from __future__ import annotations

from typing import Dict, List, Sequence

import numpy as np

from .similarity import structure_similarity


def dedup_structures(descriptor_sets: Sequence[np.ndarray],
                     threshold: float = 0.999) -> Dict:
    """Greedy clustering of structures by matched-atom SOAP similarity.

    Parameters
    ----------
    descriptor_sets : list of (N, D) matrices, one per structure.  All
        structures must share N and atom ordering (same initial model).
    threshold       : structures with similarity >= threshold to an
        already-kept representative are considered duplicates.

    Returns
    -------
    {"representatives": [indices], "clusters": {rep_index: [member indices]},
     "similarities": {(i, j): sim for kept pairs}, "threshold": float}
    """
    if not descriptor_sets:
        return {"representatives": [], "clusters": {},
                "similarities": {}, "threshold": threshold}

    reps: List[int] = []
    clusters: Dict[int, List[int]] = {}
    sims: Dict[str, float] = {}
    for i, X in enumerate(descriptor_sets):
        assigned = None
        for r in reps:
            s = structure_similarity(descriptor_sets[r], X)
            sims[f"{r},{i}"] = s
            if s >= threshold:
                assigned = r
                break
        if assigned is None:
            reps.append(i)
            clusters[i] = [i]
        else:
            clusters[assigned].append(i)

    return {"representatives": reps, "clusters": clusters,
            "similarities": sims, "threshold": threshold}


def basin_report(descriptor_sets: Sequence[np.ndarray],
                 labels: Sequence[str],
                 threshold: float = 0.999) -> List[str]:
    """Human-readable dedup report (one line per structure)."""
    res = dedup_structures(descriptor_sets, threshold=threshold)
    lines = [f"dedup threshold (matched-atom SOAP similarity): {threshold}"]
    rep_of = {}
    for rep, members in res["clusters"].items():
        for m in members:
            rep_of[m] = rep
    for i, lab in enumerate(labels):
        tag = "REPRESENTATIVE" if res["representatives"] and i in res["representatives"] \
            else f"duplicate of {labels[rep_of[i]]}"
        lines.append(f"  {lab}: {tag}")
    n_basins = len(res["representatives"])
    lines.append(f"  -> {n_basins} distinct basin(s) among {len(labels)} structure(s)")
    return lines
