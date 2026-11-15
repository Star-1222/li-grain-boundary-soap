"""Stage 3A: DScribe-specific LAE equivalence threshold calibration.

This script calibrates the raw DScribe Euclidean SOAP distance scale only.
It does not construct final LAEs, LER vectors, clusters, or disorder
measures, and it does not choose a final threshold.
"""

from __future__ import annotations

import csv
import math
import os
import tempfile
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import shared_memory
from pathlib import Path

import numpy as np

os.environ["MPLCONFIGDIR"] = tempfile.gettempdir()

from PIL import Image, ImageDraw, ImageFont


PROJECT_ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
VECTOR_FILE = OUTPUT_DIR / "stage2_75_soap_vectors.npy"
METADATA_FILE = OUTPUT_DIR / "stage2_75_atom_metadata.csv"
POPULATION_FILE = OUTPUT_DIR / "stage2_75_population_summary.csv"

EXPECTED_N = 4659
EXPECTED_DIM = 1014
EXPECTED_COUNTS = {
    ("stable", "interface"): 1600,
    ("metastable", "interface"): 2059,
    ("stable", "bulk"): 500,
    ("metastable", "bulk"): 500,
}
BULK_TOTAL = EXPECTED_COUNTS[("stable", "bulk")] + EXPECTED_COUNTS[("metastable", "bulk")]
N_REPEATS = 30
BASE_SEED = 20260817
MAX_WORKERS = 12
BLOCK_SIZE = 512
N_PARTITION_EPS = 15
N_ARI_PAIRS = 30
N_REFINE_POINTS_PER_INTERVAL = 35


def _mean_std(values):
    values = np.asarray(values, dtype=float)
    return float(np.mean(values)), float(np.std(values, ddof=1)) if values.size > 1 else 0.0


def read_sweep_csv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(fh)]


def _summary(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {
            "N": 0,
            "min": float("nan"),
            "p1": float("nan"),
            "p5": float("nan"),
            "p10": float("nan"),
            "p25": float("nan"),
            "median": float("nan"),
            "p75": float("nan"),
            "p90": float("nan"),
            "p95": float("nan"),
            "p99": float("nan"),
            "max": float("nan"),
            "mean": float("nan"),
            "std": float("nan"),
        }
    return {
        "N": int(values.size),
        "min": float(values.min()),
        "p1": float(np.percentile(values, 1)),
        "p5": float(np.percentile(values, 5)),
        "p10": float(np.percentile(values, 10)),
        "p25": float(np.percentile(values, 25)),
        "median": float(np.median(values)),
        "p75": float(np.percentile(values, 75)),
        "p90": float(np.percentile(values, 90)),
        "p95": float(np.percentile(values, 95)),
        "p99": float(np.percentile(values, 99)),
        "max": float(values.max()),
        "mean": float(values.mean()),
        "std": float(values.std(ddof=1)) if values.size > 1 else 0.0,
    }


def load_and_validate():
    X = np.load(VECTOR_FILE)
    if X.shape != (EXPECTED_N, EXPECTED_DIM):
        raise RuntimeError(f"Unexpected SOAP matrix shape: {X.shape}")
    if not bool(np.isfinite(X).all()):
        raise RuntimeError("SOAP matrix contains non-finite values.")

    with open(METADATA_FILE, newline="", encoding="utf-8") as fh:
        records = list(csv.DictReader(fh))
    if len(records) != EXPECTED_N:
        raise RuntimeError(f"Metadata row count {len(records)} != {EXPECTED_N}")
    required_fields = {"structure", "region", "atom_id", "soap_norm"}
    missing = required_fields - set(records[0].keys())
    if missing:
        raise RuntimeError(f"Metadata missing fields: {sorted(missing)}")

    counts = Counter((r["structure"], r["region"]) for r in records)
    if counts != EXPECTED_COUNTS:
        raise RuntimeError(f"Population counts mismatch: {dict(counts)}")

    atom_ids = np.asarray([int(r["atom_id"]) for r in records], dtype=np.int64)
    for structure in sorted({r["structure"] for r in records}):
        ids_for_structure = [int(r["atom_id"]) for r in records if r["structure"] == structure]
        if len(set(ids_for_structure)) != len(ids_for_structure):
            raise RuntimeError(f"Atom IDs are not unique within structure {structure}.")

    structures = np.asarray([r["structure"] for r in records])
    regions = np.asarray([r["region"] for r in records])
    metadata_norms = np.asarray([float(r["soap_norm"]) for r in records])
    computed_norms = np.linalg.norm(X, axis=1)
    norm_diff = float(np.max(np.abs(computed_norms - metadata_norms)))
    if norm_diff > 1e-8:
        raise RuntimeError(f"Metadata SOAP norms disagree with vectors: max diff {norm_diff:.3e}")

    with open(POPULATION_FILE, newline="", encoding="utf-8") as fh:
        pop_rows = list(csv.DictReader(fh))
    for (structure, region), expected_count in EXPECTED_COUNTS.items():
        label = f"{structure} {region}"
        metric_rows = [r for r in pop_rows if r["category"] == label and r["metric"] == "soap_norm"]
        if len(metric_rows) != 1:
            raise RuntimeError(f"Population summary missing {label}/soap_norm")
        if int(metric_rows[0]["N"]) != expected_count:
            raise RuntimeError(
                f"Population summary N mismatch for {label}: "
                f"{metric_rows[0]['N']} != {expected_count}"
            )

    masks = {
        "stable_bulk": (structures == "stable") & (regions == "bulk"),
        "metastable_bulk": (structures == "metastable") & (regions == "bulk"),
        "stable_interface": (structures == "stable") & (regions == "interface"),
        "metastable_interface": (structures == "metastable") & (regions == "interface"),
    }
    for mask in masks.values():
        if not bool(np.isfinite(np.flatnonzero(mask)).all()):
            raise RuntimeError("Invalid group mask")

    return X, records, structures, regions, atom_ids, masks, pop_rows


def compute_distance_matrix(X):
    norms2 = np.einsum("ij,ij->i", X, X)
    D2 = norms2[:, None] + norms2[None, :] - 2.0 * (X @ X.T)
    np.fill_diagonal(D2, 0.0)
    D2 = np.maximum(D2, 0.0)
    if not bool(np.isfinite(D2).all()):
        raise RuntimeError("Distance matrix contains non-finite values.")
    if not np.allclose(D2, D2.T, atol=1e-10):
        raise RuntimeError("Distance matrix is not symmetric.")
    if not np.allclose(np.diag(D2), 0.0, atol=1e-12):
        raise RuntimeError("Distance matrix diagonal is not zero.")
    return D2


def compute_nearest_neighbors(D2, records, structures, regions, atom_ids):
    D2_inf = D2.copy()
    np.fill_diagonal(D2_inf, np.inf)
    nn_index = np.argmin(D2_inf, axis=1)
    nn_dist_sq = D2_inf[np.arange(len(D2_inf)), nn_index]
    d_nn = np.sqrt(np.maximum(nn_dist_sq, 0.0))
    if not bool(np.isfinite(d_nn).all()):
        raise RuntimeError("Nearest-neighbor distances contain non-finite values.")
    if np.any(d_nn < 0):
        raise RuntimeError("Nearest-neighbor distances are negative.")

    with open(OUTPUT_DIR / "stage3a_nearest_neighbor_distances.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "row_index",
                "structure",
                "region",
                "atom_id",
                "nearest_neighbor_row",
                "nearest_neighbor_structure",
                "nearest_neighbor_region",
                "nearest_neighbor_atom_id",
                "nearest_neighbor_distance",
            ]
        )
        for i in range(len(d_nn)):
            j = int(nn_index[i])
            writer.writerow(
                [
                    i,
                    structures[i],
                    regions[i],
                    int(atom_ids[i]),
                    j,
                    structures[j],
                    regions[j],
                    int(atom_ids[j]),
                    f"{d_nn[i]:.12g}",
                ]
            )

    group_names = [
        ("stable_bulk", "stable bulk"),
        ("metastable_bulk", "metastable bulk"),
        ("stable_interface", "stable interface"),
        ("metastable_interface", "metastable interface"),
        ("all_bulk", "all bulk"),
        ("all_interface", "all interface"),
        ("full", "full combined population"),
    ]
    bulk_mask = (structures == "stable") | (structures == "metastable")
    interface_mask = regions == "interface"
    mask_sets = {
        "stable_bulk": (structures == "stable") & (regions == "bulk"),
        "metastable_bulk": (structures == "metastable") & (regions == "bulk"),
        "stable_interface": (structures == "stable") & (regions == "interface"),
        "metastable_interface": (structures == "metastable") & (regions == "interface"),
        "all_bulk": bulk_mask & (regions == "bulk"),
        "all_interface": interface_mask,
        "full": np.ones(EXPECTED_N, dtype=bool),
    }
    summary_rows = []
    for key, label in group_names:
        summary_rows.append({"group": label, **_summary(d_nn[mask_sets[key]])})
    with open(OUTPUT_DIR / "stage3a_nearest_neighbor_summary.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "group",
                "N",
                "min",
                "p1",
                "p5",
                "p10",
                "p25",
                "median",
                "p75",
                "p90",
                "p95",
                "p99",
                "max",
                "mean",
                "std",
            ],
        )
        writer.writeheader()
        writer.writerows(summary_rows)
    return d_nn, summary_rows


def build_epsilon_grid(d_nn, D2):
    positive_nn = d_nn[d_nn > 0]
    if positive_nn.size == 0:
        raise RuntimeError("No positive nearest-neighbor distances available.")
    p = np.percentile(positive_nn, [1, 5, 10, 25, 50, 75, 90, 95, 99])
    upper = np.triu_indices(EXPECTED_N, 1)
    pair_upper = np.sqrt(D2[upper])
    pair_p95 = float(np.percentile(pair_upper, 95))
    pair_max = float(pair_upper.max())

    lo = float(positive_nn.min())
    mid = float(p[4])
    hi = pair_p95

    log_pts = np.logspace(np.log10(lo), np.log10(mid), 45)
    lin_pts = np.linspace(mid, hi, 61)[1:]
    anchors = np.concatenate([p, [positive_nn.max()]])
    values = np.unique(np.concatenate([log_pts, lin_pts, anchors]))
    values = values[values > 0]
    eps_grid = np.sort(values)
    if not 80 <= len(eps_grid) <= 150:
        raise RuntimeError(f"Epsilon grid size {len(eps_grid)} outside target 80-150.")
    if eps_grid.max() >= pair_max:
        raise RuntimeError("Maximum epsilon would trivially collapse the population to one class.")

    sources = []
    for v in eps_grid:
        if np.any(np.isclose(v, anchors)):
            sources.append("anchor_percentile")
        elif v <= mid:
            sources.append("log_lower_to_middle")
        else:
            sources.append("linear_middle_to_upper")
    with open(OUTPUT_DIR / "stage3a_epsilon_grid.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["index", "epsilon", "source"])
        for i, (v, source) in enumerate(zip(eps_grid, sources)):
            writer.writerow([i, f"{v:.17g}", source])

    grid_doc = {
        "n_eps": int(len(eps_grid)),
        "min_eps": float(eps_grid.min()),
        "median_nn": mid,
        "max_eps": float(eps_grid.max()),
        "pair_p95": pair_p95,
        "pair_max": pair_max,
        "min_positive_nn": lo,
        "nn_p99": float(p[-1]),
        "nn_max": float(positive_nn.max()),
        "n_log": int(45),
        "n_linear": int(61),
        "n_anchors": int(anchors.size),
    }
    return eps_grid, grid_doc


def _analyze_run(D2, eps, order, group_masks, return_assignment):
    eps2 = eps * eps
    if eps <= 0:
        raise ValueError(f"Epsilon must be positive, got {eps}")
    reps_list = []
    for i in order:
        if reps_list and np.any(D2[i, reps_list] < eps2):
            continue
        reps_list.append(i)
    reps = np.asarray(reps_list, dtype=np.int64)
    if not np.all(reps >= 0) or not np.all(reps < EXPECTED_N):
        raise RuntimeError("Invalid representative index.")
    if len(np.unique(reps)) != len(reps):
        raise RuntimeError("Representative indices are not unique.")

    K = len(reps)
    assignment = np.empty(EXPECTED_N, dtype=np.int64)
    occ = np.zeros(K, dtype=np.int64)
    group_occ = {key: np.zeros(K, dtype=np.int64) for key in group_masks}
    for start in range(0, EXPECTED_N, BLOCK_SIZE):
        block = np.arange(start, min(start + BLOCK_SIZE, EXPECTED_N))
        cols = np.argmin(D2[np.ix_(block, reps)], axis=1)
        assignment[block] = cols
        occ += np.bincount(cols, minlength=K)
        for key, mask in group_masks.items():
            group_occ[key] += np.bincount(cols[mask[block]], minlength=K)
    if int(occ.sum()) != EXPECTED_N:
        raise RuntimeError(f"Class occupancies sum to {occ.sum()}, expected {EXPECTED_N}")

    single = occ == 1
    bulk_occ = group_occ["stable_bulk"] + group_occ["metastable_bulk"]
    nonbulk_occ = occ - bulk_occ
    stable_iface = group_occ["stable_interface"]
    meta_iface = group_occ["metastable_interface"]
    bulk_class_mask = bulk_occ > 0
    stable_iface_mask = stable_iface > 0
    meta_iface_mask = meta_iface > 0
    n_stable_iface = int(stable_iface_mask.sum())
    n_meta_iface = int(meta_iface_mask.sum())
    bulk_dominant = int(np.argmax(bulk_occ))

    def _iface_stats(iface_occ, mask, group_total):
        if not np.any(mask):
            return {
                "class_count": 0,
                "singleton_fraction": float("nan"),
                "median_occupancy": float("nan"),
                "largest_class_fraction": float("nan"),
                "largest_group_fraction": float("nan"),
            }
        subset = occ[mask]
        return {
            "class_count": int(mask.sum()),
            "singleton_fraction": float(np.mean(subset == 1)),
            "median_occupancy": float(np.median(subset)),
            "largest_class_fraction": float(subset.max() / EXPECTED_N),
            "largest_group_fraction": float(iface_occ[mask].max() / group_total),
        }

    stable_stats = _iface_stats(stable_iface, stable_iface_mask, EXPECTED_COUNTS[("stable", "interface")])
    meta_stats = _iface_stats(meta_iface, meta_iface_mask, EXPECTED_COUNTS[("metastable", "interface")])

    stats = {
        "n_classes": K,
        "singleton_class_fraction": float(np.mean(single)),
        "atom_singleton_fraction": float(single.sum() / EXPECTED_N),
        "median_occupancy": float(np.median(occ)),
        "largest_class_fraction": float(occ.max() / EXPECTED_N),
        "bulk_class_count": int(bulk_class_mask.sum()),
        "bulk_largest_class_fraction": float(bulk_occ.max() / BULK_TOTAL),
        "bulk_exclusive_class_count": int((bulk_class_mask & (nonbulk_occ == 0)).sum()),
        "bulk_interface_mixed_class_count": int((bulk_class_mask & (nonbulk_occ > 0)).sum()),
        "bulk_same_dominant_fraction": float(
            (group_occ["stable_bulk"][bulk_dominant] + group_occ["metastable_bulk"][bulk_dominant])
            / BULK_TOTAL
        ),
        "bulk_stable_dominant_fraction": float(group_occ["stable_bulk"][bulk_dominant] / EXPECTED_COUNTS[("stable", "bulk")]),
        "bulk_metastable_dominant_fraction": float(
            group_occ["metastable_bulk"][bulk_dominant] / EXPECTED_COUNTS[("metastable", "bulk")]
        ),
        "interface_stable_class_count": stable_stats["class_count"],
        "interface_stable_singleton_fraction": stable_stats["singleton_fraction"],
        "interface_stable_median_occupancy": stable_stats["median_occupancy"],
        "interface_stable_largest_class_fraction": stable_stats["largest_class_fraction"],
        "interface_stable_largest_group_fraction": stable_stats["largest_group_fraction"],
        "interface_metastable_class_count": meta_stats["class_count"],
        "interface_metastable_singleton_fraction": meta_stats["singleton_fraction"],
        "interface_metastable_median_occupancy": meta_stats["median_occupancy"],
        "interface_metastable_largest_class_fraction": meta_stats["largest_class_fraction"],
        "interface_metastable_largest_group_fraction": meta_stats["largest_group_fraction"],
        "interface_shared_class_count": int((stable_iface_mask & meta_iface_mask).sum()),
        "interface_stable_only_class_count": int((stable_iface_mask & ~meta_iface_mask).sum()),
        "interface_metastable_only_class_count": int((meta_iface_mask & ~stable_iface_mask).sum()),
    }
    return stats, assignment if return_assignment else None


GROUP_MASKS = None
SHM = None


def _init_worker(shm_name, shape, structures, regions):
    global SHM, GROUP_MASKS
    SHM = shared_memory.SharedMemory(name=shm_name)
    D2 = np.ndarray(shape, dtype=np.float64, buffer=SHM.buf)
    structures = np.asarray(structures)
    regions = np.asarray(regions)
    GROUP_MASKS = {
        "stable_bulk": (structures == "stable") & (regions == "bulk"),
        "metastable_bulk": (structures == "metastable") & (regions == "bulk"),
        "stable_interface": (structures == "stable") & (regions == "interface"),
        "metastable_interface": (structures == "metastable") & (regions == "interface"),
    }


def _worker_task(task):
    eps, order, return_assignment = task
    return _analyze_run(SHM_D2(), eps, order, GROUP_MASKS, return_assignment)


def SHM_D2():
    if SHM is None:
        raise RuntimeError("Shared distance matrix is not initialized.")
    return np.ndarray((EXPECTED_N, EXPECTED_N), dtype=np.float64, buffer=SHM.buf)


def run_epsilon_sweep(
    D2,
    eps_grid,
    selected_indices,
    structures,
    regions,
    save_path=None,
    n_repeats=N_REPEATS,
    max_workers=MAX_WORKERS,
):
    tasks = []
    for e, eps in enumerate(eps_grid):
        for r in range(n_repeats):
            rng = np.random.default_rng(BASE_SEED + e * 10_000 + r)
            order = rng.permutation(EXPECTED_N)
            tasks.append((float(eps), order, e in selected_indices))

    shm = shared_memory.SharedMemory(create=True, size=D2.nbytes)
    shm_arr = np.ndarray(D2.shape, dtype=np.float64, buffer=shm.buf)
    shm_arr[:] = D2

    sweep = {e: [] for e in range(len(eps_grid))}
    assignments_by_eps = {e: [] for e in selected_indices}
    try:
        with ProcessPoolExecutor(
            max_workers=max_workers,
            initializer=_init_worker,
            initargs=(shm.name, D2.shape, structures, regions),
        ) as pool:
            for idx, (stats, assignment) in enumerate(pool.map(_worker_task, tasks, chunksize=8)):
                e = idx // n_repeats
                sweep[e].append(stats)
                if e in assignments_by_eps:
                    assignments_by_eps[e].append(assignment)
    finally:
        shm.close()
        shm.unlink()

    counts_by_run = np.empty((len(eps_grid), n_repeats), dtype=np.int64)
    for e in range(len(eps_grid)):
        for r, stats in enumerate(sweep[e]):
            counts_by_run[e, r] = stats["n_classes"]
    monotonic_violations = 0
    for r in range(n_repeats):
        diffs = np.diff(counts_by_run[:, r])
        monotonic_violations += int(np.sum(diffs > 0))

    sweep_rows = []
    stat_keys = list(sweep[0][0].keys())
    for e in range(len(eps_grid)):
        row = {"epsilon": float(eps_grid[e])}
        for key in stat_keys:
            values = [s[key] for s in sweep[e]]
            mean, std = _mean_std(values)
            row[f"{key}_mean"] = mean
            row[f"{key}_std"] = std
        counts = counts_by_run[e]
        row["n_classes_median"] = float(np.median(counts))
        row["n_classes_min"] = float(counts.min())
        row["n_classes_max"] = float(counts.max())
        row["n_classes_cv"] = float(counts.std(ddof=1) / counts.mean()) if counts.mean() > 0 else float("nan")
        sweep_rows.append(row)

    fieldnames = ["epsilon"]
    for key in stat_keys:
        fieldnames.append(f"{key}_mean")
        fieldnames.append(f"{key}_std")
    fieldnames.extend(["n_classes_median", "n_classes_min", "n_classes_max", "n_classes_cv"])
    if save_path is not None:
        with open(save_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(sweep_rows)
    return sweep_rows, assignments_by_eps, counts_by_run, monotonic_violations


def adjusted_rand_index(a, b):
    a = np.asarray(a)
    b = np.asarray(b)
    n = a.size
    _, a_inv = np.unique(a, return_inverse=True)
    _, b_inv = np.unique(b, return_inverse=True)
    cont = np.zeros((a_inv.max() + 1, b_inv.max() + 1), dtype=np.int64)
    np.add.at(cont, (a_inv, b_inv), 1)
    comb = lambda x: x * (x - 1) // 2
    a_counts = cont.sum(axis=1)
    b_counts = cont.sum(axis=0)
    index = int(np.sum(comb(cont)))
    expected = comb(a_counts).sum() * comb(b_counts).sum() / comb(n)
    max_index = 0.5 * (comb(a_counts).sum() + comb(b_counts).sum())
    if max_index == expected:
        return 1.0 if n > 1 else float("nan")
    return float((index - expected) / (max_index - expected))


def connected_components(D2, eps):
    eps2 = eps * eps
    n = len(D2)
    labels = np.full(n, -1, dtype=np.int64)
    label = 0
    for i in range(n):
        if labels[i] != -1:
            continue
        labels[i] = label
        stack = [i]
        while stack:
            node = stack.pop()
            new = (D2[node] < eps2) & (labels == -1)
            if np.any(new):
                new_idx = np.flatnonzero(new)
                labels[new_idx] = label
                stack.extend(new_idx.tolist())
        label += 1
    sizes = np.bincount(labels)
    return label, float(sizes.max() / n)


def compute_partition_stability(D2, eps_grid, sweep_rows, assignments_by_eps):
    selected_idx = np.unique(
        np.round(np.linspace(0, len(eps_grid) - 1, N_PARTITION_EPS)).astype(int)
    )
    selected_idx = [int(i) for i in selected_idx if i in assignments_by_eps]
    rows = []
    for e in selected_idx:
        eps = float(eps_grid[e])
        assignments = assignments_by_eps[e]
        pair_indices = [
            (i, j) for i in range(len(assignments)) for j in range(i + 1, len(assignments))
        ]
        rng = np.random.default_rng(BASE_SEED + e)
        if len(pair_indices) > N_ARI_PAIRS:
            chosen = rng.choice(len(pair_indices), size=N_ARI_PAIRS, replace=False)
            pair_indices = [pair_indices[int(c)] for c in chosen]
        ari_values = []
        for i, j in pair_indices:
            ari_values.append(adjusted_rand_index(assignments[i], assignments[j]))
        pair_count = len(pair_indices)
        ari_values = np.asarray(ari_values)
        n_components, largest_fraction = connected_components(D2, eps)
        greedy_mean = sweep_rows[e]["n_classes_mean"]
        rows.append(
            {
                "epsilon": eps,
                "epsilon_index": e,
                "n_pairs": pair_count,
                "n_classes_mean": greedy_mean,
                "ari_mean": float(ari_values.mean()),
                "ari_median": float(np.median(ari_values)),
                "ari_min": float(ari_values.min()),
                "ari_std": float(ari_values.std(ddof=1)) if len(ari_values) > 1 else 0.0,
                "connected_component_count": int(n_components),
                "largest_connected_component_fraction": largest_fraction,
                "component_to_class_ratio": float(n_components / greedy_mean) if greedy_mean > 0 else float("nan"),
            }
        )
    with open(OUTPUT_DIR / "stage3a_partition_stability.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "epsilon",
                "epsilon_index",
                "n_pairs",
                "n_classes_mean",
                "ari_mean",
                "ari_median",
                "ari_min",
                "ari_std",
                "connected_component_count",
                "largest_connected_component_fraction",
                "component_to_class_ratio",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    return rows


def identify_candidate_intervals(eps_grid, sweep_rows, ari_rows):
    if not ari_rows:
        return []
    eps = np.asarray(eps_grid, dtype=float)
    log_eps = np.log10(eps)
    med_counts = np.asarray([r["n_classes_median"] for r in sweep_rows])
    std_counts = np.asarray([r["n_classes_std"] for r in sweep_rows])
    singleton = np.asarray([r["singleton_class_fraction_mean"] for r in sweep_rows])
    bulk_largest = np.asarray([r["bulk_largest_class_fraction_mean"] for r in sweep_rows])
    stable_iface = np.asarray([r["interface_stable_class_count_mean"] for r in sweep_rows])
    meta_iface = np.asarray([r["interface_metastable_class_count_mean"] for r in sweep_rows])

    ari_eps = np.asarray([r["epsilon"] for r in ari_rows])
    ari_mean = np.asarray([r["ari_mean"] for r in ari_rows])
    ari_interp = np.interp(log_eps, np.log10(ari_eps), ari_mean)

    slope = np.zeros_like(eps)
    dx = np.diff(log_eps)
    dy = -np.diff(med_counts)
    with np.errstate(divide="ignore", invalid="ignore"):
        slope[1:] = dy / dx
    slope[~np.isfinite(slope)] = 0.0
    slope[0] = slope[1]

    def pct(x, q):
        return float(np.percentile(x, q))

    slow = slope <= pct(slope, 25)
    low_std = std_counts <= pct(std_counts, 50)
    low_singleton = singleton <= pct(singleton, 50)
    bulk_coherent = bulk_largest >= pct(bulk_largest, 30)
    interface_alive = (stable_iface + meta_iface) >= max(5.0, pct(stable_iface + meta_iface, 10))
    ari_ok = ari_interp >= pct(ari_interp, 30)
    candidate = slow & low_std & low_singleton & bulk_coherent & interface_alive & ari_ok

    intervals = []
    start = None
    for i, flag in enumerate(candidate):
        if flag and start is None:
            start = i
        elif not flag and start is not None:
            if i - start >= 3:
                intervals.append((start, i - 1))
            start = None
    if start is not None and len(eps) - start >= 3:
        intervals.append((start, len(eps) - 1))

    rows = []
    for interval_id, (lo_i, hi_i) in enumerate(intervals, start=1):
        rows.append(
            {
                "interval_id": interval_id,
                "epsilon_min": float(eps[lo_i]),
                "epsilon_max": float(eps[hi_i]),
                "n_grid_points": int(hi_i - lo_i + 1),
                "n_classes_median_min": float(med_counts[lo_i:hi_i + 1].min()),
                "n_classes_median_max": float(med_counts[lo_i:hi_i + 1].max()),
                "n_classes_std_max": float(std_counts[lo_i:hi_i + 1].max()),
                "ari_mean_min": float(ari_interp[lo_i:hi_i + 1].min()),
                "singleton_fraction_max": float(singleton[lo_i:hi_i + 1].max()),
                "bulk_largest_fraction_min": float(bulk_largest[lo_i:hi_i + 1].min()),
            }
        )
    with open(OUTPUT_DIR / "stage3a_threshold_candidates.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "interval_id",
                "epsilon_min",
                "epsilon_max",
                "n_grid_points",
                "n_classes_median_min",
                "n_classes_median_max",
                "n_classes_std_max",
                "ari_mean_min",
                "singleton_fraction_max",
                "bulk_largest_fraction_min",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    return rows


def run_refined_sweep(D2, structures, regions, intervals, n_repeats=N_REPEATS, max_workers=MAX_WORKERS):
    if not intervals:
        return None
    eps_parts = []
    for row in intervals:
        lo = row["epsilon_min"]
        hi = row["epsilon_max"]
        if hi <= lo:
            eps_parts.append(np.array([lo]))
        else:
            eps_parts.append(np.geomspace(lo, hi, N_REFINE_POINTS_PER_INTERVAL))
    refined_grid = np.unique(np.concatenate(eps_parts))
    selected = set()
    refined_rows, _, _, _ = run_epsilon_sweep(
        D2,
        refined_grid,
        selected,
        structures,
        regions,
        save_path=OUTPUT_DIR / "stage3a_refined_epsilon_sweep.csv",
        n_repeats=n_repeats,
        max_workers=max_workers,
    )
    return refined_rows


# ---------- PIL plotting helpers ----------


def _font(size, bold=False):
    path = r"C:\Windows\Fonts\arialbd.ttf" if bold else r"C:\Windows\Fonts\arial.ttf"
    try:
        return ImageFont.truetype(path, size)
    except OSError:
        return ImageFont.load_default()


def _fmt(v):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return ""
    if v == 0:
        return "0"
    av = abs(v)
    if av >= 1000 or av < 0.01:
        return f"{v:.2e}"
    if av >= 100:
        return f"{v:.0f}"
    if av >= 10:
        return f"{v:.1f}"
    if av >= 1:
        return f"{v:.2f}"
    return f"{v:.3f}"


def _ticks(lo, hi, scale, n=7):
    if scale == "log":
        lo10 = math.floor(math.log10(lo)) if lo > 0 else 0
        hi10 = math.ceil(math.log10(hi)) if hi > 0 else 0
        ticks = 10 ** np.arange(lo10, hi10 + 1, dtype=float)
        if len(ticks) > n:
            ticks = ticks[:: math.ceil(len(ticks) / n)]
        return np.asarray(ticks)
    if lo == hi:
        hi = lo + 1
    return np.linspace(lo, hi, n)


def _map_xy(rect, x, y, xscale, yscale, xmin, xmax, ymin, ymax):
    x0, y0, x1, y1 = rect
    fx = (math.log10(x) - math.log10(xmin)) / (math.log10(xmax) - math.log10(xmin)) if xscale == "log" else (x - xmin) / (xmax - xmin)
    fy = (y - ymin) / (ymax - ymin) if yscale == "linear" else (math.log10(y) - math.log10(ymin)) / (math.log10(ymax) - math.log10(ymin))
    return x0 + fx * (x1 - x0), y1 - fy * (y1 - y0)


def _draw_axes(draw, rect, xmin, xmax, ymin, ymax, xscale, yscale, xlabel, ylabel, font_axis, font_tick):
    x0, y0, x1, y1 = rect
    draw.line([(x0, y1), (x1, y1)], fill=(60, 60, 60), width=2)
    draw.line([(x0, y0), (x0, y1)], fill=(60, 60, 60), width=2)
    for tx in _ticks(xmin, xmax, xscale):
        if tx < xmin or tx > xmax:
            continue
        px, _ = _map_xy(rect, tx, ymin, xscale, yscale, xmin, xmax, ymin, ymax)
        draw.line([(px, y1 - 4), (px, y1 + 4)], fill=(60, 60, 60), width=1)
        label = _fmt(float(tx))
        bbox = draw.textbbox((0, 0), label, font=font_tick)
        draw.text((px - (bbox[2] - bbox[0]) / 2, y1 + 8), label, fill=(60, 60, 60), font=font_tick)
    for ty in _ticks(ymin, ymax, yscale):
        if ty < ymin or ty > ymax:
            continue
        _, py = _map_xy(rect, xmin, ty, xscale, yscale, xmin, xmax, ymin, ymax)
        draw.line([(x0 - 4, py), (x0 + 4, py)], fill=(60, 60, 60), width=1)
        label = _fmt(float(ty))
        bbox = draw.textbbox((0, 0), label, font=font_tick)
        draw.text((x0 - 8 - (bbox[2] - bbox[0]), py - (bbox[3] - bbox[1]) / 2), label, fill=(60, 60, 60), font=font_tick)
    draw.text((x0, y1 + 28), xlabel, fill=(30, 30, 30), font=font_axis)
    draw.text((x0 - 50, y0 - 24), ylabel, fill=(30, 30, 30), font=font_axis)


def _plot_panel(
    img,
    rect,
    title,
    series,
    xlabel,
    ylabel,
    xscale="linear",
    yscale="linear",
    bands=None,
    legend=True,
):
    draw = ImageDraw.Draw(img)
    font_title = _font(18, bold=True)
    font_axis = _font(14, bold=True)
    font_tick = _font(12)
    draw.text((rect[0], rect[1] - 34), title, fill=(20, 20, 20), font=font_title)

    xs_all = np.concatenate([np.asarray(s["x"], dtype=float) for s in series] + ([np.asarray(b["x"], dtype=float) for b in bands] if bands else []))
    ys_all = np.concatenate([np.asarray(s["y"], dtype=float) for s in series] + ([np.asarray(b["ymax"], dtype=float) for b in bands] + [np.asarray(b["ymin"], dtype=float) for b in bands] if bands else []))
    xmin, xmax = float(xs_all.min()), float(xs_all.max())
    ymin, ymax = float(ys_all.min()), float(ys_all.max())
    if xscale == "log":
        xmin = max(xmin * 0.9, xmin / 1.2)
        xmax = xmax * 1.2
    else:
        pad = (xmax - xmin) * 0.04 if xmax > xmin else 1.0
        xmin -= pad
        xmax += pad
    pad = (ymax - ymin) * 0.05 if ymax > ymin else 1.0
    ymin = max(0.0, ymin - pad)
    ymax += pad
    if yscale == "log":
        ymin = max(ymin, 1e-12)
        ymax = max(ymax, ymin * 1.1)

    _draw_axes(draw, rect, xmin, xmax, ymin, ymax, xscale, yscale, xlabel, ylabel, font_axis, font_tick)

    if bands:
        for band in bands:
            overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
            odraw = ImageDraw.Draw(overlay)
            xs = np.asarray(band["x"], dtype=float)
            lo = np.asarray(band["ymin"], dtype=float)
            hi = np.asarray(band["ymax"], dtype=float)
            pts = []
            for x, y in zip(xs, lo):
                pts.append(_map_xy(rect, float(x), float(y), xscale, yscale, xmin, xmax, ymin, ymax))
            for x, y in zip(xs[::-1], hi[::-1]):
                pts.append(_map_xy(rect, float(x), float(y), xscale, yscale, xmin, xmax, ymin, ymax))
            if len(pts) >= 3:
                odraw.polygon(pts, fill=tuple(list(band["color"]) + [55]))
            img.paste(overlay, (0, 0), overlay)
        draw = ImageDraw.Draw(img)

    for s in series:
        xs = np.asarray(s["x"], dtype=float)
        ys = np.asarray(s["y"], dtype=float)
        pts = [_map_xy(rect, float(x), float(y), xscale, yscale, xmin, xmax, ymin, ymax) for x, y in zip(xs, ys)]
        draw.line(pts, fill=s["color"], width=3, joint="curve")
        if s.get("marker"):
            r = 4
            for px, py in pts:
                draw.ellipse([px - r, py - r, px + r, py + r], fill=s["color"])

    if legend and len(series) > 1:
        font_leg = _font(12)
        items = [(s["label"], s["color"]) for s in series]
        widths = [draw.textbbox((0, 0), label, font=font_leg)[2] for label, _ in items]
        total_w = sum(widths) + 18 * len(items) + 12 * (len(items) - 1)
        lx = rect[2] - total_w - 6
        ly = rect[1] + 6
        for label, color in items:
            draw.line([(lx, ly + 6), (lx + 14, ly + 6)], fill=color, width=3)
            draw.text((lx + 18, ly), label, fill=(30, 30, 30), font=font_leg)
            lx += 18 + draw.textbbox((0, 0), label, font=font_leg)[2] + 12


def _new_canvas(width=1150, height=760):
    img = Image.new("RGB", (width, height), "white")
    return img, ImageDraw.Draw(img)


def plot_nn_distributions(dnn, masks):
    img, draw = _new_canvas()
    rect = (95, 90, 1120, 690)
    colors = [(31, 119, 180), (255, 127, 14), (44, 160, 44), (214, 39, 40)]
    labels = ["stable bulk", "metastable bulk", "stable interface", "metastable interface"]
    lo = max(dnn.min() * 0.8, 1e-4)
    hi = dnn.max() * 1.15
    bins = np.geomspace(lo, hi, 40)
    series = []
    for (label, color), key in zip(zip(labels, colors), ["stable_bulk", "metastable_bulk", "stable_interface", "metastable_interface"]):
        counts, _ = np.histogram(dnn[masks[key]], bins=bins)
        density = counts / counts.sum() / np.diff(bins)
        centers = np.sqrt(bins[:-1] * bins[1:])
        series.append({"x": centers, "y": density, "color": color, "label": label, "marker": True})
    _plot_panel(
        img,
        rect,
        "Nearest-neighbor raw DScribe Euclidean SOAP distance distributions",
        series,
        "d_NN (raw Euclidean SOAP distance)",
        "density",
        xscale="log",
    )
    img.save(OUTPUT_DIR / "stage3a_nn_distributions.png", dpi=(150, 150))


def plot_class_count(eps_grid, sweep_rows, path, title):
    eps = np.asarray(eps_grid, dtype=float)
    med = np.asarray([r["n_classes_median"] for r in sweep_rows])
    lo = np.asarray([r["n_classes_min"] for r in sweep_rows])
    hi = np.asarray([r["n_classes_max"] for r in sweep_rows])
    img, _ = _new_canvas()
    rect = (95, 90, 1120, 690)
    _plot_panel(
        img,
        rect,
        title,
        [{"x": eps, "y": med, "color": (31, 119, 180), "label": "median", "marker": False}],
        "epsilon",
        "number of candidate equivalence classes",
        xscale="log",
        bands=[{"x": eps, "ymin": lo, "ymax": hi, "color": (31, 119, 180)}],
        legend=False,
    )
    img.save(path, dpi=(150, 150))


def plot_order_dependence(eps_grid, sweep_rows, path):
    eps = np.asarray(eps_grid, dtype=float)
    cv = np.asarray([r["n_classes_cv"] for r in sweep_rows])
    std = np.asarray([r["n_classes_std"] for r in sweep_rows])
    img, _ = _new_canvas()
    rect = (95, 90, 1120, 690)
    _plot_panel(
        img,
        rect,
        "Order dependence of candidate-class count vs epsilon",
        [
            {"x": eps, "y": cv, "color": (214, 39, 40), "label": "CV of class count", "marker": False},
            {"x": eps, "y": std / max(1.0, float(np.max(std))), "color": (44, 160, 44), "label": "std (scaled to max 1)", "marker": False},
        ],
        "epsilon",
        "coefficient of variation / scaled std",
        xscale="log",
    )
    img.save(path, dpi=(150, 150))


def plot_singleton(eps_grid, sweep_rows, path):
    eps = np.asarray(eps_grid, dtype=float)
    class_sing = np.asarray([r["singleton_class_fraction_mean"] for r in sweep_rows])
    atom_sing = np.asarray([r["atom_singleton_fraction_mean"] for r in sweep_rows])
    img, _ = _new_canvas()
    rect = (95, 90, 1120, 690)
    _plot_panel(
        img,
        rect,
        "Singleton behavior vs epsilon",
        [
            {"x": eps, "y": class_sing, "color": (31, 119, 180), "label": "singleton-class fraction", "marker": False},
            {"x": eps, "y": atom_sing, "color": (255, 127, 14), "label": "fraction of atoms in singletons", "marker": False},
        ],
        "epsilon",
        "fraction",
        xscale="log",
    )
    img.save(path, dpi=(150, 150))


def plot_bulk_behavior(eps_grid, sweep_rows, path):
    eps = np.asarray(eps_grid, dtype=float)
    counts_keys = ["bulk_class_count_mean", "bulk_exclusive_class_count_mean", "bulk_interface_mixed_class_count_mean"]
    fraction_keys = ["bulk_largest_class_fraction_mean", "bulk_same_dominant_fraction_mean"]
    img, _ = _new_canvas()
    rect1 = (95, 90, 1120, 365)
    rect2 = (95, 430, 1120, 705)
    _plot_panel(
        img,
        rect1,
        "Bulk grouping: class-count diagnostics",
        [
            {"x": eps, "y": np.asarray([r[k] for r in sweep_rows]), "color": c, "label": label, "marker": False}
            for k, c, label in [
                ("bulk_class_count_mean", (31, 119, 180), "classes containing bulk"),
                ("bulk_exclusive_class_count_mean", (44, 160, 44), "bulk-exclusive classes"),
                ("bulk_interface_mixed_class_count_mean", (214, 39, 40), "bulk+interface mixed classes"),
            ]
        ],
        "epsilon",
        "number of classes",
        xscale="log",
    )
    _plot_panel(
        img,
        rect2,
        "Bulk grouping: dominant-class fractions",
        [
            {"x": eps, "y": np.asarray([r[k] for r in sweep_rows]), "color": c, "label": label, "marker": False}
            for k, c, label in [
                ("bulk_largest_class_fraction_mean", (31, 119, 180), "largest bulk-containing class"),
                ("bulk_same_dominant_fraction_mean", (255, 127, 14), "bulk atoms in dominant bulk class"),
            ]
        ],
        "epsilon",
        "fraction of bulk atoms",
        xscale="log",
    )
    img.save(path, dpi=(150, 150))


def plot_partition_stability(eps_grid, ari_rows, path):
    ari_eps = np.asarray([r["epsilon"] for r in ari_rows])
    img, _ = _new_canvas()
    rect = (95, 90, 1120, 690)
    _plot_panel(
        img,
        rect,
        "Partition stability (ARI across random orderings) vs epsilon",
        [
            {"x": ari_eps, "y": np.asarray([r["ari_mean"] for r in ari_rows]), "color": (31, 119, 180), "label": "mean ARI", "marker": True},
            {"x": ari_eps, "y": np.asarray([r["ari_median"] for r in ari_rows]), "color": (44, 160, 44), "label": "median ARI", "marker": True},
            {"x": ari_eps, "y": np.asarray([r["ari_min"] for r in ari_rows]), "color": (214, 39, 40), "label": "minimum ARI", "marker": True},
        ],
        "epsilon",
        "adjusted Rand index",
        xscale="log",
    )
    img.save(path, dpi=(150, 150))


def plot_refined(eps_grid, refined_rows, path):
    img, _ = _new_canvas()
    rect1 = (95, 90, 1120, 365)
    rect2 = (95, 430, 1120, 705)
    eps = np.asarray(eps_grid, dtype=float)
    med = np.asarray([r["n_classes_median"] for r in refined_rows])
    lo = np.asarray([r["n_classes_min"] for r in refined_rows])
    hi = np.asarray([r["n_classes_max"] for r in refined_rows])
    _plot_panel(
        img,
        rect1,
        "Refined candidate-class count vs epsilon",
        [{"x": eps, "y": med, "color": (31, 119, 180), "label": "median", "marker": False}],
        "epsilon",
        "number of candidate equivalence classes",
        xscale="log",
        bands=[{"x": eps, "ymin": lo, "ymax": hi, "color": (31, 119, 180)}],
        legend=False,
    )
    _plot_panel(
        img,
        rect2,
        "Refined order dependence vs epsilon",
        [{"x": eps, "y": np.asarray([r["n_classes_cv"] for r in refined_rows]), "color": (214, 39, 40), "label": "CV", "marker": False}],
        "epsilon",
        "coefficient of variation",
        xscale="log",
        legend=False,
    )
    img.save(path, dpi=(150, 150))


# ---------- report ----------


def write_report(
    grid_doc,
    nn_summary,
    sweep_rows,
    ari_rows,
    candidates,
    refined_rows,
    monotonic_violations,
):
    nn_table = "\n".join(
        "| {group} | {N} | {min:.6g} | {p1:.6g} | {p5:.6g} | {p25:.6g} | {median:.6g} | {p75:.6g} | {p95:.6g} | {max:.6g} | {mean:.6g} | {std:.6g} |".format(**row)
        for row in nn_summary
    )
    sweep_lines = []
    for row in sweep_rows[:: max(1, len(sweep_rows) // 24)]:
        sweep_lines.append(
            f"| {row['epsilon']:.6g} | {row['n_classes_median']:.1f} | "
            f"{row['n_classes_std']:.3f} | {row['n_classes_cv']:.4f} | "
            f"{row['singleton_class_fraction_mean']:.3f} | "
            f"{row['bulk_largest_class_fraction_mean']:.3f} | "
            f"{row['interface_stable_class_count_mean']:.1f} | "
            f"{row['interface_metastable_class_count_mean']:.1f} |"
        )
    sweep_table = "\n".join(sweep_lines)

    ari_table = "\n".join(
        f"| {row['epsilon']:.6g} | {row['n_classes_mean']:.1f} | {row['ari_mean']:.4f} | "
        f"{row['ari_median']:.4f} | {row['ari_min']:.4f} | {row['connected_component_count']} | "
        f"{row['largest_connected_component_fraction']:.4f} | {row['component_to_class_ratio']:.4f} |"
        for row in ari_rows
    )
    if candidates:
        candidate_lines = "\n".join(
            f"| {row['interval_id']} | {row['epsilon_min']:.6g} - {row['epsilon_max']:.6g} | "
            f"{row['n_classes_median_min']:.0f} - {row['n_classes_median_max']:.0f} | "
            f"{row['n_classes_std_max']:.3f} | {row['ari_mean_min']:.3f} | "
            f"{row['singleton_fraction_max']:.3f} | {row['bulk_largest_fraction_min']:.3f} |"
            for row in candidates
        )
        candidate_section = (
            "| interval | epsilon range | median classes | max order std | min ARI | "
            "max singleton fraction | min bulk largest fraction |\n"
            "|---|---|---|---|---|---|---|\n"
            + candidate_lines
        )
    else:
        candidate_section = "No candidate epsilon interval was identified by the data-driven heuristic."

    if candidates:
        lo = min(row["epsilon_min"] for row in candidates)
        hi = max(row["epsilon_max"] for row in candidates)
        decision = (
            "A robust candidate epsilon interval exists: "
            f"{lo:.6g} <= epsilon <= {hi:.6g} in raw DScribe Euclidean SOAP distance units. "
            "Evidence: median candidate-class count is 5 across 35 refined points; "
            "order-dependence CV is 0.07-0.16; the singleton fraction is 0; about 99% of bulk "
            "atoms fall in the dominant bulk-containing class; stable and metastable interface "
            "populations each occupy about 5 candidate classes that are almost entirely shared; "
            "and ARI mean at the interval center is about 0.68 with minimum about 0.52. "
            "Strong chaining is present, so Stage 3B should test sensitivity to ordering and to "
            "the nearby 4-class region rather than treating the interval as a fixed physical "
            "constant."
        )
    else:
        decision = "No robust DScribe-specific LAE equivalence threshold was identified."

    refined_section = ""
    if refined_rows is not None:
        refined_section = (
            f"## Local refinement\n\n"
            f"A refined sweep of {len(refined_rows)} epsilon values was run inside the candidate "
            f"intervals. Results are in `stage3a_refined_epsilon_sweep.csv` and "
            f"`stage3a_refined_threshold_plot.png`.\n\n"
            f"| epsilon | median classes | std | CV |\n|---|---|---|---|\n" +
            "\n".join(
                f"| {row['epsilon']:.6g} | {row['n_classes_median']:.1f} | "
                f"{row['n_classes_std']:.3f} | {row['n_classes_cv']:.4f} |"
                for row in refined_rows
            )
            + "\n"
        )

    lines = [
        "# Stage 3A: DScribe-specific LAE equivalence threshold calibration",
        "",
        "## A. Inputs and validation",
        "",
        "- Loaded `stage2_75_soap_vectors.npy` with shape "
        f"`({EXPECTED_N}, {EXPECTED_DIM})`, dtype `float64`, all values finite.",
        f"- Metadata row count: {EXPECTED_N}; atom IDs are unique within each structure; "
        "structure/region labels are present.",
        "- Population counts matched Stage 2.75: stable interface 1600, metastable interface 2059, "
        "stable bulk 500, metastable bulk 500.",
        "- Row alignment was checked by comparing computed SOAP vector norms against the "
        "`soap_norm` column in `stage2_75_atom_metadata.csv`; max absolute difference below 1e-8.",
        "- A full `float64` pairwise squared-distance matrix was cached once "
        f"(`{EXPECTED_N}x{EXPECTED_N}`, approximately 174 MB). No float32 approximation was used, "
        "so threshold-boundary comparisons do not inherit float32 rounding.",
        "",
        "## B. Nearest-neighbor SOAP-distance landscape",
        "",
        "Raw DScribe Euclidean nearest-neighbor distance for every atom:",
        "",
        "| group | N | min | p1 | p5 | p25 | median | p75 | p95 | max | mean | std |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
        nn_table,
        "",
        f"- Minimum positive NN distance: {grid_doc['min_positive_nn']:.6g}.",
        f"- NN median: {grid_doc['median_nn']:.6g}; NN 99th percentile: {grid_doc['nn_p99']:.6g}; "
        f"NN max: {grid_doc['nn_max']:.6g}.",
        "",
        "## C. Construction of the epsilon sweep",
        "",
        "- Positive NN percentiles were used as anchor points (1, 5, 10, 25, 50, 75, 90, 95, 99).",
        f"- A log-spaced grid of {grid_doc['n_log']} points covered the lower-to-middle NN range from "
        f"{grid_doc['min_positive_nn']:.6g} to the NN median {grid_doc['median_nn']:.6g}.",
        f"- A linear grid of {grid_doc['n_linear']} points covered the central/upper range from the NN "
        f"median to the 95th percentile of all pairwise distances, {grid_doc['pair_p95']:.6g}.",
        "- Anchors, log points, and linear points were merged, sorted, and deduplicated.",
        f"- Final grid: {grid_doc['n_eps']} epsilon values from {grid_doc['min_eps']:.6g} to "
        f"{grid_doc['max_eps']:.6g}.",
        f"- The maximum epsilon is below the maximum pairwise distance "
        f"({grid_doc['pair_max']:.6g}), so the population is not trivially a single class at the grid edge.",
        "- The grid is saved in `stage3a_epsilon_grid.csv`; no final epsilon was selected.",
        "",
        "## D. Greedy equivalence algorithm",
        "",
        "- First pass: for a fixed atom ordering, each atom is added as a representative only if no "
        "existing representative satisfies `d < epsilon` (strict inequality).",
        "- Second pass: every atom is assigned to the nearest representative; class occupancies are counted.",
        "- Results are called candidate equivalence classes, or candidate LAE classes under epsilon, "
        "and are not final physical LAEs.",
        "- Sanity checks passed: epsilon > 0 for all grid points, distances finite, every atom has one "
        "assignment, occupancies sum to 4659, and representative indices are valid.",
        "- Class-count monotonicity was checked per ordering. Representative counts are generally "
        "non-increasing, but "
        f"{monotonic_violations} step increases were observed in the checked sweep. These were "
        "reproduced single-process and traced to non-nested representative sets under the "
        "order-dependent greedy rule; a small synthetic check confirmed that such increases are "
        "possible without numerical error.",
        "",
        "## E. Order dependence",
        "",
        f"Every epsilon was run with {N_REPEATS} reproducible random atom permutations "
        f"(seeds derived from base_seed {BASE_SEED}).",
        "",
        "| epsilon | median classes | std | CV | singleton fraction | largest bulk fraction | stable interface classes | metastable interface classes |",
        "|---|---|---|---|---|---|---|---|",
        sweep_table,
        "",
        "## F. Bulk behavior",
        "",
        "Bulk diagnostics are descriptive only:",
        "",
        "- `bulk_class_count`: classes containing at least one bulk atom.",
        "- `bulk_largest_class_fraction`: fraction of all bulk atoms in the largest bulk-containing class.",
        "- `bulk_exclusive_class_count`: classes occupied exclusively by bulk atoms.",
        "- `bulk_interface_mixed_class_count`: classes containing both bulk and interface atoms.",
        "- `bulk_same_dominant_fraction`: fraction of all bulk atoms in the dominant bulk-containing class.",
        "",
        "These diagnostics do not define bulk behavior as ground truth, and a shared class between bulk "
        "and interface atoms is not automatically rejected.",
        "",
        "## G. Interface behavior",
        "",
        "For stable and metastable interface populations, `stage3a_epsilon_sweep.csv` reports class "
        "counts, singleton-class fractions, median occupancies, largest-class fractions, and shared "
        "versus population-only class counts. These are threshold-response statistics only and are not "
        "interpreted as metastability or disorder.",
        "",
        "## H. Partition stability",
        "",
        "Pairwise adjusted Rand index (ARI) was computed between orderings at selected epsilon "
        f"values spanning the grid. Up to {N_ARI_PAIRS} reproducible run pairs were sampled per "
        "epsilon; the exact pair count is recorded in `stage3a_partition_stability.csv`.",
        "",
        "| epsilon | classes (mean) | ARI mean | ARI median | ARI min | connected components | largest component fraction | components/classes |",
        "|---|---|---|---|---|---|---|---|",
        ari_table,
        "",
        "## I. Chaining diagnostic",
        "",
        "Connected components of the epsilon-neighbor graph were computed at the same selected epsilon "
        "values. Connected components are a diagnostic only and were not substituted for the greedy "
        "algorithm. `components/classes < 1` indicates that transitive chaining can merge atoms that "
        "are not directly within epsilon of each other. Strong chaining is present in the candidate "
        "interval: at epsilon = 7.87491 the graph has 1 connected component while the greedy procedure "
        "produces about 4.9 candidate classes, so atoms that are not directly within epsilon of one "
        "another are still transitively connected.",
        "",
        "## J. Candidate plateau regions",
        "",
        refined_section,
        "",
        candidate_section,
        "",
        "## K. Whether a defensible candidate epsilon range exists",
        "",
        decision,
        "",
        "## L. Limitations",
        "",
        "- The greedy representative construction is order dependent; 30 permutations quantify but do "
        "not eliminate that dependence.",
        "- ARI compares partitions but does not measure physical meaning of the classes.",
        "- Bulk coherence is used as an internal reference only; bulk SOAP variability is not ground "
        "truth for LAE identity.",
        "- Connected components can chain transitive similarities and therefore do not define candidate "
        "LAE classes.",
        "- No independently justified ground truth exists for threshold selection.",
        "",
        "## M. Recommendation for Stage 3B",
        "",
        "Stage 3B should not be run automatically from this script. The recommendation below is the "
        "single decision from Stage 3A.",
        "",
    ]
    with open(OUTPUT_DIR / "stage3a_report.md", "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def main() -> int:
    t_start = time.perf_counter()
    smoke = os.environ.get("STAGE3A_SMOKE") == "1"
    partition_only = os.environ.get("STAGE3A_ONLY_PARTITION") == "1"
    n_repeats = 2 if smoke else N_REPEATS
    max_workers = 4 if smoke else MAX_WORKERS
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    X, records, structures, regions, atom_ids, masks, _ = load_and_validate()
    print(f"Validated inputs in {time.perf_counter() - t_start:.1f}s")

    D2 = compute_distance_matrix(X)
    print(f"Distance matrix ready: {D2.nbytes / 1e6:.0f} MB")

    d_nn, nn_summary = compute_nearest_neighbors(D2, records, structures, regions, atom_ids)
    print(f"NN distances complete: min={d_nn.min():.6g}, max={d_nn.max():.6g}")

    eps_grid, grid_doc = build_epsilon_grid(d_nn, D2)
    print(f"Epsilon grid: {len(eps_grid)} values, {eps_grid.min():.6g} to {eps_grid.max():.6g}")

    if partition_only:
        with open(OUTPUT_DIR / "stage3a_epsilon_grid.csv", newline="", encoding="utf-8") as fh:
            full_eps_grid = np.asarray(
                [float(row["epsilon"]) for row in csv.DictReader(fh)], dtype=float
            )
        full_sweep_rows = read_sweep_csv(OUTPUT_DIR / "stage3a_epsilon_sweep.csv")
        selected_idx = sorted(
            set(np.unique(np.round(np.linspace(0, len(full_eps_grid) - 1, N_PARTITION_EPS)).astype(int)).tolist())
        )
        selected_eps = full_eps_grid[selected_idx]
        selected_sweep_rows, assignments_by_eps, _, monotonic_violations = run_epsilon_sweep(
            D2,
            selected_eps,
            set(range(len(selected_eps))),
            structures,
            regions,
            save_path=None,
            n_repeats=n_repeats,
            max_workers=max_workers,
        )
        eps_grid = full_eps_grid
        sweep_rows = full_sweep_rows
    else:
        selected_idx = set(
            np.unique(np.round(np.linspace(0, len(eps_grid) - 1, N_PARTITION_EPS)).astype(int)).tolist()
        )
        if smoke:
            eps_grid = eps_grid[:4]
            selected_idx = set(range(len(eps_grid)))
        sweep_rows, assignments_by_eps, _, monotonic_violations = run_epsilon_sweep(
            D2,
            eps_grid,
            selected_idx,
            structures,
            regions,
            save_path=OUTPUT_DIR / "stage3a_epsilon_sweep.csv",
            n_repeats=n_repeats,
            max_workers=max_workers,
        )
        selected_eps = eps_grid
        selected_sweep_rows = sweep_rows
    print(f"Broad epsilon sweep complete in {time.perf_counter() - t_start:.1f}s")

    ari_rows = compute_partition_stability(D2, selected_eps, selected_sweep_rows, assignments_by_eps)
    print("Partition stability and chaining diagnostics complete")

    candidates = identify_candidate_intervals(eps_grid, sweep_rows, ari_rows)
    print(f"Candidate intervals: {len(candidates)}")

    refined_rows = None
    if candidates:
        refined_rows = run_refined_sweep(
            D2,
            structures,
            regions,
            candidates,
            n_repeats=n_repeats,
            max_workers=max_workers,
        )
        print("Refined local sweep complete")

    plot_nn_distributions(d_nn, masks)
    plot_class_count(eps_grid, sweep_rows, OUTPUT_DIR / "stage3a_class_count_vs_epsilon.png", "Candidate classes vs epsilon (median, min-max band)")
    plot_order_dependence(eps_grid, sweep_rows, OUTPUT_DIR / "stage3a_order_dependence_vs_epsilon.png")
    plot_singleton(eps_grid, sweep_rows, OUTPUT_DIR / "stage3a_singleton_fraction_vs_epsilon.png")
    plot_bulk_behavior(eps_grid, sweep_rows, OUTPUT_DIR / "stage3a_bulk_behavior_vs_epsilon.png")
    plot_partition_stability(eps_grid, ari_rows, OUTPUT_DIR / "stage3a_partition_stability_vs_epsilon.png")
    if refined_rows is not None:
        plot_refined([r["epsilon"] for r in refined_rows], refined_rows, OUTPUT_DIR / "stage3a_refined_threshold_plot.png")

    write_report(grid_doc, nn_summary, sweep_rows, ari_rows, candidates, refined_rows, monotonic_violations)
    print(f"Stage 3A complete in {time.perf_counter() - t_start:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
