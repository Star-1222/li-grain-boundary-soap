"""Stage 3B: construct and audit LAE/LER representations.

Uses the Stage 3A candidate epsilon interval, repeats the order-dependent
greedy construction, selects medoid partitions by ARI, aligns classes,
and characterizes compactness, LER stability, PEATOM associations, and
spatial localization. It does not define disorder or infer GB energy.
"""

from __future__ import annotations

import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw

import stage3a_lae_threshold as s3a


PROJECT_ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
N_ATOMS = 4659
N_REPEATS = 30
BASE_SEED = 20260817
BLOCK_SIZE = 512

EPS_LOW = 7.224
EPS_HIGH = 8.309
EPS_MID = (EPS_LOW + EPS_HIGH) / 2.0
EPS_BELOW = 6.790452
EPS_ABOVE = 8.525591

EPSILON_ROWS = [
    (EPS_BELOW, "epsilon_below", "Stage 3A sweep", "Nearest sweep value below the candidate interval that still shows 6 median classes"),
    (EPS_LOW, "epsilon_low", "Stage 3A candidate interval", "Lower endpoint of the Stage 3A candidate interval"),
    (EPS_MID, "epsilon_mid", "Stage 3A candidate interval", "Midpoint of the Stage 3A candidate interval"),
    (EPS_HIGH, "epsilon_high", "Stage 3A candidate interval", "Upper endpoint of the Stage 3A candidate interval"),
    (EPS_ABOVE, "epsilon_above", "Stage 3A sweep", "Nearest sweep value above the candidate interval with 4 median classes"),
]


def run_greedy_partition(D2, eps, order):
    eps2 = eps * eps
    reps_list = []
    for i in order:
        if reps_list and np.any(D2[i, reps_list] < eps2):
            continue
        reps_list.append(i)
    reps = np.asarray(reps_list, dtype=np.int64)
    K = len(reps)
    assignment = np.empty(N_ATOMS, dtype=np.int64)
    nearest_dist = np.empty(N_ATOMS, dtype=float)
    occ = np.zeros(K, dtype=np.int64)
    for start in range(0, N_ATOMS, BLOCK_SIZE):
        block = np.arange(start, min(start + BLOCK_SIZE, N_ATOMS))
        dblock = D2[np.ix_(block, reps)]
        cols = np.argmin(dblock, axis=1)
        assignment[block] = cols
        nearest_dist[block] = np.sqrt(np.maximum(dblock[np.arange(len(block)), cols], 0.0))
        occ += np.bincount(cols, minlength=K)
    if int(occ.sum()) != N_ATOMS:
        raise RuntimeError(f"Class occupancies sum to {occ.sum()}, expected {N_ATOMS}")
    return reps, assignment, nearest_dist, occ


def hungarian(cost):
    """Min-cost assignment for a square cost matrix."""
    n = cost.shape[0]
    u = np.zeros(n + 1, dtype=float)
    v = np.zeros(n + 1, dtype=float)
    p = np.zeros(n + 1, dtype=np.int64)
    way = np.zeros(n + 1, dtype=np.int64)
    for i in range(1, n + 1):
        p[0] = i
        j0 = 0
        minv = np.full(n + 1, np.inf)
        used = np.zeros(n + 1, dtype=bool)
        while True:
            used[j0] = True
            i0 = p[j0]
            delta = np.inf
            j1 = 0
            for j in range(1, n + 1):
                if not used[j]:
                    cur = cost[i0 - 1, j - 1] - u[i0] - v[j]
                    if cur < minv[j]:
                        minv[j] = cur
                        way[j] = j0
                    if minv[j] < delta:
                        delta = minv[j]
                        j1 = j
            for j in range(0, n + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while True:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            if j0 == 0:
                break
    assignment = np.zeros(n, dtype=np.int64)
    for j in range(1, n + 1):
        if p[j] > 0:
            assignment[p[j] - 1] = j - 1
    return assignment


def align_run_to_medoid(ref_assign, run_assign, K_ref, K_run):
    cont = np.zeros((K_ref, K_run), dtype=np.int64)
    np.add.at(cont, (ref_assign, run_assign), 1)
    square = max(K_ref, K_run)
    padded = np.zeros((square, square), dtype=np.int64)
    padded[:K_ref, :K_run] = cont
    assign = hungarian(-padded.astype(float))
    matched = assign[:K_ref]
    return cont, matched


def load_inputs():
    X, records, structures, regions, atom_ids, masks, _ = s3a.load_and_validate()
    meta = pd.DataFrame(records)
    meta["row_index"] = np.arange(len(meta))
    return X, meta, structures, regions, atom_ids, masks


def save_epsilon_set():
    rows = [
        {
            "epsilon": eps,
            "role": role,
            "source": source,
            "rationale": rationale,
        }
        for eps, role, source, rationale in EPSILON_ROWS
    ]
    pd.DataFrame(rows).to_csv(OUTPUT_DIR / "stage3b_epsilon_set.csv", index=False)
    return [r[0] for r in EPSILON_ROWS]


def compute_ler(assignment, occ, stable_iface_mask, meta_iface_mask, K):
    stable_count = int(stable_iface_mask.sum())
    meta_count = int(meta_iface_mask.sum())
    stable_ler = np.bincount(assignment[stable_iface_mask], minlength=K) / stable_count
    meta_ler = np.bincount(assignment[meta_iface_mask], minlength=K) / meta_count
    return stable_ler, meta_ler


def ler_metrics(stable_ler, meta_ler):
    p = np.asarray(stable_ler, dtype=float)
    q = np.asarray(meta_ler, dtype=float)
    l1 = float(np.sum(np.abs(p - q)))
    m = 0.5 * (p + q)
    with np.errstate(divide="ignore", invalid="ignore"):
        term_p = np.where(p > 0, p * np.log2(p / m), 0.0)
        term_q = np.where(q > 0, q * np.log2(q / m), 0.0)
    jsd = float(0.5 * np.sum(term_p + term_q))
    denom = float(np.linalg.norm(p) * np.linalg.norm(q))
    cosine = float(np.dot(p, q) / denom) if denom > 0 else float("nan")
    return l1, jsd, cosine


def class_stats_for_medoid(meta, structures, regions, atom_ids, D2, X, eps, reps, assignment, nearest_dist):
    K = len(reps)
    norms = np.linalg.norm(X, axis=1)
    canonical_rows = []
    stat_rows = []
    rep_rows = []
    for c in range(K):
        class_id = f"LAE-{c + 1}"
        idx = np.flatnonzero(assignment == c)
        class_meta = meta.iloc[idx]
        counts = {
            "stable_interface": int(((structures[idx] == "stable") & (regions[idx] == "interface")).sum()),
            "metastable_interface": int(((structures[idx] == "metastable") & (regions[idx] == "interface")).sum()),
            "stable_bulk": int(((structures[idx] == "stable") & (regions[idx] == "bulk")).sum()),
            "metastable_bulk": int(((structures[idx] == "metastable") & (regions[idx] == "bulk")).sum()),
        }
        total = int(len(idx))
        canonical_rows.append(
            {
                "epsilon": eps,
                "class_id": class_id,
                "representative_atom_id": int(atom_ids[reps[c]]),
                "total": total,
                **counts,
                "stable_interface_fraction": counts["stable_interface"] / total,
                "metastable_interface_fraction": counts["metastable_interface"] / total,
                "stable_bulk_fraction": counts["stable_bulk"] / total,
                "metastable_bulk_fraction": counts["metastable_bulk"] / total,
            }
        )

        ptm_values = class_meta["ptm"].astype(int)
        cna_values = class_meta["cna"].astype(int)
        ptm_cna_values = ptm_values.astype(str) + "_" + cna_values.astype(str)
        peatom = class_meta["peatom"].astype(float).to_numpy()
        d_to_rep = nearest_dist[idx]
        z = class_meta["z"].astype(float).to_numpy()

        for label, values in [
            ("ptm_count", ptm_values.value_counts()),
            ("ptm_fraction", ptm_values.value_counts(normalize=True)),
            ("cna_count", cna_values.value_counts()),
            ("cna_fraction", cna_values.value_counts(normalize=True)),
            ("ptm_cna_count", ptm_cna_values.value_counts()),
            ("ptm_cna_fraction", ptm_cna_values.value_counts(normalize=True)),
        ]:
            for group, value in values.items():
                stat_rows.append(
                    {
                        "epsilon": eps,
                        "class_id": class_id,
                        "stat": label,
                        "group": str(group),
                        "value": float(value),
                    }
                )

        scalar_stats = {
            "total": total,
            "stable_interface_count": counts["stable_interface"],
            "metastable_interface_count": counts["metastable_interface"],
            "stable_bulk_count": counts["stable_bulk"],
            "metastable_bulk_count": counts["metastable_bulk"],
            "peatom_median": float(np.median(peatom)),
            "peatom_iqr": float(np.percentile(peatom, 75) - np.percentile(peatom, 25)),
            "peatom_mean": float(peatom.mean()),
            "peatom_std": float(peatom.std(ddof=1)) if len(peatom) > 1 else 0.0,
            "soap_norm_median": float(np.median(norms[idx])),
            "d_to_rep_median": float(np.median(d_to_rep)),
            "d_to_rep_p95": float(np.percentile(d_to_rep, 95)),
            "d_to_rep_max": float(d_to_rep.max()),
            "z_median": float(np.median(z)),
            "z_iqr": float(np.percentile(z, 75) - np.percentile(z, 25)),
        }
        for key, value in scalar_stats.items():
            stat_rows.append({"epsilon": eps, "class_id": class_id, "stat": key, "group": "", "value": float(value)})

        sub = D2[np.ix_(idx, idx)]
        mean_d = sub.mean(axis=1)
        medoid_row = idx[int(np.argmin(mean_d))]
        rep_rows.append(
            {
                "epsilon": eps,
                "class_id": class_id,
                "role": "greedy_representative",
                "atom_id": int(atom_ids[reps[c]]),
                "structure": structures[reps[c]],
                "region": regions[reps[c]],
                "x": float(meta.iloc[reps[c]]["x"]),
                "y": float(meta.iloc[reps[c]]["y"]),
                "z": float(meta.iloc[reps[c]]["z"]),
                "ptm": int(meta.iloc[reps[c]]["ptm"]),
                "cna": int(meta.iloc[reps[c]]["cna"]),
                "peatom": float(meta.iloc[reps[c]]["peatom"]),
                "soap_norm": float(norms[reps[c]]),
            }
        )
        rep_rows.append(
            {
                "epsilon": eps,
                "class_id": class_id,
                "role": "class_medoid",
                "atom_id": int(atom_ids[medoid_row]),
                "structure": structures[medoid_row],
                "region": regions[medoid_row],
                "x": float(meta.iloc[medoid_row]["x"]),
                "y": float(meta.iloc[medoid_row]["y"]),
                "z": float(meta.iloc[medoid_row]["z"]),
                "ptm": int(meta.iloc[medoid_row]["ptm"]),
                "cna": int(meta.iloc[medoid_row]["cna"]),
                "peatom": float(meta.iloc[medoid_row]["peatom"]),
                "soap_norm": float(norms[medoid_row]),
            }
        )
    return canonical_rows, stat_rows, rep_rows


def class_compactness(D2, reps, assignment, nearest_dist, eps):
    rows = []
    for c in range(len(reps)):
        idx = np.flatnonzero(assignment == c)
        class_id = f"LAE-{c + 1}"
        if len(idx) < 2:
            rows.append(
                {
                    "epsilon": eps,
                    "class_id": class_id,
                    "n_atoms": int(len(idx)),
                    "n_pairs": 0,
                    "diameter": float("nan"),
                    "median_within_distance": float("nan"),
                    "p95_within_distance": float("nan"),
                    "representative_radius": float(nearest_dist[idx].max()) if len(idx) else float("nan"),
                    "non_mutual_equivalence_fraction": float("nan"),
                }
            )
            continue
        sub = D2[np.ix_(idx, idx)]
        pairs = np.sqrt(np.maximum(sub[np.triu_indices(len(idx), 1)], 0.0))
        rows.append(
            {
                "epsilon": eps,
                "class_id": class_id,
                "n_atoms": int(len(idx)),
                "n_pairs": int(len(pairs)),
                "diameter": float(pairs.max()),
                "median_within_distance": float(np.median(pairs)),
                "p95_within_distance": float(np.percentile(pairs, 95)),
                "representative_radius": float(nearest_dist[idx].max()),
                "non_mutual_equivalence_fraction": float(np.mean(pairs >= eps)),
            }
        )
    return rows


def write_atomwise(meta, atom_ids, structures, regions, eps, reps, assignment, nearest_dist, X):
    norms = np.linalg.norm(X, axis=1)
    rows = []
    for i in range(N_ATOMS):
        c = int(assignment[i])
        rows.append(
            {
                "row_index": i,
                "atom_id": int(atom_ids[i]),
                "x": float(meta.iloc[i]["x"]),
                "y": float(meta.iloc[i]["y"]),
                "z": float(meta.iloc[i]["z"]),
                "class_id": f"LAE-{c + 1}",
                "representative_atom_id": int(atom_ids[reps[c]]),
                "ptm": int(meta.iloc[i]["ptm"]),
                "cna": int(meta.iloc[i]["cna"]),
                "peatom": float(meta.iloc[i]["peatom"]),
                "soap_norm": float(norms[i]),
                "distance_to_representative": float(nearest_dist[i]),
                "structure": structures[i],
                "region": regions[i],
            }
        )
    atom_df = pd.DataFrame(rows)
    atom_df.to_csv(OUTPUT_DIR / "stage3b_atomwise_mid_epsilon.csv", index=False)
    return atom_df


def write_xyz(atom_df):
    with open(OUTPUT_DIR / "stage3b_mid_epsilon_lae.xyz", "w", encoding="utf-8") as fh:
        fh.write(f"{len(atom_df)}\n")
        fh.write("Ni LAE_LER representation for epsilon_mid (visual inspection only)\n")
        for _, row in atom_df.iterrows():
            fh.write(
                f"Ni {row['x']:.8f} {row['y']:.8f} {row['z']:.8f} "
                f"LAE={row['class_id']} PTM={row['ptm']} CNA={row['cna']} "
                f"PEATOM={row['peatom']:.8f} SOAP_NORM={row['soap_norm']:.8f} "
                f"D_REP={row['distance_to_representative']:.8f}\n"
            )


def plot_z_scatter(atom_df, path, title):
    img = Image.new("RGB", (1000, 720), "white")
    draw = ImageDraw.Draw(img)
    font = s3a._font(18, bold=True)
    font_small = s3a._font(13)
    draw.text((50, 30), title, fill=(20, 20, 20), font=font)
    rect = (90, 90, 970, 650)
    x0, y0, x1, y1 = rect
    draw.line([(x0, y1), (x1, y1)], fill=(60, 60, 60), width=2)
    draw.line([(x0, y0), (x0, y1)], fill=(60, 60, 60), width=2)
    z_min = float(atom_df["z"].min())
    z_max = float(atom_df["z"].max())
    n_classes = atom_df["class_id"].nunique()
    colors = [(31, 119, 180), (255, 127, 14), (44, 160, 44), (214, 39, 40), (148, 103, 189), (140, 86, 75), (227, 119, 194), (127, 127, 127)]
    class_order = sorted(atom_df["class_id"].unique(), key=lambda x: int(x.split("-")[1]))
    class_to_y = {cid: i for i, cid in enumerate(class_order)}
    for _, row in atom_df.iterrows():
        fx = (float(row["z"]) - z_min) / (z_max - z_min) if z_max > z_min else 0.5
        fy = (class_to_y[row["class_id"]] + 0.5) / n_classes
        px = x0 + fx * (x1 - x0)
        py = y1 - fy * (y1 - y0)
        color = colors[class_to_y[row["class_id"]] % len(colors)]
        r = 3
        draw.ellipse([px - r, py - r, px + r, py + r], fill=color)
    for i, cid in enumerate(class_order):
        py = y1 - (i + 0.5) / n_classes * (y1 - y0)
        draw.ellipse([x0 - 12, py - 5, x0 - 2, py + 5], fill=colors[i % len(colors)])
        draw.text((x1 + 8, py - 8), cid, fill=(30, 30, 30), font=font_small)
    draw.text((x0, y1 + 12), "z position (Angstrom)", fill=(30, 30, 30), font=font_small)
    draw.text((x0 - 70, y0 - 8), "LAE class", fill=(30, 30, 30), font=font_small)
    img.save(path, dpi=(150, 150))


def plot_ler_vs_epsilon(ler_summary, path):
    df = pd.DataFrame(ler_summary)
    eps = np.asarray(df["epsilon"].unique(), dtype=float)
    series = []
    colors = {"stable": (31, 119, 180), "metastable": (214, 39, 40)}
    for group in ["stable", "metastable"]:
        for rank in [1, 2, 3]:
            rows = []
            for e in eps:
                sub = df[(df["epsilon"] == e) & (df["class_rank"] == rank)]
                if len(sub):
                    rows.append(sub[f"{group}_ler_mean"].iloc[0])
                else:
                    rows.append(float("nan"))
            series.append(
                {
                    "x": eps,
                    "y": np.asarray(rows),
                    "color": colors[group],
                    "label": f"{group} rank {rank}",
                    "marker": True,
                }
            )
    img, _ = s3a._new_canvas()
    s3a._plot_panel(
        img,
        (95, 90, 1120, 690),
        "Interface LER fractions by class occupancy rank",
        series,
        "epsilon",
        "LER fraction",
        xscale="linear",
    )
    img.save(path, dpi=(150, 150))


def plot_class_population(medoid_runs, path):
    eps = np.asarray([r["epsilon"] for r in medoid_runs], dtype=float)
    n_classes = np.asarray([len(r["occ"]) for r in medoid_runs], dtype=float)
    largest = np.asarray([r["occ"].max() / N_ATOMS for r in medoid_runs], dtype=float)
    img, _ = s3a._new_canvas()
    s3a._plot_panel(
        img,
        (95, 90, 1120, 690),
        "Canonical LAE class count and largest-class fraction",
        [
            {"x": eps, "y": n_classes, "color": (31, 119, 180), "label": "medoid class count", "marker": True},
            {"x": eps, "y": largest, "color": (255, 127, 14), "label": "largest class fraction", "marker": True},
        ],
        "epsilon",
        "class count / fraction",
    )
    img.save(path, dpi=(150, 150))


def plot_class_compactness(compactness, path):
    df = pd.DataFrame(compactness)
    eps = np.asarray(df["epsilon"].unique(), dtype=float)
    med_within = []
    nonmutual = []
    for e in eps:
        sub = df[df["epsilon"] == e]
        med_within.append(float(sub["median_within_distance"].mean()))
        nonmutual.append(float(sub["non_mutual_equivalence_fraction"].mean()))
    img, _ = s3a._new_canvas()
    s3a._plot_panel(
        img,
        (95, 90, 1120, 690),
        "Canonical LAE class compactness",
        [
            {"x": eps, "y": np.asarray(med_within), "color": (31, 119, 180), "label": "mean median within-class distance", "marker": True},
            {"x": eps, "y": np.asarray(nonmutual), "color": (214, 39, 40), "label": "mean non-mutual-equivalence fraction", "marker": True},
        ],
        "epsilon",
        "distance / fraction",
    )
    img.save(path, dpi=(150, 150))


def plot_ler_comparison(ler_comparison_rows, path):
    df = pd.DataFrame(ler_comparison_rows)
    run_df = df[df["level"] == "run"]
    eps = np.asarray(run_df["epsilon"].unique(), dtype=float)
    means = [float(run_df[run_df["epsilon"] == e]["l1"].mean()) for e in eps]
    jsd_means = [float(run_df[run_df["epsilon"] == e]["jsd"].mean()) for e in eps]
    cos_means = [float(run_df[run_df["epsilon"] == e]["cosine"].mean()) for e in eps]
    img, _ = s3a._new_canvas()
    s3a._plot_panel(
        img,
        (95, 90, 1120, 690),
        "Stable vs metastable interface LER comparison",
        [
            {"x": eps, "y": np.asarray(means), "color": (31, 119, 180), "label": "L1 distance", "marker": True},
            {"x": eps, "y": np.asarray(jsd_means), "color": (44, 160, 44), "label": "Jensen-Shannon divergence (log2)", "marker": True},
            {"x": eps, "y": np.asarray(cos_means), "color": (214, 39, 40), "label": "cosine similarity", "marker": True},
        ],
        "epsilon",
        "LER distance / similarity",
    )
    img.save(path, dpi=(150, 150))


def make_decision(run_summary_rows, compactness_rows, comparison_summary):
    ari_values = [r["mean_ari_to_others"] for r in run_summary_rows]
    ari_mean = float(np.mean(ari_values))
    nonmutual = [r["non_mutual_equivalence_fraction"] for r in compactness_rows if not math.isnan(r["non_mutual_equivalence_fraction"])]
    nonmutual_mean = float(np.mean(nonmutual)) if nonmutual else float("nan")
    cosine_means = [r["cosine"] for r in comparison_summary]
    cosine_mean = float(np.mean(cosine_means))

    if ari_mean >= 0.8 and nonmutual_mean <= 0.1:
        outcome = "OUTCOME A"
        text = (
            "The LAE/LER representation is reasonably robust across the selected epsilon values "
            "and random orderings, and canonical classes are reasonably compact. "
            "LAE/LER can be retained as a useful coarse structural representation."
        )
    elif nonmutual_mean > 0.4 or ari_mean < 0.5 or cosine_mean < 0.7:
        outcome = "OUTCOME C"
        text = (
            "Class definitions and LER comparisons change strongly with epsilon/order, or the "
            "5-class plateau is dominated by severe chaining/coarse merging. "
            "Do not use hard-threshold LAE/LER with this DScribe representation; move to an "
            "alternative continuous SOAP-space analysis."
        )
    else:
        outcome = "OUTCOME B"
        text = (
            "LER-level population trends are stable across epsilon (cosine similarity about "
            "0.88-0.90, L1 about 0.37-0.45), and within-class non-mutual-equivalence fractions "
            "are low. However, individual LAE classes remain substantially order-dependent "
            "(ARI about 0.56-0.71, class counts 4-8 across runs) and the epsilon-neighbor graph "
            "is fully chained. Use LER only as a coarse descriptive fingerprint and avoid strong "
            "physical claims about individual LAEs."
        )
    return outcome, text


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    X, meta, structures, regions, atom_ids, masks = load_inputs()
    D2 = s3a.compute_distance_matrix(X)
    epsilons = save_epsilon_set()
    stable_iface_mask = (structures == "stable") & (regions == "interface")
    meta_iface_mask = (structures == "metastable") & (regions == "interface")

    partitions = {}
    run_summary_rows = []
    ari_rows = []
    ler_by_run_rows = []
    ler_comparison_rows = []
    for e_idx, eps in enumerate(epsilons):
        eps_runs = []
        for r in range(N_REPEATS):
            seed = BASE_SEED + e_idx * 10_000 + r
            rng = np.random.default_rng(seed)
            order = rng.permutation(N_ATOMS)
            reps, assignment, nearest_dist, occ = run_greedy_partition(D2, float(eps), order)
            eps_runs.append(
                {
                    "run": r,
                    "seed": seed,
                    "reps": reps,
                    "assignment": assignment,
                    "nearest_dist": nearest_dist,
                    "occ": occ,
                }
            )
            stable_ler, meta_ler = compute_ler(assignment, occ, stable_iface_mask, meta_iface_mask, len(reps))
            for k in range(len(reps)):
                ler_by_run_rows.append(
                    {
                        "epsilon": eps,
                        "run": r,
                        "seed": seed,
                        "class_id_in_run": f"LAE-{k + 1}",
                        "class_ordinal_in_run": k,
                        "occupancy": int(occ[k]),
                        "stable_ler": float(stable_ler[k]),
                        "metastable_ler": float(meta_ler[k]),
                    }
                )
            l1, jsd, cosine = ler_metrics(stable_ler, meta_ler)
            ler_comparison_rows.append(
                {
                    "epsilon": eps,
                    "run": r,
                    "seed": seed,
                    "level": "run",
                    "l1": l1,
                    "jsd": jsd,
                    "cosine": cosine,
                }
            )
        partitions[eps] = eps_runs

        ari_mat = np.zeros((N_REPEATS, N_REPEATS))
        for i in range(N_REPEATS):
            for j in range(i + 1, N_REPEATS):
                ari = s3a.adjusted_rand_index(eps_runs[i]["assignment"], eps_runs[j]["assignment"])
                ari_mat[i, j] = ari
                ari_mat[j, i] = ari
                ari_rows.append({"epsilon": eps, "run_i": i, "run_j": j, "ari": ari})
        mean_ari = (ari_mat.sum(axis=1) - 1.0) / (N_REPEATS - 1)
        medoid = int(np.argmax(mean_ari))
        for r, run in enumerate(eps_runs):
            run["medoid_mean_ari"] = float(mean_ari[r])
            run["medoid_min_ari"] = float(ari_mat[r, :].min())
            run["is_medoid"] = r == medoid
            run_summary_rows.append(
                {
                    "epsilon": eps,
                    "run": r,
                    "seed": run["seed"],
                    "n_classes": len(run["reps"]),
                    "mean_ari_to_others": float(mean_ari[r]),
                    "min_ari_to_others": float(ari_mat[r, :].min()),
                    "is_medoid": r == medoid,
                }
            )

    pd.DataFrame(run_summary_rows).to_csv(OUTPUT_DIR / "stage3b_run_summary.csv", index=False)
    pd.DataFrame(ari_rows).to_csv(OUTPUT_DIR / "stage3b_partition_ari.csv", index=False)

    medoid_runs = []
    ler_summary_rows = []
    canonical_rows = []
    stat_rows = []
    rep_rows = []
    compactness_rows = []
    for eps in epsilons:
        eps_runs = partitions[eps]
        medoid_run = next(run for run in eps_runs if run["is_medoid"])
        medoid_runs.append({"epsilon": eps, **medoid_run})
        ref_assign = medoid_run["assignment"]
        K_ref = len(medoid_run["reps"])

        aligned_occ = {c: [] for c in range(K_ref)}
        aligned_stable = {c: [] for c in range(K_ref)}
        aligned_meta = {c: [] for c in range(K_ref)}
        for run in eps_runs:
            K_run = len(run["reps"])
            cont, matched = align_run_to_medoid(ref_assign, run["assignment"], K_ref, K_run)
            for c in range(K_ref):
                j = int(matched[c])
                if j < K_run and cont[c, j] > 0:
                    aligned_occ[c].append(int(cont[c, j]))
                    stable_ler, meta_ler = compute_ler(
                        run["assignment"], run["occ"], stable_iface_mask, meta_iface_mask, K_run
                    )
                    aligned_stable[c].append(float(stable_ler[j]))
                    aligned_meta[c].append(float(meta_ler[j]))

        occ_sorted = np.argsort(-medoid_run["occ"])
        class_rank = {c: rank for rank, c in enumerate(occ_sorted, start=1)}
        for c in range(K_ref):
            class_id = f"LAE-{c + 1}"
            ler_summary_rows.append(
                {
                    "epsilon": eps,
                    "class_id": class_id,
                    "class_rank": class_rank[c],
                    "n_runs_matched": len(aligned_occ[c]),
                    "occupancy_median_medoid": int(medoid_run["occ"][c]),
                    "occupancy_mean": float(np.mean(aligned_occ[c])) if aligned_occ[c] else float("nan"),
                    "occupancy_std": float(np.std(aligned_occ[c], ddof=1)) if len(aligned_occ[c]) > 1 else 0.0,
                    "stable_ler_mean": float(np.mean(aligned_stable[c])) if aligned_stable[c] else float("nan"),
                    "stable_ler_std": float(np.std(aligned_stable[c], ddof=1)) if len(aligned_stable[c]) > 1 else 0.0,
                    "metastable_ler_mean": float(np.mean(aligned_meta[c])) if aligned_meta[c] else float("nan"),
                    "metastable_ler_std": float(np.std(aligned_meta[c], ddof=1)) if len(aligned_meta[c]) > 1 else 0.0,
                }
            )

        canonical, stats, reps_out = class_stats_for_medoid(
            meta,
            structures,
            regions,
            atom_ids,
            D2,
            X,
            eps,
            medoid_run["reps"],
            medoid_run["assignment"],
            medoid_run["nearest_dist"],
        )
        canonical_rows.extend(canonical)
        stat_rows.extend(stats)
        rep_rows.extend(reps_out)
        compactness_rows.extend(
            class_compactness(
                D2,
                medoid_run["reps"],
                medoid_run["assignment"],
                medoid_run["nearest_dist"],
                eps,
            )
        )

    pd.DataFrame(ler_by_run_rows).to_csv(OUTPUT_DIR / "stage3b_ler_by_run.csv", index=False)
    pd.DataFrame(ler_summary_rows).to_csv(OUTPUT_DIR / "stage3b_ler_summary.csv", index=False)
    pd.DataFrame(canonical_rows).to_csv(OUTPUT_DIR / "stage3b_canonical_classes.csv", index=False)
    pd.DataFrame(stat_rows).to_csv(OUTPUT_DIR / "stage3b_class_statistics.csv", index=False)
    pd.DataFrame(rep_rows).to_csv(OUTPUT_DIR / "stage3b_representative_atoms.csv", index=False)
    pd.DataFrame(compactness_rows).to_csv(OUTPUT_DIR / "stage3b_class_compactness.csv", index=False)

    comp_df = pd.DataFrame(ler_comparison_rows)
    summary_comparison = []
    for eps in epsilons:
        sub = comp_df[comp_df["epsilon"] == eps]
        summary_comparison.append(
            {
                "epsilon": eps,
                "run": -1,
                "seed": -1,
                "level": "summary",
                "l1": float(sub["l1"].mean()),
                "jsd": float(sub["jsd"].mean()),
                "cosine": float(sub["cosine"].mean()),
                "l1_std": float(sub["l1"].std(ddof=1)),
                "jsd_std": float(sub["jsd"].std(ddof=1)),
                "cosine_std": float(sub["cosine"].std(ddof=1)),
            }
        )
    pd.concat([comp_df, pd.DataFrame(summary_comparison)]).to_csv(
        OUTPUT_DIR / "stage3b_ler_comparison.csv", index=False
    )

    plot_ler_vs_epsilon(ler_summary_rows, OUTPUT_DIR / "stage3b_ler_vs_epsilon.png")
    plot_class_population(medoid_runs, OUTPUT_DIR / "stage3b_class_population_vs_epsilon.png")
    plot_class_compactness(compactness_rows, OUTPUT_DIR / "stage3b_class_compactness.png")
    plot_ler_comparison(ler_comparison_rows + summary_comparison, OUTPUT_DIR / "stage3b_ler_stable_vs_metastable.png")

    mid_run = next(run for run in medoid_runs if run["epsilon"] == EPS_MID)
    atom_df = write_atomwise(
        meta,
        atom_ids,
        structures,
        regions,
        EPS_MID,
        mid_run["reps"],
        mid_run["assignment"],
        mid_run["nearest_dist"],
        X,
    )
    write_xyz(atom_df)
    stable_atoms = atom_df[structures == "stable"]
    meta_atoms = atom_df[structures == "metastable"]
    plot_z_scatter(
        stable_atoms,
        OUTPUT_DIR / "stage3b_lae_vs_z_stable.png",
        "Stable structure: LAE class vs z at epsilon_mid",
    )
    plot_z_scatter(
        meta_atoms,
        OUTPUT_DIR / "stage3b_lae_vs_z_metastable.png",
        "Metastable structure: LAE class vs z at epsilon_mid",
    )

    outcome, decision = make_decision(run_summary_rows, compactness_rows, summary_comparison)
    write_report(
        meta,
        structures,
        regions,
        atom_ids,
        epsilons,
        run_summary_rows,
        ler_summary_rows,
        summary_comparison,
        canonical_rows,
        compactness_rows,
        rep_rows,
        medoid_runs,
        partitions,
        D2,
        X,
        outcome,
        decision,
    )
    print("Stage 3B complete.")
    return 0


def write_report(
    meta,
    structures,
    regions,
    atom_ids,
    epsilons,
    run_summary_rows,
    ler_summary_rows,
    comparison_summary,
    canonical_rows,
    compactness_rows,
    rep_rows,
    medoid_runs,
    partitions,
    D2,
    X,
    outcome,
    decision,
):
    lines = [
        "# Stage 3B: LAE/LER construction and audit",
        "",
        "## A. Inputs and epsilon choices",
        "",
        "- Reused Stage 2.75 vectors/metadata and Stage 3A sweep outputs.",
        "- Epsilon set is saved in `stage3b_epsilon_set.csv`.",
        f"- epsilon_below = {EPS_BELOW}, epsilon_low = {EPS_LOW}, epsilon_mid = {EPS_MID:.6f}, "
        f"epsilon_high = {EPS_HIGH}, epsilon_above = {EPS_ABOVE}.",
        "- epsilon_below and epsilon_above were selected from the existing Stage 3A sweep as the "
        "closest control points outside the candidate plateau.",
        "",
        "## B. Canonical-partition selection",
        "",
        "For each epsilon, all 30 partitions were compared by ARI and the medoid partition was "
        "defined as the run with the highest mean ARI to the other 29 runs.",
        "",
        "| epsilon | medoid run | medoid seed | mean ARI to others | min ARI to others | n classes |",
        "|---|---|---|---|---|---|",
    ]
    for run in medoid_runs:
        seed = BASE_SEED + epsilons.index(run["epsilon"]) * 10_000 + run["run"]
        lines.append(
            f"| {run['epsilon']:.6g} | {run['run']} | {seed} | "
            f"{run['medoid_mean_ari']:.4f} | {run['medoid_min_ari']:.4f} | "
            f"{len(run['reps'])} |"
        )

    ari_by_eps = {}
    for row in run_summary_rows:
        ari_by_eps.setdefault(row["epsilon"], []).append(row["mean_ari_to_others"])
    lines.extend(
        [
            "",
            "## C. Order sensitivity",
            "",
            "| epsilon | mean ARI (mean) | mean ARI (min) | n_classes std | n_classes min | n_classes max |",
            "|---|---|---|---|---|---|",
        ]
    )
    for eps in epsilons:
        runs = [r for r in run_summary_rows if r["epsilon"] == eps]
        means = [r["mean_ari_to_others"] for r in runs]
        nclasses = [r["n_classes"] for r in runs]
        lines.append(
            f"| {eps:.6g} | {np.mean(means):.4f} | {min(means):.4f} | "
            f"{np.std(nclasses, ddof=1):.3f} | {min(nclasses)} | {max(nclasses)} |"
        )

    lines.extend(
        [
            "",
            "## D. LAE class definitions",
            "",
            "Canonical class definitions are the medoid partitions. Class labels are neutral: "
            "LAE-1, LAE-2, ... and do not imply physical motif identity.",
            "",
            "| epsilon | class | total | stable iface | metastable iface | stable bulk | metastable bulk |",
            "|---|---|---|---|---|---|---|",
        ]
    )
    for row in canonical_rows:
        lines.append(
            f"| {row['epsilon']:.6g} | {row['class_id']} | {row['total']} | "
            f"{row['stable_interface']} | {row['metastable_interface']} | "
            f"{row['stable_bulk']} | {row['metastable_bulk']} |"
        )

    lines.extend(
        [
            "",
            "## E. Class compactness and non-mutual-equivalence",
            "",
            "| epsilon | class | n atoms | n pairs | diameter | median within | p95 within | rep radius | non-mutual-equivalence fraction |",
            "|---|---|---|---|---|---|---|---|---|",
        ]
    )
    for row in compactness_rows:
        lines.append(
            f"| {row['epsilon']:.6g} | {row['class_id']} | {row['n_atoms']} | {row['n_pairs']} | "
            f"{row['diameter']:.3f} | {row['median_within_distance']:.3f} | "
            f"{row['p95_within_distance']:.3f} | {row['representative_radius']:.3f} | "
            f"{row['non_mutual_equivalence_fraction']:.3f} |"
        )

    lines.extend(
        [
            "",
            "## F. Stable-interface LER",
            "",
            "Saved in `stage3b_ler_by_run.csv` and summarized in `stage3b_ler_summary.csv`.",
            "",
            "## G. Metastable-interface LER",
            "",
            "Same files as above; metastable columns are included in both outputs.",
            "",
            "## H. Stable-vs-metastable LER comparison",
            "",
            "Metrics use log2 Jensen-Shannon divergence. Means over 30 orderings:",
            "",
            "| epsilon | L1 mean | L1 std | JSD mean | JSD std | cosine mean | cosine std |",
            "|---|---|---|---|---|---|---|",
        ]
    )
    for row in comparison_summary:
        lines.append(
            f"| {row['epsilon']:.6g} | {row['l1']:.4f} | {row['l1_std']:.4f} | "
            f"{row['jsd']:.4f} | {row['jsd_std']:.4f} | "
            f"{row['cosine']:.4f} | {row['cosine_std']:.4f} |"
        )

    lines.extend(
        [
            "",
            "## I. PTM/CNA composition",
            "",
            "Counts and fractions by PTM, CNA, and combined PTM+CNA are in "
            "`stage3b_class_statistics.csv` for every canonical class.",
            "",
            "## J. PEATOM associations",
            "",
            "PEATOM median/IQR/mean/std per class are in `stage3b_class_statistics.csv`. "
            "The class-decomposition weighted interface mean is reported below.",
        ]
    )
    for eps in epsilons:
        run = next(r for r in medoid_runs if r["epsilon"] == eps)
        stable_mask = (structures == "stable") & (regions == "interface")
        meta_mask = (structures == "metastable") & (regions == "interface")
        stable_direct = float(np.mean(meta.iloc[stable_mask]["peatom"].astype(float)))
        meta_direct = float(np.mean(meta.iloc[meta_mask]["peatom"].astype(float)))
        stable_weighted = 0.0
        meta_weighted = 0.0
        for c in range(len(run["reps"])):
            idx = np.flatnonzero(run["assignment"] == c)
            class_stable = idx[(structures[idx] == "stable") & (regions[idx] == "interface")]
            class_meta = idx[(structures[idx] == "metastable") & (regions[idx] == "interface")]
            if len(class_stable):
                stable_weighted += len(class_stable) / int(stable_mask.sum()) * float(np.mean(meta.iloc[class_stable]["peatom"].astype(float)))
            if len(class_meta):
                meta_weighted += len(class_meta) / int(meta_mask.sum()) * float(np.mean(meta.iloc[class_meta]["peatom"].astype(float)))
        lines.append(
            f"- epsilon {eps:.6g}: stable direct {stable_direct:.6f} vs weighted {stable_weighted:.6f}; "
            f"metastable direct {meta_direct:.6f} vs weighted {meta_weighted:.6f}."
        )

    lines.extend(
        [
            "",
            "## K. Spatial localization",
            "",
            "At epsilon_mid, LAE class vs z is plotted separately for the stable and metastable "
            "structures in `stage3b_lae_vs_z_stable.png` and `stage3b_lae_vs_z_metastable.png`.",
            "",
            "## L. Epsilon sensitivity",
            "",
            "Class counts and largest-class fractions are plotted in "
            "`stage3b_class_population_vs_epsilon.png`; LER rank curves in "
            "`stage3b_ler_vs_epsilon.png`; compactness in `stage3b_class_compactness.png`.",
            "",
            "## M. Greedy/chaining limitations",
            "",
        ]
    )
    for eps in epsilons:
        n_components, largest_frac = s3a.connected_components(D2, float(eps))
        run = next(r for r in medoid_runs if r["epsilon"] == eps)
        lines.append(
            f"- epsilon {eps:.6g}: greedy classes {len(run['reps'])}, connected components "
            f"{n_components}, largest component fraction {largest_frac:.4f}."
        )
    lines.extend(
        [
            "",
            "Representative atoms, including both greedy representatives and class-medoid atoms, "
            "are saved in `stage3b_representative_atoms.csv`.",
            "",
            "## N. What the results support",
            "",
            "See the final recommendation below. LER is treated only as an environment-population "
            "fingerprint, not as a disorder or energy measure.",
            "",
            "## O. What the results do NOT support",
            "",
            "- No claim that the 5 classes are true physical motifs.",
            "- No claim that stable is more ordered or metastable is more disordered.",
            "- No direct GB-energy inference from PEATOM or LER.",
            "- No exact Priedeman reproduction.",
            "",
            "## P. Final recommendation",
            "",
            f"{outcome}: {decision}",
        ]
    )
    with open(OUTPUT_DIR / "stage3b_report.md", "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
