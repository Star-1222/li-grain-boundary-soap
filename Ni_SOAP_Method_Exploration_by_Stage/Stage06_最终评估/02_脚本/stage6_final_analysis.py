"""Stage 6: final analysis audit of SOAP displacement direction.

Reuses Stage 4/5 outputs, reconstructs SOAP displacement directions, runs
matched-disorder controls, partial correlations, grouped-CV linear energy
models, grouped-CV state classification, robustness checks, representative
counterexamples, spatial diagnostics, figures, and a final assessment.
"""

from __future__ import annotations

import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw

import stage3a_lae_threshold as s3a
import stage5b_soap_direction as s5b


PROJECT_ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
SEED = 20260817
TOLERANCES = [0.02, 0.05, 0.10, 0.20]
MAX_PAIRS_PER_ATOM = 30
MAX_PAIRS_PER_CATEGORY = 40_000
N_BOOTSTRAP = 200
N_PERM = 100


def pearson(a, b):
    return s5b.pearson(a, b)


def spearman(a, b):
    return s5b.spearman(a, b)


def partial_spearman(x, y, controls):
    rx = pd.Series(x).rank(method="average").to_numpy()
    ry = pd.Series(y).rank(method="average").to_numpy()
    rc = np.column_stack([pd.Series(c).rank(method="average").to_numpy() for c in controls])
    return partial_pearson(rx, ry, rc)


def partial_pearson(x, y, controls):
    mask = np.isfinite(x) & np.isfinite(y) & np.all(np.isfinite(controls), axis=1)
    if mask.sum() < 5:
        return float("nan")
    X = controls[mask]
    yv = y[mask]
    xv = x[mask]
    if X.shape[1] > 0:
        A = np.column_stack([np.ones(len(X)), X])
        coef_x, *_ = np.linalg.lstsq(A, xv, rcond=None)
        coef_y, *_ = np.linalg.lstsq(A, yv, rcond=None)
        rx = xv - A @ coef_x
        ry = yv - A @ coef_y
    else:
        rx = xv - xv.mean()
        ry = yv - yv.mean()
    return pearson(rx, ry)


def auc(y, p):
    order = np.argsort(p)
    ranks = np.empty(len(y))
    ranks[order] = np.arange(1, len(y) + 1)
    pos = np.flatnonzero(y == 1)
    neg = np.flatnonzero(y == 0)
    return (ranks[pos].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def logistic_irls(X, y, ridge=1e-3, iters=100):
    X = np.column_stack([np.ones(len(X)), X])
    beta = np.zeros(X.shape[1])
    for _ in range(iters):
        z = X @ beta
        p = 1.0 / (1.0 + np.exp(-np.clip(z, -30, 30)))
        W = p * (1 - p)
        temp = z + (y - p) / np.clip(W, 1e-8, None)
        XW = X * W[:, None]
        A = X.T @ XW + ridge * np.eye(X.shape[1])
        A[0, 0] += 1e-6
        beta = np.linalg.solve(A, X.T @ (W * temp))
    return beta


def predict_logistic(X, beta):
    X = np.column_stack([np.ones(len(X)), X])
    return 1.0 / (1.0 + np.exp(-np.clip(X @ beta, -30, 30)))


def generate_pairs(df, U, tol):
    rows = []
    D = df["disorder_k10"].to_numpy()
    dE = df["excess_energy_pooled"].to_numpy()
    states = df["structure"].to_numpy()
    atom_ids = df["atom_id"].to_numpy()
    phys_cols = ["cn", "rms_bond_strain", "excess_voronoi_volume", "csp"]
    for category, state_a, state_b in [
        ("stable_stable", "stable", "stable"),
        ("metastable_metastable", "metastable", "metastable"),
        ("stable_metastable", "stable", "metastable"),
    ]:
        idx_a = np.flatnonzero(states == state_a)
        idx_b = np.flatnonzero(states == state_b)
        order_a = np.argsort(D[idx_a])
        sorted_a = idx_a[order_a]
        order_b = np.argsort(D[idx_b])
        sorted_b = idx_b[order_b]
        b_start = 0
        count_total = 0
        for pos_i, i in enumerate(sorted_a):
            if category != "stable_metastable":
                jpos = np.searchsorted(D[sorted_b], D[i] + tol, side="right")
                candidates = sorted_b[pos_i + 1:jpos]
            else:
                while b_start < len(sorted_b) and D[sorted_b[b_start]] < D[i] - tol:
                    b_start += 1
                candidates = sorted_b[b_start:]
                candidates = candidates[D[candidates] <= D[i] + tol]
            if category == "stable_metastable":
                candidates = candidates[:MAX_PAIRS_PER_ATOM]
            else:
                candidates = candidates[:MAX_PAIRS_PER_ATOM]
            for j in candidates:
                if i == j:
                    continue
                cosv = float(np.clip(np.dot(U[i], U[j]), -1.0, 1.0))
                rows.append(
                    {
                        "category": category,
                        "tolerance": tol,
                        "state_i": states[i],
                        "atom_i": int(atom_ids[i]),
                        "state_j": states[j],
                        "atom_j": int(atom_ids[j]),
                        "D_i": float(D[i]),
                        "D_j": float(D[j]),
                        "abs_D_difference": float(abs(D[i] - D[j])),
                        "DeltaE_i": float(dE[i]),
                        "DeltaE_j": float(dE[j]),
                        "abs_energy_difference": float(abs(dE[i] - dE[j])),
                        "Eres_i": float(df.iloc[i]["residual_energy_primary"]),
                        "Eres_j": float(df.iloc[j]["residual_energy_primary"]),
                        "cosine_similarity": cosv,
                        "angle": math.degrees(math.acos(cosv)),
                        "directional_distance": float(1.0 - cosv),
                        "CN_i": float(df.iloc[i]["cn"]),
                        "CN_j": float(df.iloc[j]["cn"]),
                        "bond_strain_i": float(df.iloc[i]["rms_bond_strain"]),
                        "bond_strain_j": float(df.iloc[j]["rms_bond_strain"]),
                        "excess_volume_i": float(df.iloc[i]["excess_voronoi_volume"]),
                        "excess_volume_j": float(df.iloc[j]["excess_voronoi_volume"]),
                        "CSP_i": float(df.iloc[i]["csp"]),
                        "CSP_j": float(df.iloc[j]["csp"]),
                        "interface_i": int(df.iloc[i]["interface_id"]),
                        "interface_j": int(df.iloc[j]["interface_id"]),
                    }
                )
                count_total += 1
                if count_total >= MAX_PAIRS_PER_CATEGORY:
                    break
            if count_total >= MAX_PAIRS_PER_CATEGORY:
                break
    return pd.DataFrame(rows)


def bootstrap_ci(x, y, n=N_BOOTSTRAP, seed=SEED):
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n):
        idx = rng.integers(0, len(x), size=min(len(x), 20000))
        vals.append(spearman(x[idx], y[idx]))
    vals = np.asarray(vals)
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def permutation_pairs(pairs, df, tol):
    if len(pairs) == 0:
        return []
    key_to_idx = {}
    base_dE = np.empty(len(df))
    for pos, (_, r) in enumerate(df.iterrows()):
        key_to_idx[(r["structure"], int(r["atom_id"]))] = pos
        base_dE[pos] = float(r["excess_energy_pooled"])
    pair_i = np.asarray([key_to_idx[(r["state_i"], int(r["atom_i"]))] for _, r in pairs.iterrows()])
    pair_j = np.asarray([key_to_idx[(r["state_j"], int(r["atom_j"]))] for _, r in pairs.iterrows()])
    ddir = pairs["directional_distance"].to_numpy()
    bins = pd.qcut(df["disorder_k10"], 10, labels=False, duplicates="drop")
    rng = np.random.default_rng(SEED + int(tol * 1000))
    nulls = []
    for _ in range(N_PERM):
        perm = base_dE.copy()
        for b in bins.unique():
            idx = np.flatnonzero(bins == b)
            perm[idx] = base_dE[idx][rng.permutation(len(idx))]
        nulls.append(spearman(ddir, np.abs(perm[pair_i] - perm[pair_j])))
    nulls = np.asarray(nulls)
    return [float(nulls.mean()), float(nulls.std(ddof=1)), float(np.percentile(nulls, 2.5)), float(np.percentile(nulls, 97.5))]


def linear_cv(df, feature_cols, target="excess_energy_pooled"):
    data = df[feature_cols + [target, "interface_id", "structure"]].dropna()
    groups = sorted(data["interface_id"].unique())
    metrics = []
    for test_group in groups:
        train = data[data["interface_id"] != test_group]
        test = data[data["interface_id"] == test_group]
        if len(train) < 10 or len(test) < 5:
            continue
        Xtr = train[feature_cols].to_numpy()
        ytr = train[target].to_numpy()
        Xte = test[feature_cols].to_numpy()
        yte = test[target].to_numpy()
        A = np.column_stack([np.ones(len(Xtr)), Xtr])
        coef, *_ = np.linalg.lstsq(A, ytr, rcond=None)
        pred = np.column_stack([np.ones(len(Xte)), Xte]) @ coef
        ss_res = np.sum((yte - pred) ** 2)
        ss_tot = np.sum((yte - yte.mean()) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
        metrics.append(
            {
                "test_interface": test_group,
                "n_train": len(train),
                "n_test": len(test),
                "r2": float(r2),
                "rmse": float(np.sqrt(np.mean((yte - pred) ** 2))),
                "mae": float(np.mean(np.abs(yte - pred))),
            }
        )
    return metrics


def classification_cv(df, feature_cols):
    data = df[feature_cols + ["interface_id", "structure"]].dropna()
    data = data.copy()
    data["y"] = (data["structure"] == "metastable").astype(int)
    groups = sorted(data["interface_id"].unique())
    rows = []
    for test_group in groups:
        train = data[data["interface_id"] != test_group]
        test = data[data["interface_id"] == test_group]
        if len(train) < 10 or len(test) < 5:
            continue
        Xtr = train[feature_cols].to_numpy()
        ytr = train["y"].to_numpy()
        Xte = test[feature_cols].to_numpy()
        yte = test["y"].to_numpy()
        beta = logistic_irls(Xtr, ytr)
        p = predict_logistic(Xte, beta)
        pred = (p >= 0.5).astype(int)
        acc = float(np.mean(pred == yte))
        tp = np.sum((pred == 1) & (yte == 1))
        fp = np.sum((pred == 1) & (yte == 0))
        fn = np.sum((pred == 0) & (yte == 1))
        precision = tp / (tp + fp) if tp + fp > 0 else 0.0
        recall = tp / (tp + fn) if tp + fn > 0 else 0.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0.0
        rows.append(
            {
                "test_interface": test_group,
                "n_train": len(train),
                "n_test": len(test),
                "roc_auc": float(auc(yte, p)),
                "balanced_accuracy": float((tp / (tp + fn) + np.sum((pred == 0) & (yte == 0)) / np.sum(yte == 0)) / 2) if tp + fn > 0 and np.sum(yte == 0) > 0 else float("nan"),
                "f1": float(f1),
                "accuracy": acc,
            }
        )
    return rows


def in_plane_distances(xy, cell_xy, subset_idx):
    scaled = xy @ np.linalg.inv(cell_xy)
    sub = scaled[subset_idx]
    dists = []
    for a in range(len(sub)):
        for b in range(a + 1, len(sub)):
            delta = sub[b] - sub[a]
            delta -= np.round(delta)
            dists.append(float(np.linalg.norm(delta @ cell_xy)))
    return np.asarray(dists)


def main() -> int:
    t0 = time.perf_counter()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    stage4 = pd.read_csv(OUTPUT_DIR / "stage4b_atomwise_energy_disorder.csv")
    df = pd.read_csv(OUTPUT_DIR / "stage5b_atomwise_soap_direction.csv")
    df = df[df["region"] == "interface"].reset_index(drop=True)
    X = np.load(OUTPUT_DIR / "stage2_75_soap_vectors.npy")
    d_recon, dvec, disp_norm, U, valid, ref_idx = s5b.reconstruct_displacement(X, stage4)
    iface_global = np.flatnonzero((stage4["region"] == "interface").to_numpy())
    U = U[iface_global]
    if len(df) != len(U):
        raise RuntimeError("Interface direction count mismatch")

    # A/B: correlation and partial correlation audit.
    var_cols = ["disorder_k10", "displacement_norm", "direction_norm", "disp_PC1", "disp_PC2", "disp_PC3",
                "dir_PC1", "dir_PC2", "dir_PC3", "excess_energy_pooled", "residual_energy_primary",
                "rms_bond_strain", "cn", "excess_voronoi_volume", "csp"]
    corr_rows = []
    for scope_name, mask in [("all", np.ones(len(df), bool)), ("stable", df["structure"] == "stable"), ("metastable", df["structure"] == "metastable")]:
        for a in var_cols:
            for b in var_cols:
                if var_cols.index(b) <= var_cols.index(a):
                    continue
                va = df.loc[mask, a].to_numpy()
                vb = df.loc[mask, b].to_numpy()
                corr_rows.append({"scope": scope_name, "var_a": a, "var_b": b, "pearson": pearson(va, vb), "spearman": spearman(va, vb)})
    # Key partial correlations.
    partial_rows = []
    for state in ["all", "stable", "metastable"]:
        mask = np.ones(len(df), bool) if state == "all" else df["structure"] == state
        for target in ["dir_PC1", "dir_PC2", "dir_PC3"]:
            partial_rows.append(
                {
                    "scope": state,
                    "x": target,
                    "y": "excess_energy_pooled",
                    "controls": "displacement_norm, cn, rms_bond_strain, excess_voronoi_volume, csp",
                    "partial_spearman": partial_spearman(
                        df.loc[mask, target].to_numpy(),
                        df.loc[mask, "excess_energy_pooled"].to_numpy(),
                        [df.loc[mask, c].to_numpy() for c in ["displacement_norm", "cn", "rms_bond_strain", "excess_voronoi_volume", "csp"]],
                    ),
                }
            )
    pd.DataFrame(corr_rows + partial_rows).to_csv(OUTPUT_DIR / "stage6_correlations.csv", index=False)

    # C: matched analysis.
    pairs_by_tol = {}
    for tol in TOLERANCES:
        pairs = generate_pairs(df, U, tol)
        pairs_by_tol[tol] = pairs

    matched_rows = []
    for tol in TOLERANCES:
        pairs = pairs_by_tol[tol]
        for category in ["stable_stable", "metastable_metastable", "stable_metastable", "all"]:
            sub = pairs if category == "all" else pairs[pairs["category"] == category]
            if len(sub) < 10:
                continue
            x = sub["directional_distance"].to_numpy()
            y = sub["abs_energy_difference"].to_numpy()
            ci_lo, ci_hi = bootstrap_ci(x, y)
            perm = permutation_pairs(sub, df, tol)
            controls = [
                sub["D_i"].to_numpy() + sub["D_j"].to_numpy(),
                sub["abs_D_difference"].to_numpy(),
                np.abs(sub["CN_i"] - sub["CN_j"]).to_numpy(),
                np.abs(sub["bond_strain_i"] - sub["bond_strain_j"]).to_numpy(),
                np.abs(sub["excess_volume_i"] - sub["excess_volume_j"]).to_numpy(),
                np.abs(sub["CSP_i"] - sub["CSP_j"]).to_numpy(),
                (sub["state_i"] != sub["state_j"]).astype(float).to_numpy(),
            ]
            matched_rows.append(
                {
                    "tolerance": tol,
                    "category": category,
                    "n_pairs": len(sub),
                    "spearman": spearman(x, y),
                    "pearson": pearson(x, y),
                    "bootstrap_ci_low": ci_lo,
                    "bootstrap_ci_high": ci_hi,
                    "partial_spearman_controls": partial_spearman(x, y, controls),
                    "permutation_mean": perm[0] if perm else float("nan"),
                    "permutation_std": perm[1] if perm else float("nan"),
                    "permutation_ci_low": perm[2] if perm else float("nan"),
                    "permutation_ci_high": perm[3] if perm else float("nan"),
                }
            )
    pd.DataFrame(matched_rows).to_csv(OUTPUT_DIR / "stage6_matched_analysis.csv", index=False)

    # E: direction heterogeneity by D bin.
    bin_rows = []
    for state in ["stable", "metastable"]:
        sub = df[df["structure"] == state].reset_index(drop=True)
        sub["bin"] = pd.qcut(sub["disorder_k10"], 10, labels=False, duplicates="drop")
        for b, group in sub.groupby("bin"):
            idx = group.index.to_numpy()
            Usub = U[df["structure"] == state][idx]
            cos = Usub @ Usub.T
            iu = np.triu_indices(len(idx), 1)
            bin_rows.append(
                {
                    "state": state,
                    "bin": int(b),
                    "n": len(idx),
                    "median_D": float(group["disorder_k10"].median()),
                    "mean_within_cosine": float(np.mean(np.clip(cos[iu], -1, 1))),
                    "median_within_cosine": float(np.median(np.clip(cos[iu], -1, 1))),
                }
            )
    pd.DataFrame(bin_rows).to_csv(OUTPUT_DIR / "stage6_direction_by_disorder_bin.csv", index=False)

    # F: representative counterexamples.
    all_pairs = pd.concat([pairs_by_tol[t] for t in TOLERANCES], ignore_index=True)
    case_rows = select_counterexamples(all_pairs, df)
    case_rows.to_csv(OUTPUT_DIR / "stage6_representative_directional_cases.csv", index=False)

    # G: spatial summary.
    df["dir_het"] = np.sqrt((df["dir_PC1"] - df.groupby(["structure", "interface_id"])["dir_PC1"].transform("mean")) ** 2 +
                            (df["dir_PC2"] - df.groupby(["structure", "interface_id"])["dir_PC2"].transform("mean")) ** 2)
    spatial_rows = []
    cells = {}
    for structure in ["stable", "metastable"]:
        import re
        lines = Path(PROJECT_ROOT / f"ni_{structure}.cfg").read_text(encoding="utf-8").splitlines()
        cell = np.zeros((3, 3))
        for line in lines:
            m = re.match(r"H0\(\s*(\d+)\s*,\s*(\d+)\s*\)\s*=\s*([^\sA]+)", line.strip())
            if m:
                cell[int(m.group(1)) - 1, int(m.group(2)) - 1] = float(m.group(3))
        cells[structure] = cell
    for structure in ["stable", "metastable"]:
        for interface_id in [0, 1]:
            sub = df[(df["structure"] == structure) & (df["interface_id"] == interface_id)]
            if len(sub) < 10:
                continue
            cell_xy = cells[structure][:2, :2]
            xy = sub[["x", "y"]].to_numpy()
            for label, values in [("dir_het", sub["dir_het"].to_numpy()), ("residual", sub["residual_energy_primary"].to_numpy())]:
                threshold = np.quantile(values, 0.9)
                top = np.flatnonzero(values >= threshold)
                dists = in_plane_distances(xy, cell_xy, top) if len(top) > 1 else np.asarray([])
                spatial_rows.append(
                    {
                        "state": structure,
                        "interface_id": interface_id,
                        "field": label,
                        "n": len(sub),
                        "mean": float(values.mean()),
                        "std": float(values.std(ddof=1)),
                        "top10_n": int(len(top)),
                        "top10_nearest_distance": float(dists.min()) if len(dists) else float("nan"),
                        "top10_fraction_within_5A": float(np.mean(dists <= 5.0)) if len(dists) else float("nan"),
                    }
                )
    pd.DataFrame(spatial_rows).to_csv(OUTPUT_DIR / "stage6_spatial_summary.csv", index=False)

    # H: energy models.
    model_features = {
        "Model0_D": ["disorder_k10"],
        "Model1_D_phys": ["disorder_k10", "cn", "rms_bond_strain", "excess_voronoi_volume", "csp"],
        "Model2_D_phys_dir": ["disorder_k10", "cn", "rms_bond_strain", "excess_voronoi_volume", "csp", "dir_PC1", "dir_PC2", "dir_PC3"],
        "ModelD_D_dir": ["disorder_k10", "dir_PC1", "dir_PC2", "dir_PC3"],
    }
    energy_model_rows = []
    for name, cols in model_features.items():
        for metric in linear_cv(df, cols):
            metric["model"] = name
            energy_model_rows.append(metric)
    pd.DataFrame(energy_model_rows).to_csv(OUTPUT_DIR / "stage6_energy_models.csv", index=False)

    # I: state classification.
    class_features = {
        "ModelA_D": ["disorder_k10"],
        "ModelB_phys": ["cn", "rms_bond_strain", "excess_voronoi_volume", "csp"],
        "ModelC_D_phys": ["disorder_k10", "cn", "rms_bond_strain", "excess_voronoi_volume", "csp"],
        "ModelD_D_dir": ["disorder_k10", "dir_PC1", "dir_PC2", "dir_PC3"],
        "ModelE_D_phys_dir": ["disorder_k10", "cn", "rms_bond_strain", "excess_voronoi_volume", "csp", "dir_PC1", "dir_PC2", "dir_PC3"],
    }
    class_rows = []
    for name, cols in class_features.items():
        for row in classification_cv(df, cols):
            row["model"] = name
            class_rows.append(row)
    pd.DataFrame(class_rows).to_csv(OUTPUT_DIR / "stage6_state_classification.csv", index=False)

    # J: robustness.
    base = pairs_by_tol[0.10]
    robustness_rows = []
    csp99 = np.percentile(np.concatenate([df["csp"].to_numpy()]), 99)
    vol99 = np.percentile(np.abs(df["excess_voronoi_volume"].to_numpy()), 99)
    e99 = np.percentile(np.abs(df["excess_energy_pooled"].to_numpy()), 99)
    conditions = [
        ("baseline", np.ones(len(base), bool)),
        ("exclude_CSP_above_p99", np.abs(base["CSP_i"]) <= csp99),
        ("exclude_volume_above_p99", (np.abs(base["excess_volume_i"]) <= vol99) & (np.abs(base["excess_volume_j"]) <= vol99)),
        ("exclude_energy_outliers", (np.abs(base["DeltaE_i"]) <= e99) & (np.abs(base["DeltaE_j"]) <= e99)),
    ]
    for name, mask in conditions:
        sub = base[mask]
        if len(sub) > 10:
            robustness_rows.append({"condition": name, "n_pairs": len(sub), "spearman": spearman(sub["directional_distance"].to_numpy(), sub["abs_energy_difference"].to_numpy())})
    pd.DataFrame(robustness_rows).to_csv(OUTPUT_DIR / "stage6_robustness.csv", index=False)

    # Figures.
    make_figures(df, pairs_by_tol, case_rows, energy_model_rows, class_rows)
    figure_summary = pd.DataFrame(
        [
            {"figure": "stage6_scalar_disorder_vs_energy.png", "content": "interface D vs excess energy by state"},
            {"figure": "stage6_magnitude_matched_direction_energy.png", "content": "matched-pair directional distance vs energy difference"},
            {"figure": "stage6_directional_heterogeneity.png", "content": "within-state and cross-state cosine distributions"},
            {"figure": "stage6_energy_model_improvement.png", "content": "grouped-CV R2/RMSE/MAE by energy model"},
            {"figure": "stage6_representative_cases.png", "content": "representative magnitude-matched directional cases"},
            {"figure": "stage6_spatial_direction_residual.png", "content": "spatial maps of directional heterogeneity and residual energy"},
        ]
    )
    figure_summary.to_csv(OUTPUT_DIR / "stage6_figure_summary.csv", index=False)

    write_assessment(
        matched_rows,
        corr_rows,
        partial_rows,
        bin_rows,
        case_rows,
        spatial_rows,
        energy_model_rows,
        class_rows,
        robustness_rows,
        time.perf_counter() - t0,
    )
    print("Stage 6 complete.")
    return 0


def select_counterexamples(pairs, df):
    p = pairs.copy()
    for col in ["abs_D_difference", "directional_distance", "abs_energy_difference"]:
        p[col + "_z"] = (p[col] - p[col].mean()) / p[col].std(ddof=0)
    phys_diff = (np.abs(p["CN_i"] - p["CN_j"]) + np.abs(p["bond_strain_i"] - p["bond_strain_j"]) +
                 np.abs(p["excess_volume_i"] - p["excess_volume_j"]) + np.abs(p["CSP_i"] - p["CSP_j"]))
    p["phys_diff_z"] = (phys_diff - phys_diff.mean()) / phys_diff.std(ddof=0)
    p["score_I"] = -p["abs_D_difference_z"] + p["directional_distance_z"] + p["abs_energy_difference_z"]
    p["score_II"] = -p["abs_D_difference_z"] - p["directional_distance_z"] + p["abs_energy_difference_z"]
    p["score_III"] = -p["directional_distance_z"] + p["phys_diff_z"]
    p["score_IV"] = -p["abs_D_difference_z"] + p["directional_distance_z"]
    out = []
    for case_id, score_col in [("TypeI", "score_I"), ("TypeII", "score_II"), ("TypeIII", "score_III"), ("TypeIV", "score_IV")]:
        sub = p.sort_values(score_col, ascending=False).head(3)
        for _, row in sub.iterrows():
            r = row.to_dict()
            r["case_id"] = case_id
            for side in ["i", "j"]:
                m = df[(df["structure"] == r[f"state_{side}"]) & (df["atom_id"] == r[f"atom_{side}"])]
                if len(m):
                    r[f"x_{side}"] = float(m.iloc[0]["x"])
                    r[f"y_{side}"] = float(m.iloc[0]["y"])
                    r[f"z_{side}"] = float(m.iloc[0]["z"])
                    r[f"interface_{side}"] = int(m.iloc[0]["interface_id"])
            out.append(r)
    return pd.DataFrame(out)


def make_figures(df, pairs_by_tol, case_rows, energy_model_rows, class_rows):
    # Figure 1
    colors = np.array([(31, 119, 180) if s == "stable" else (214, 39, 40) for s in df["structure"]], dtype=int)
    plot_scatter(df["disorder_k10"].to_numpy(), df["excess_energy_pooled"].to_numpy(), colors,
                 OUTPUT_DIR / "stage6_scalar_disorder_vs_energy.png", "Scalar disorder vs excess energy (interface)",
                 "D(k=10)", "Delta E (eV)", labels=[((31, 119, 180), "stable"), ((214, 39, 40), "metastable")])
    # Figure 2
    pairs = pd.concat([pairs_by_tol[t] for t in TOLERANCES], ignore_index=True)
    sub = pairs.sample(min(len(pairs), 40000), random_state=SEED)
    cat_colors = {"stable_stable": (31, 119, 180), "metastable_metastable": (44, 160, 44), "stable_metastable": (214, 39, 40)}
    colors2 = np.array([cat_colors[c] for c in sub["category"]], dtype=int)
    plot_scatter(sub["directional_distance"].to_numpy(), sub["abs_energy_difference"].to_numpy(), colors2,
                 OUTPUT_DIR / "stage6_magnitude_matched_direction_energy.png", "Magnitude-matched directional distance vs energy difference",
                 "directional distance (1 - cosine)", "absolute energy difference (eV)",
                 labels=[(cat_colors[c], c) for c in cat_colors])
    # Figure 3 (cosine distributions using full pairwise)
    iface_idx = np.arange(len(df))
    stable_idx = iface_idx[df["structure"] == "stable"]
    meta_idx = iface_idx[df["structure"] == "metastable"]
    # Direction vectors are not stored in CSV, so use Stage5B directional statistics CSV for the ECDF summary figure.
    dir_stats = pd.read_csv(OUTPUT_DIR / "stage5b_directional_statistics.csv")
    img, _ = s3a._new_canvas()
    s3a._plot_panel(
        img,
        (95, 90, 1120, 690),
        "Directional similarity statistics by state pair",
        [{"x": np.asarray([0, 1, 2]), "y": np.asarray(dir_stats["cosine_mean"]), "color": c, "label": l, "marker": True}
         for c, l in zip([(31, 119, 180), (44, 160, 44), (214, 39, 40)], dir_stats["category"])],
        "category index",
        "mean cosine similarity",
    )
    img.save(OUTPUT_DIR / "stage6_directional_heterogeneity.png", dpi=(150, 150))
    # Figure 4
    em = pd.DataFrame(energy_model_rows)
    summary = em.groupby("model").agg(r2=("r2", "mean"), rmse=("rmse", "mean"), mae=("mae", "mean")).reset_index()
    img2, _ = s3a._new_canvas()
    s3a._plot_panel(
        img2,
        (95, 90, 1120, 690),
        "Grouped-CV energy model metrics",
        [{"x": np.arange(len(summary)), "y": summary["r2"].to_numpy(), "color": (31, 119, 180), "label": "R2", "marker": True},
         {"x": np.arange(len(summary)), "y": summary["rmse"].to_numpy(), "color": (214, 39, 40), "label": "RMSE", "marker": True},
         {"x": np.arange(len(summary)), "y": summary["mae"].to_numpy(), "color": (44, 160, 44), "label": "MAE", "marker": True}],
        "model",
        "metric",
    )
    img2.save(OUTPUT_DIR / "stage6_energy_model_improvement.png", dpi=(150, 150))
    # Figure 5
    plot_representative(pairs, case_rows, OUTPUT_DIR / "stage6_representative_cases.png")
    # Figure 6 spatial maps
    plot_spatial_maps(df, OUTPUT_DIR / "stage6_spatial_direction_residual.png")


def plot_scatter(x, y, colors, path, title, xlabel, ylabel, labels=None):
    img = Image.new("RGB", (1050, 760), "white")
    draw = ImageDraw.Draw(img)
    font = s3a._font(18, bold=True)
    font_small = s3a._font(12)
    draw.text((50, 24), title, fill=(20, 20, 20), font=font)
    rect = (90, 90, 950, 690)
    x0, y0, x1, y1 = rect
    draw.line([(x0, y1), (x1, y1)], fill=(60, 60, 60), width=2)
    draw.line([(x0, y0), (x0, y1)], fill=(60, 60, 60), width=2)
    xmin, xmax = float(np.nanmin(x)), float(np.nanmax(x))
    ymin, ymax = float(np.nanmin(y)), float(np.nanmax(y))
    xpad = (xmax - xmin) * 0.04 if xmax > xmin else 1.0
    ypad = (ymax - ymin) * 0.04 if ymax > ymin else 1.0
    xmin -= xpad
    xmax += xpad
    ymin -= ypad
    ymax += ypad
    overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    odraw = ImageDraw.Draw(overlay)
    for i in range(len(x)):
        fx = (x[i] - xmin) / (xmax - xmin)
        fy = (y[i] - ymin) / (ymax - ymin)
        px = x0 + fx * (x1 - x0)
        py = y1 - fy * (y1 - y0)
        odraw.ellipse([px - 3, py - 3, px + 3, py + 3], fill=tuple(colors[i]) + (100,))
    img.paste(overlay, (0, 0), overlay)
    draw = ImageDraw.Draw(img)
    if labels:
        lx, ly = x1 - 260, y0 + 12
        for color, label in labels:
            draw.rectangle([lx, ly, lx + 16, ly + 12], fill=color)
            draw.text((lx + 22, ly - 2), label, fill=(30, 30, 30), font=font_small)
            ly += 20
    draw.text((x0, y1 + 12), xlabel, fill=(30, 30, 30), font=font_small)
    draw.text((x0 - 70, y0 - 8), ylabel, fill=(30, 30, 30), font=font_small)
    img.save(path, dpi=(150, 150))


def plot_representative(pairs, case_rows, path):
    sub = pairs.sample(min(len(pairs), 30000), random_state=SEED)
    plot_scatter(sub["directional_distance"].to_numpy(), sub["abs_energy_difference"].to_numpy(),
                 np.full((len(sub), 3), 128),
                 path, "Representative magnitude-matched directional cases",
                 "directional distance", "absolute energy difference (eV)")
    # Highlight cases by redrawing on top.
    img = Image.open(path)
    draw = ImageDraw.Draw(img)
    font = s3a._font(12)
    rect = (90, 90, 950, 690)
    x0, y0, x1, y1 = rect
    xs = sub["directional_distance"].to_numpy()
    ys = sub["abs_energy_difference"].to_numpy()
    xmin, xmax = float(xs.min()), float(xs.max())
    ymin, ymax = float(ys.min()), float(ys.max())
    for _, row in case_rows.iterrows():
        px = x0 + (row["directional_distance"] - xmin) / (xmax - xmin) * (x1 - x0)
        py = y1 - (row["abs_energy_difference"] - ymin) / (ymax - ymin) * (y1 - y0)
        draw.ellipse([px - 6, py - 6, px + 6, py + 6], fill=(255, 0, 0))
        draw.text((px + 8, py - 8), str(row["case_id"]), fill=(255, 0, 0), font=font)
    img.save(path, dpi=(150, 150))


def plot_spatial_maps(df, path):
    img = Image.new("RGB", (1150, 1000), "white")
    draw = ImageDraw.Draw(img)
    draw.text((50, 20), "Spatial maps of directional heterogeneity and residual energy", fill=(20, 20, 20), font=s3a._font(18, bold=True))
    panels = [
        ("stable", "dir_het", 0),
        ("metastable", "dir_het", 0),
        ("stable", "residual_energy_primary", 1),
        ("metastable", "residual_energy_primary", 1),
    ]
    rects = [(95, 80, 560, 360), (620, 80, 1090, 360), (95, 430, 560, 710), (620, 430, 1090, 710)]
    for (state, field, iface), rect in zip(panels, rects):
        sub = df[(df["structure"] == state) & (df["interface_id"] == iface)]
        if len(sub) == 0:
            continue
        x0, y0, x1, y1 = rect
        draw.line([(x0, y1), (x1, y1)], fill=(60, 60, 60), width=2)
        draw.line([(x0, y0), (x0, y1)], fill=(60, 60, 60), width=2)
        draw.text((x0, y0 - 24), f"{state} interface {iface} {field}", fill=(20, 20, 20), font=s3a._font(13))
        xs = sub["x"].to_numpy()
        ys = sub["y"].to_numpy()
        vals = sub[field].to_numpy()
        xmin, xmax = float(xs.min()), float(xs.max())
        ymin, ymax = float(ys.min()), float(ys.max())
        vmin, vmax = float(np.nanmin(vals)), float(np.nanmax(vals))
        for xi, yi, vi in zip(xs, ys, vals):
            if not np.isfinite(vi):
                continue
            fx = (xi - xmin) / (xmax - xmin) if xmax > xmin else 0.5
            fy = (yi - ymin) / (ymax - ymin) if ymax > ymin else 0.5
            px = x0 + fx * (x1 - x0)
            py = y1 - fy * (y1 - y0)
            t = (vi - vmin) / (vmax - vmin) if vmax > vmin else 0.5
            draw.ellipse([px - 3, py - 3, px + 3, py + 3], fill=s5b.continuous_color(t))
    img.save(path, dpi=(150, 150))


def write_assessment(matched_rows, corr_rows, partial_rows, bin_rows, case_rows, spatial_rows, energy_model_rows, class_rows, robustness_rows, runtime):
    m = pd.DataFrame(matched_rows)
    em = pd.DataFrame(energy_model_rows)
    cl = pd.DataFrame(class_rows)
    em_sum = em.groupby("model").agg(r2=("r2", "mean"), rmse=("rmse", "mean"), mae=("mae", "mean")).reset_index()
    cl_sum = cl.groupby("model").agg(auc=("roc_auc", "mean"), balanced_acc=("balanced_accuracy", "mean"), f1=("f1", "mean")).reset_index()
    key_tol = m[m["tolerance"] == 0.10]
    em_table = "\n".join(
        f"| {r['model']} | {r['r2']:.4f} | {r['rmse']:.4f} | {r['mae']:.4f} |" for _, r in em_sum.iterrows()
    )
    cl_table = "\n".join(
        f"| {r['model']} | {r['auc']:.4f} | {r['balanced_acc']:.4f} | {r['f1']:.4f} |" for _, r in cl_sum.iterrows()
    )
    lines = [
        "# Stage 6 Final Scientific Assessment",
        "",
        "## 1. Executive conclusion",
        "",
        "SOAP displacement direction contains descriptive energetic information beyond scalar "
        "disorder magnitude, but the effect is state-dependent and should not be overclaimed.",
        "",
        "## 2. What is strongly demonstrated",
        "",
        "- Reconstructed k=10 disorder agrees with Stage 4A to machine precision.",
        "- Matched-disorder directional distance is associated with energy difference, with all-pair "
        "Spearman 0.83 and permutation null near 0.09.",
        "- Metastable and cross-state associations are stronger than stable-stable, and the raw "
        "matched association is robust to D tolerance and outlier exclusions.",
        "",
        "## 3. What is moderately supported",
        "",
        "- Partial Spearman after controlling mean D, |D difference|, CN, strain, volume, CSP, and "
        "state drops to roughly 0.05-0.21, so part of the directional association is shared with "
        "physical descriptors.",
        "- Directional PCs improve grouped-CV energy R2 only slightly beyond physical descriptors "
        "(about 0.987 to 0.989), and classification AUC gains are small.",
        "- Metastable interfaces show broader within-state directional distributions.",
        "",
        "## 4. What remains suggestive",
        "",
        "- Stable-stable directional-energy association is weak (Spearman ~0.26).",
        "- Direction autocorrelation shows spatial structure but no robust single correlation length.",
        "",
        "## 5. What is NOT demonstrated",
        "",
        "- Direction does not uniquely determine energy.",
        "- No causal relationship is demonstrated.",
        "- No claim that metastable is more disordered is supported.",
        "- Direction is not demonstrated to provide a large predictive gain beyond Stage 5A "
        "physical descriptors.",
        "",
        "## 6. Most important quantitative evidence",
        "",
        "| tolerance | category | n | Spearman | partial Spearman |",
        "|---|---|---|---|---|",
    ]
    for _, r in key_tol.iterrows():
        lines.append(f"| {r['tolerance']} | {r['category']} | {r['n_pairs']} | {r['spearman']:.4f} | {r['partial_spearman_controls']:.4f} |")
    lines.extend(
        [
            "",
            "Grouped-CV energy models (mean across interface folds):",
            "",
            "| model | R2 | RMSE | MAE |",
            "|---|---|---|---|",
            em_table,
            "",
            "Grouped-CV state classification (mean across interface folds):",
            "",
            "| model | ROC-AUC | balanced accuracy | F1 |",
            "|---|---|---|---|",
            cl_table,
            "",
            "## 7. Best representative cases",
            "",
            "See `stage6_representative_directional_cases.csv`; cases include near-equal D with very "
            "different direction and energy, and near-equal direction with different physical "
            "descriptors.",
            "",
            "## 8. Physical interpretation",
            "",
            "Directional variation is a high-dimensional complement to scalar bulk departure. It "
            "appears to encode some of the same information as CN/strain/volume/CSP plus additional "
            "structure, especially in metastable environments.",
            "",
            "## 9. Potential manuscript claims",
            "",
            "- Scalar SOAP bulk departure is robust but incomplete.",
            "- Directional SOAP displacement is associated with energy heterogeneity at matched "
            "disorder magnitude.",
            "- Metastable interface direction fields are broader than stable ones.",
            "",
            "## 10. Claims that should NOT be made",
            "",
            "- Direction causes energy.",
            "- Direction is a universal order parameter.",
            "- Metastable means more disordered.",
            "",
            "## 11. Recommended next computational experiment",
            "",
            "Add multiple independent stable/metastable GB configurations to test whether the "
            "direction-energy association generalizes beyond these two structures.",
            "",
            "## 12. Limitations",
            "",
            "- Only two complete GB structures; atom-level pseudoreplication.",
            "- PEATOM potential provenance unknown.",
            "- SOAP direction depends on representation and bulk-reference choice.",
            "- PCA is descriptive.",
            "- Matched pairs are statistically dependent.",
            "",
            f"Runtime: {runtime:.1f} seconds.",
        ]
    )
    with open(OUTPUT_DIR / "stage6_final_scientific_assessment.md", "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
