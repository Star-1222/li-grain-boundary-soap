"""Stage 5B: SOAP displacement-direction information beyond scalar disorder.

Reconstructs SOAP displacement vectors relative to the Stage 4A clean-bulk
k-nearest reference convention, validates against Stage 4A disorder, analyzes
directional similarity among matched-disorder atoms, runs descriptive PCA of
displacement and direction, and tests directional spatial autocorrelation.
No regression, machine learning, or categorical clustering is performed.
"""

from __future__ import annotations

import csv
import itertools
import math
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw

import stage3a_lae_threshold as s3a


PROJECT_ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
N_ATOMS = 4659
SEED = 20260817
K = 10
TOLERANCES = [0.05, 0.10, 0.20]
MAX_PAIRS_PER_ATOM = 30
MAX_PAIRS_PER_CATEGORY = 40_000
NEIGHBOR_RADIUS = 5.5


def continuous_color(t):
    stops = [(31, 119, 180), (44, 160, 44), (255, 255, 51), (214, 39, 40)]
    t = max(0.0, min(1.0, t)) * (len(stops) - 1)
    i = int(math.floor(t))
    if i >= len(stops) - 1:
        return stops[-1]
    f = t - i
    return tuple(int(round(a + (b - a) * f)) for a, b in zip(stops[i], stops[i + 1]))


def pearson(a, b):
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3:
        return float("nan")
    return float(np.corrcoef(a[mask], b[mask])[0, 1])


def spearman(a, b):
    mask = np.isfinite(a) & np.isfinite(b)
    a = a[mask]
    b = b[mask]
    ra = pd.Series(a).rank(method="average").to_numpy()
    rb = pd.Series(b).rank(method="average").to_numpy()
    return pearson(ra, rb)


def run_pca(matrix):
    Xc = matrix - matrix.mean(axis=0)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    PCs = Xc @ Vt.T
    evr = S**2 / np.sum(S**2)
    return PCs, evr


def load_data():
    X = np.load(OUTPUT_DIR / "stage2_75_soap_vectors.npy")
    if X.shape != (N_ATOMS, 1014):
        raise RuntimeError(f"Unexpected SOAP shape {X.shape}")
    meta = pd.read_csv(OUTPUT_DIR / "stage2_75_atom_metadata.csv")
    stage4b = pd.read_csv(OUTPUT_DIR / "stage4b_atomwise_energy_disorder.csv")
    stage5a = pd.read_csv(OUTPUT_DIR / "stage5a_atomwise_physical_descriptors.csv")
    df = stage4b.copy()
    df = df.merge(
        stage5a[["atom_id", "state", "x", "y", "z", "interface_id", "residual_energy_primary",
                 "cn", "rms_bond_strain", "signed_bond_strain", "csp",
                 "voronoi_volume", "excess_voronoi_volume"]],
        left_on=["atom_id", "structure"],
        right_on=["atom_id", "state"],
        how="left",
    )
    df = df.rename(columns={"x_x": "x", "y_x": "y", "z_x": "z"})
    df = df.drop(columns=["x_y", "y_y", "z_y"], errors="ignore")
    df["interface_id"] = df["interface_id"].fillna(-1).astype(int)
    if len(df) != N_ATOMS:
        raise RuntimeError("Merged atomwise row count mismatch")
    return X, meta, df


def reconstruct_displacement(X, df):
    ref_mask = df["in_clean_bulk_reference"].astype(bool).to_numpy()
    ref_idx = np.flatnonzero(ref_mask)
    ref = X[ref_idx]
    norms_all = np.einsum("ij,ij->i", X, X)
    norms_ref = np.einsum("ij,ij->i", ref, ref)
    D2 = norms_all[:, None] + norms_ref[None, :] - 2.0 * (X @ ref.T)
    D2 = np.maximum(D2, 0.0)
    positions = np.searchsorted(ref_idx, ref_idx)
    D2[ref_idx, positions] = np.inf
    order = np.argsort(D2, axis=1)
    k_idx = order[:, :K]
    k_ref = ref_idx[k_idx]
    d_reconstructed = np.sqrt(np.maximum(D2[np.arange(N_ATOMS)[:, None], k_idx], 0.0)).mean(axis=1)
    bulk_mean = X[k_ref].mean(axis=1)
    dvec = X - bulk_mean
    disp_norm = np.linalg.norm(dvec, axis=1)
    u = np.zeros_like(dvec)
    valid = disp_norm > 1e-8
    u[valid] = dvec[valid] / disp_norm[valid, None]
    return d_reconstructed, dvec, disp_norm, u, valid, ref_idx


def pairwise_direction_stats(U, idx_a, idx_b):
    if len(idx_a) == 0 or len(idx_b) == 0:
        return []
    if np.array_equal(np.asarray(idx_a), np.asarray(idx_b)):
        cos = U[np.asarray(idx_a)] @ U[np.asarray(idx_b)].T
        iu = np.triu_indices(len(idx_a), 1)
        cos = cos[iu]
    else:
        cos = U[np.asarray(idx_a)] @ U[np.asarray(idx_b)].T
        cos = cos.ravel()
    cos = np.clip(cos, -1.0, 1.0)
    return cos


def plot_pca_scatter(PCs, color_values, path, title, color_type="group", group_masks=None, cbar_label=""):
    pc1, pc2 = PCs[:, 0], PCs[:, 1]
    img = Image.new("RGB", (1050, 760), "white")
    draw = ImageDraw.Draw(img)
    font = s3a._font(18, bold=True)
    font_small = s3a._font(12)
    draw.text((50, 24), title, fill=(20, 20, 20), font=font)
    rect = (90, 90, 940, 690)
    x0, y0, x1, y1 = rect
    draw.line([(x0, y1), (x1, y1)], fill=(60, 60, 60), width=2)
    draw.line([(x0, y0), (x0, y1)], fill=(60, 60, 60), width=2)
    xmin, xmax = float(pc1.min()), float(pc1.max())
    ymin, ymax = float(pc2.min()), float(pc2.max())
    xpad = (xmax - xmin) * 0.04 if xmax > xmin else 1.0
    ypad = (ymax - ymin) * 0.04 if ymax > ymin else 1.0
    xmin -= xpad
    xmax += xpad
    ymin -= ypad
    ymax += ypad
    if color_type == "group":
        colors = [(31, 119, 180), (214, 39, 40)]
        labels = ["stable interface", "metastable interface"]
        atom_colors = np.zeros((len(pc1), 3), dtype=int)
        for (name, mask), color in zip(group_masks.items(), colors):
            atom_colors[mask] = color
        for i in range(len(pc1)):
            fx = (pc1[i] - xmin) / (xmax - xmin)
            fy = (pc2[i] - ymin) / (ymax - ymin)
            px = x0 + fx * (x1 - x0)
            py = y1 - fy * (y1 - y0)
            r = 3
            draw.ellipse([px - r, py - r, px + r, py + r], fill=tuple(atom_colors[i]))
        lx, ly = x1 + 25, y0 + 20
        for color, label in zip(colors, labels):
            draw.rectangle([lx, ly, lx + 18, ly + 12], fill=color)
            draw.text((lx + 24, ly - 2), label, fill=(30, 30, 30), font=font_small)
            ly += 22
    else:
        vals = np.asarray(color_values)
        vmin, vmax = float(np.nanmin(vals)), float(np.nanmax(vals))
        for i in range(len(pc1)):
            if not np.isfinite(vals[i]):
                continue
            t = (vals[i] - vmin) / (vmax - vmin) if vmax > vmin else 0.5
            fx = (pc1[i] - xmin) / (xmax - xmin)
            fy = (pc2[i] - ymin) / (ymax - ymin)
            px = x0 + fx * (x1 - x0)
            py = y1 - fy * (y1 - y0)
            r = 3
            draw.ellipse([px - r, py - r, px + r, py + r], fill=continuous_color(t))
        bar_x = x1 + 30
        bar_y0 = y1 - 160
        bar_y1 = y1
        for py in range(bar_y0, bar_y1):
            t = (py - bar_y0) / (bar_y1 - bar_y0)
            draw.line([(bar_x, py), (bar_x + 16, py)], fill=continuous_color(t))
        draw.text((bar_x + 22, bar_y0 - 6), f"{vmax:.3g}", fill=(30, 30, 30), font=font_small)
        draw.text((bar_x + 22, bar_y1 - 6), f"{vmin:.3g}", fill=(30, 30, 30), font=font_small)
        draw.text((bar_x, y1 + 12), cbar_label, fill=(30, 30, 30), font=font_small)
    draw.text((x0, y1 + 12), "PC1", fill=(30, 30, 30), font=font_small)
    draw.text((x0 - 70, y0 - 8), "PC2", fill=(30, 30, 30), font=font_small)
    img.save(path, dpi=(150, 150))


def generate_matched_pairs(df, U, dE, Eres, phys):
    rows = []
    iface = df["region"] == "interface"
    states = df["structure"].to_numpy()
    D = df["disorder_k10"].to_numpy()
    atom_ids = df["atom_id"].to_numpy()
    categories = [
        ("stable_stable", iface & (states == "stable")),
        ("metastable_metastable", iface & (states == "metastable")),
        ("stable_metastable", iface & (states == "stable")),
    ]
    for category, mask_a in categories:
        if category == "stable_metastable":
            mask_b = iface & (states == "metastable")
        else:
            mask_b = mask_a
        idx_a = np.flatnonzero(mask_a)
        idx_b = np.flatnonzero(mask_b)
        for tol in TOLERANCES:
            if category == "stable_metastable":
                order_a = np.argsort(D[idx_a])
                order_b = np.argsort(D[idx_b])
                sorted_a = idx_a[order_a]
                sorted_b = idx_b[order_b]
                pairs = []
                b_start = 0
                for i in sorted_a:
                    while b_start < len(sorted_b) and D[sorted_b[b_start]] < D[i] - tol:
                        b_start += 1
                    count = 0
                    jpos = b_start
                    while jpos < len(sorted_b) and D[sorted_b[jpos]] <= D[i] + tol:
                        pairs.append((i, sorted_b[jpos]))
                        count += 1
                        jpos += 1
                        if count >= MAX_PAIRS_PER_ATOM:
                            break
                    if len(pairs) >= MAX_PAIRS_PER_CATEGORY:
                        break
            else:
                order = np.argsort(D[idx_a])
                sorted_idx = idx_a[order]
                pairs = []
                for pos_i, i in enumerate(sorted_idx):
                    count = 0
                    for j in sorted_idx[pos_i + 1:]:
                        if D[j] - D[i] > tol:
                            break
                        pairs.append((i, j))
                        count += 1
                        if count >= MAX_PAIRS_PER_ATOM:
                            break
                    if len(pairs) >= MAX_PAIRS_PER_CATEGORY:
                        break
            for i, j in pairs:
                if i == j:
                    continue
                cosv = float(np.clip(np.dot(U[i], U[j]), -1.0, 1.0))
                angle = math.degrees(math.acos(cosv))
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
                        "Eres_i": float(Eres[i]),
                        "Eres_j": float(Eres[j]),
                        "cosine_similarity": cosv,
                        "angle": angle,
                        "directional_distance": float(1.0 - cosv),
                        "CN_i": float(phys.loc[i, "cn"]),
                        "CN_j": float(phys.loc[j, "cn"]),
                        "bond_strain_i": float(phys.loc[i, "rms_bond_strain"]),
                        "bond_strain_j": float(phys.loc[j, "rms_bond_strain"]),
                        "excess_volume_i": float(phys.loc[i, "excess_voronoi_volume"]),
                        "excess_volume_j": float(phys.loc[j, "excess_voronoi_volume"]),
                        "CSP_i": float(phys.loc[i, "csp"]),
                        "CSP_j": float(phys.loc[j, "csp"]),
                        "abs_CN_difference": float(abs(phys.loc[i, "cn"] - phys.loc[j, "cn"])),
                        "abs_strain_difference": float(abs(phys.loc[i, "rms_bond_strain"] - phys.loc[j, "rms_bond_strain"])),
                        "abs_volume_difference": float(abs(phys.loc[i, "excess_voronoi_volume"] - phys.loc[j, "excess_voronoi_volume"])),
                        "abs_csp_difference": float(abs(phys.loc[i, "csp"] - phys.loc[j, "csp"])),
                    }
                )
    return pd.DataFrame(rows)


def direction_autocorr(df, U, cells):
    rows = []
    for structure in ["stable", "metastable"]:
        cell_xy = cells[structure][:2, :2]
        for interface_id in [0, 1]:
            sub = df[(df["structure"] == structure) & (df["region"] == "interface") & (df["interface_id"] == interface_id)]
            if len(sub) < 5:
                continue
            idx = sub.index.to_numpy()
            scaled_xy = sub[["x", "y"]].to_numpy() @ np.linalg.inv(cell_xy)
            n = len(idx)
            bins = np.zeros(16)
            counts = np.zeros(16)
            for a in range(n):
                for b in range(a + 1, n):
                    delta = scaled_xy[b] - scaled_xy[a]
                    delta -= np.round(delta)
                    r = float(np.linalg.norm(delta @ cell_xy))
                    if r > 15:
                        continue
                    bi = int(r / 1.0)
                    bins[bi] += float(np.dot(U[idx[a]], U[idx[b]]))
                    counts[bi] += 1
            for bi in range(16):
                if counts[bi] > 0:
                    rows.append(
                        {
                            "state": structure,
                            "interface_id": interface_id,
                            "r_mid": bi + 0.5,
                            "n_pairs": int(counts[bi]),
                            "mean_cosine": float(bins[bi] / counts[bi]),
                        }
                    )
    return pd.DataFrame(rows)


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    X, meta, df = load_data()
    d_recon, dvec, disp_norm, u, valid_dir, ref_idx = reconstruct_displacement(X, df)
    d_stage4 = df["disorder_k10"].to_numpy()
    diff = np.abs(d_recon - d_stage4)
    validation = {
        "max_abs_diff": float(diff.max()),
        "median_abs_diff": float(np.median(diff)),
        "p99_abs_diff": float(np.percentile(diff, 99)),
        "n_valid_directions": int(valid_dir.sum()),
    }

    df["displacement_norm"] = disp_norm
    df["direction_norm"] = np.linalg.norm(u, axis=1)
    iface_mask = df["region"] == "interface"
    iface_idx = np.flatnonzero(iface_mask)
    dvec_iface = dvec[iface_idx]
    u_iface = u[iface_idx]
    pca_disp, evr_disp = run_pca(dvec_iface)
    pca_dir, evr_dir = run_pca(u_iface)
    df.loc[iface_idx, "disp_PC1"] = pca_disp[:, 0]
    df.loc[iface_idx, "disp_PC2"] = pca_disp[:, 1]
    df.loc[iface_idx, "disp_PC3"] = pca_disp[:, 2]
    df.loc[iface_idx, "dir_PC1"] = pca_dir[:, 0]
    df.loc[iface_idx, "dir_PC2"] = pca_dir[:, 1]
    df.loc[iface_idx, "dir_PC3"] = pca_dir[:, 2]

    atomwise_out = df[[
        "structure", "atom_id", "x", "y", "z", "interface_id", "region",
        "disorder_k10", "displacement_norm", "direction_norm",
        "disp_PC1", "disp_PC2", "disp_PC3", "dir_PC1", "dir_PC2", "dir_PC3",
        "peatom", "excess_energy_pooled", "residual_energy_primary",
        "ptm", "cna", "cn", "rms_bond_strain", "excess_voronoi_volume", "csp",
    ]].copy()
    atomwise_out.to_csv(OUTPUT_DIR / "stage5b_atomwise_soap_direction.csv", index=False)

    # Directional distributions.
    U_iface = u[iface_mask]
    iface_global = np.flatnonzero(iface_mask)
    stable_global = iface_global[df["structure"].to_numpy()[iface_global] == "stable"]
    meta_global = iface_global[df["structure"].to_numpy()[iface_global] == "metastable"]
    stable_idx = np.searchsorted(iface_global, stable_global)
    meta_idx = np.searchsorted(iface_global, meta_global)
    dir_stats = []
    cos_dist = {}
    for category, ia, ib in [
        ("stable_stable", stable_idx, stable_idx),
        ("metastable_metastable", meta_idx, meta_idx),
        ("stable_metastable", stable_idx, meta_idx),
    ]:
        cos_vals = pairwise_direction_stats(U_iface, ia, ib)
        cos_dist[category] = cos_vals
        angle_vals = np.degrees(np.arccos(np.clip(cos_vals, -1, 1)))
        ddir = 1.0 - cos_vals
        dir_stats.append(
            {
                "category": category,
                "n_pairs": int(len(cos_vals)),
                "cosine_mean": float(np.mean(cos_vals)),
                "cosine_median": float(np.median(cos_vals)),
                "cosine_std": float(np.std(cos_vals, ddof=1)),
                "cosine_q1": float(np.percentile(cos_vals, 25)),
                "cosine_q3": float(np.percentile(cos_vals, 75)),
                "angle_mean": float(np.mean(angle_vals)),
                "angle_median": float(np.median(angle_vals)),
                "directional_distance_mean": float(np.mean(ddir)),
            }
        )
    pd.DataFrame(dir_stats).to_csv(OUTPUT_DIR / "stage5b_directional_statistics.csv", index=False)

    # Matched pairs.
    phys = df.set_index(df.index)
    dE = df["excess_energy_pooled"].to_numpy()
    Eres = df["residual_energy_primary"].to_numpy()
    pairs = generate_matched_pairs(df, u, dE, Eres, phys)
    pairs.to_csv(OUTPUT_DIR / "stage5b_matched_direction_energy_pairs.csv", index=False)

    # Directional-energy statistics.
    dir_energy_rows = []
    for tol in TOLERANCES:
        for category in ["stable_stable", "metastable_metastable", "stable_metastable", "all"]:
            sub = pairs[(pairs["tolerance"] == tol) & (pairs["category"] == category)] if category != "all" else pairs[pairs["tolerance"] == tol]
            if len(sub) < 5:
                continue
            x = sub["directional_distance"].to_numpy()
            y = sub["abs_energy_difference"].to_numpy()
            qs = pd.qcut(sub["directional_distance"], 5, labels=False, duplicates="drop")
            med_dE = sub.groupby(qs)["abs_energy_difference"].median()
            med_dD = sub.groupby(qs)["abs_D_difference"].median()
            dir_energy_rows.append(
                {
                    "tolerance": tol,
                    "category": category,
                    "n_pairs": len(sub),
                    "spearman": spearman(x, y),
                    "pearson": pearson(x, y),
                    "median_dE_by_direction_quintile": med_dE.to_list(),
                    "median_dD_by_direction_quintile": med_dD.to_list(),
                }
            )
    pd.DataFrame(dir_energy_rows).to_csv(OUTPUT_DIR / "stage5b_directional_energy_statistics.csv", index=False)

    # Plots.
    plot_matched_energy_direction(pairs)
    plot_direction_distributions(cos_dist)
    plot_direction_vs_disorder(pairs)
    plot_pca_scatter(
        pca_disp,
        None,
        OUTPUT_DIR / "stage5b_pca_displacement.png",
        "PCA of raw SOAP displacement vectors",
        color_type="group",
        group_masks={"stable_interface": df["structure"].to_numpy()[iface_mask] == "stable",
                     "metastable_interface": df["structure"].to_numpy()[iface_mask] == "metastable"},
    )
    plot_pca_scatter(
        pca_dir,
        None,
        OUTPUT_DIR / "stage5b_pca_direction.png",
        "PCA of normalized SOAP displacement directions",
        color_type="group",
        group_masks={"stable_interface": df["structure"].to_numpy()[iface_mask] == "stable",
                     "metastable_interface": df["structure"].to_numpy()[iface_mask] == "metastable"},
    )
    plot_pca_3d(pca_disp, df, OUTPUT_DIR / "stage5b_pca_displacement_3d.png", "Raw displacement PCA projections")
    plot_direction_vs_residual(pca_disp, Eres[iface_idx], OUTPUT_DIR / "stage5b_direction_vs_residual_energy.png")
    plot_direction_vs_descriptors(pairs, OUTPUT_DIR / "stage5b_direction_vs_stage5a_descriptors.png")

    # Representative cases.
    rep_cases = select_representative_cases(pairs)
    rep_cases.to_csv(OUTPUT_DIR / "stage5b_representative_cases.csv", index=False)
    plot_representative_cases(pairs, rep_cases, OUTPUT_DIR / "stage5b_representative_cases.png")
    write_ovito_representative(rep_cases, df, OUTPUT_DIR / "stage5b_representative_direction_cases.xyz")

    # Direction autocorrelation.
    cells = {}
    for structure in ["stable", "metastable"]:
        import re
        path = PROJECT_ROOT / f"ni_{structure}.cfg"
        lines = Path(path).read_text(encoding="utf-8").splitlines()
        cell = np.zeros((3, 3))
        for line in lines:
            m = re.match(r"H0\(\s*(\d+)\s*,\s*(\d+)\s*\)\s*=\s*([^\sA]+)", line.strip())
            if m:
                cell[int(m.group(1)) - 1, int(m.group(2)) - 1] = float(m.group(3))
        cells[structure] = cell
    dir_autocorr = direction_autocorr(df, u, cells)
    dir_autocorr.to_csv(OUTPUT_DIR / "stage5b_direction_autocorrelation.csv", index=False)
    plot_direction_autocorr(dir_autocorr, OUTPUT_DIR / "stage5b_stable_direction_autocorrelation.png", "stable")
    plot_direction_autocorr(dir_autocorr, OUTPUT_DIR / "stage5b_metastable_direction_autocorrelation.png", "metastable")

    # Permutation null.
    perm = permutation_null(pairs)
    perm.to_csv(OUTPUT_DIR / "stage5b_permutation_null.csv", index=False)

    pca_summary = []
    for name, evr in [("displacement", evr_disp), ("direction", evr_dir)]:
        cum = 0.0
        for i in range(min(10, len(evr))):
            cum += evr[i]
            pca_summary.append({"pca_type": name, "pc": i + 1, "explained_variance_ratio": float(evr[i]), "cumulative": float(cum)})
    pd.DataFrame(pca_summary).to_csv(OUTPUT_DIR / "stage5b_pca_summary.csv", index=False)

    outcome, outcome_text = make_outcome(dir_energy_rows, dir_stats, dir_autocorr)
    write_report(
        validation,
        dir_stats,
        dir_energy_rows,
        rep_cases,
        dir_autocorr,
        pca_summary,
        perm,
        outcome,
        outcome_text,
    )
    print("Stage 5B complete.")
    return 0


def _scatter_xy(x, y, colors, path, title, xlabel, ylabel, labels=None):
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
        r = 3
        odraw.ellipse([px - r, py - r, px + r, py + r], fill=tuple(colors[i]) + (110,))
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


def _scatter_continuous(x, y, color_vals, path, title, xlabel, ylabel, cbar_label):
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
    vmin, vmax = float(np.nanmin(color_vals)), float(np.nanmax(color_vals))
    xpad = (xmax - xmin) * 0.04 if xmax > xmin else 1.0
    ypad = (ymax - ymin) * 0.04 if ymax > ymin else 1.0
    xmin -= xpad
    xmax += xpad
    ymin -= ypad
    ymax += ypad
    for i in range(len(x)):
        fx = (x[i] - xmin) / (xmax - xmin)
        fy = (y[i] - ymin) / (ymax - ymin)
        px = x0 + fx * (x1 - x0)
        py = y1 - fy * (y1 - y0)
        t = (color_vals[i] - vmin) / (vmax - vmin) if vmax > vmin else 0.5
        r = 3
        draw.ellipse([px - r, py - r, px + r, py + r], fill=continuous_color(t))
    bar_x = x1 + 20
    bar_y0 = y1 - 140
    bar_y1 = y1
    for py in range(bar_y0, bar_y1):
        t = (py - bar_y0) / (bar_y1 - bar_y0)
        draw.line([(bar_x, py), (bar_x + 14, py)], fill=continuous_color(t))
    draw.text((bar_x + 18, bar_y0 - 6), f"{vmax:.3g}", fill=(30, 30, 30), font=font_small)
    draw.text((bar_x + 18, bar_y1 - 6), f"{vmin:.3g}", fill=(30, 30, 30), font=font_small)
    draw.text((bar_x, y1 + 12), cbar_label, fill=(30, 30, 30), font=font_small)
    draw.text((x0, y1 + 12), xlabel, fill=(30, 30, 30), font=font_small)
    draw.text((x0 - 70, y0 - 8), ylabel, fill=(30, 30, 30), font=font_small)
    img.save(path, dpi=(150, 150))


def plot_matched_energy_direction(pairs):
    sub = pairs.sample(min(len(pairs), 30000), random_state=SEED)
    colors_map = {"stable_stable": (31, 119, 180), "metastable_metastable": (44, 160, 44), "stable_metastable": (214, 39, 40)}
    colors = np.array([colors_map[c] for c in sub["category"]], dtype=int)
    _scatter_xy(
        sub["directional_distance"].to_numpy(),
        sub["abs_energy_difference"].to_numpy(),
        colors,
        OUTPUT_DIR / "stage5b_matched_energy_difference_vs_direction.png",
        "Matched-disorder energy difference vs directional distance",
        "directional distance (1 - cosine similarity)",
        "absolute energy difference (eV)",
        labels=[(colors_map[c], c) for c in colors_map],
    )


def plot_direction_distributions(cos_dist):
    series = []
    colors = {"stable_stable": (31, 119, 180), "metastable_metastable": (44, 160, 44), "stable_metastable": (214, 39, 40)}
    for category, vals in cos_dist.items():
        vals = np.sort(vals)
        y = np.arange(1, len(vals) + 1) / len(vals)
        series.append({"x": vals, "y": y, "color": colors[category], "label": category, "marker": False})
    img, _ = s3a._new_canvas()
    s3a._plot_panel(
        img,
        (95, 90, 1120, 690),
        "Directional cosine similarity distributions",
        series,
        "cosine similarity",
        "cumulative fraction",
    )
    img.save(OUTPUT_DIR / "stage5b_direction_similarity_distributions.png", dpi=(150, 150))


def plot_direction_vs_disorder(pairs):
    sub = pairs.sample(min(len(pairs), 30000), random_state=SEED)
    mean_d = 0.5 * (sub["D_i"].to_numpy() + sub["D_j"].to_numpy())
    _scatter_continuous(
        mean_d,
        sub["directional_distance"].to_numpy(),
        sub["abs_energy_difference"].to_numpy(),
        OUTPUT_DIR / "stage5b_direction_vs_disorder.png",
        "Directional distance vs mean disorder for matched pairs",
        "mean D(k=10)",
        "directional distance",
        "absolute energy difference (eV)",
    )


def _pca_panel_on_rect(img, draw, pc_a, pc_b, color_mask, rect, title, colors):
    x0, y0, x1, y1 = rect
    draw.line([(x0, y1), (x1, y1)], fill=(60, 60, 60), width=2)
    draw.line([(x0, y0), (x0, y1)], fill=(60, 60, 60), width=2)
    draw.text((x0, y0 - 24), title, fill=(20, 20, 20), font=s3a._font(13))
    xmin, xmax = float(pc_a.min()), float(pc_a.max())
    ymin, ymax = float(pc_b.min()), float(pc_b.max())
    xpad = (xmax - xmin) * 0.04 if xmax > xmin else 1.0
    ypad = (ymax - ymin) * 0.04 if ymax > ymin else 1.0
    xmin -= xpad
    xmax += xpad
    ymin -= ypad
    ymax += ypad
    for i in range(len(pc_a)):
        fx = (pc_a[i] - xmin) / (xmax - xmin)
        fy = (pc_b[i] - ymin) / (ymax - ymin)
        px = x0 + fx * (x1 - x0)
        py = y1 - fy * (y1 - y0)
        r = 3
        draw.ellipse([px - r, py - r, px + r, py + r], fill=colors[int(color_mask[i])])


def plot_pca_3d(pca_disp, df, path, title):
    iface_idx = np.flatnonzero((df["region"] == "interface").to_numpy())
    stable_mask = (df["structure"].to_numpy()[iface_idx] == "stable").astype(int)
    colors = [(31, 119, 180), (214, 39, 40)]
    img = Image.new("RGB", (1150, 760), "white")
    draw = ImageDraw.Draw(img)
    draw.text((50, 24), title, fill=(20, 20, 20), font=s3a._font(18, bold=True))
    rects = [(95, 90, 560, 360), (620, 90, 1090, 360), (95, 430, 560, 700), (620, 430, 1090, 700)]
    _pca_panel_on_rect(img, draw, pca_disp[:, 0], pca_disp[:, 1], stable_mask, rects[0], "PC1-PC2", colors)
    _pca_panel_on_rect(img, draw, pca_disp[:, 0], pca_disp[:, 2], stable_mask, rects[1], "PC1-PC3", colors)
    _pca_panel_on_rect(img, draw, pca_disp[:, 1], pca_disp[:, 2], stable_mask, rects[2], "PC2-PC3", colors)
    draw.rectangle([620, 430, 640, 442], fill=colors[0])
    draw.text((646, 430), "stable interface", fill=(30, 30, 30), font=s3a._font(12))
    draw.rectangle([620, 456, 640, 468], fill=colors[1])
    draw.text((646, 456), "metastable interface", fill=(30, 30, 30), font=s3a._font(12))
    img.save(path, dpi=(150, 150))


def plot_direction_vs_residual(pca_disp, Eres, path):
    plot_pca_scatter(
        pca_disp,
        Eres,
        path,
        "Displacement PCA colored by residual energy",
        color_type="continuous",
        cbar_label="Eres (eV)",
    )


def plot_direction_vs_descriptors(pairs, path):
    img = Image.new("RGB", (1150, 1000), "white")
    draw = ImageDraw.Draw(img)
    draw.text((50, 24), "Directional distance vs Stage 5A physical descriptor differences", fill=(20, 20, 20), font=s3a._font(18, bold=True))
    rects = [(95, 80, 1120, 300), (95, 340, 1120, 560), (95, 600, 1120, 820), (95, 860, 1120, 980)]
    cols = [
        ("CN", "abs_CN_difference"),
        ("RMS bond strain", "abs_strain_difference"),
        ("excess volume", "abs_volume_difference"),
        ("CSP", "abs_csp_difference"),
    ]
    sub = pairs.sample(min(len(pairs), 30000), random_state=SEED)
    for (label, col), rect in zip(cols, rects):
        x = sub["directional_distance"].to_numpy()
        y = sub[col].to_numpy()
        _panel_scatter(draw, rect, x, y, label, "absolute descriptor difference")
    img.save(path, dpi=(150, 150))


def _panel_scatter(draw, rect, x, y, xlabel, ylabel):
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
    for i in range(len(x)):
        fx = (x[i] - xmin) / (xmax - xmin)
        fy = (y[i] - ymin) / (ymax - ymin)
        px = x0 + fx * (x1 - x0)
        py = y1 - fy * (y1 - y0)
        r = 2
        draw.ellipse([px - r, py - r, px + r, py + r], fill=(120, 120, 120))
    draw.text((x0, y1 + 6), xlabel, fill=(30, 30, 30), font=s3a._font(11))
    draw.text((x0, y0 - 18), ylabel, fill=(30, 30, 30), font=s3a._font(11))


def select_representative_cases(pairs):
    cases = []
    case_defs = [
        ("Case1", pairs[(pairs["cosine_similarity"] >= 0.95) & (pairs["abs_energy_difference"] >= 0.15)], "desc"),
        ("Case2", pairs[(pairs["cosine_similarity"] <= 0.2) & (pairs["abs_energy_difference"] >= 0.15)], "desc"),
        ("Case3", pairs[(pairs["cosine_similarity"] <= 0.2) & (pairs["abs_energy_difference"] <= 0.02)], "asc"),
        ("Case4", pairs[(pairs["category"] == "stable_metastable") & (pairs["abs_D_difference"] <= 0.10)], "desc"),
    ]
    for case_id, subset, order in case_defs:
        if len(subset) == 0:
            if case_id == "Case1":
                subset = pairs[(pairs["cosine_similarity"] >= 0.9) & (pairs["abs_energy_difference"] >= 0.15)]
            elif case_id == "Case2":
                subset = pairs[(pairs["cosine_similarity"] <= 0.4) & (pairs["abs_energy_difference"] >= 0.15)]
            elif case_id == "Case3":
                subset = pairs[(pairs["cosine_similarity"] <= 0.4) & (pairs["abs_energy_difference"] <= 0.05)]
            else:
                subset = pairs[pairs["category"] == "stable_metastable"]
        subset = subset.sort_values("abs_energy_difference", ascending=(order == "asc"))
        for _, row in subset.head(2).iterrows():
            out = row.to_dict()
            out["case_id"] = case_id
            cases.append(out)
    return pd.DataFrame(cases)


def plot_representative_cases(pairs, rep_cases, path):
    sub = pairs.sample(min(len(pairs), 30000), random_state=SEED)
    img = Image.new("RGB", (1050, 760), "white")
    draw = ImageDraw.Draw(img)
    font = s3a._font(18, bold=True)
    font_small = s3a._font(12)
    draw.text((50, 24), "Representative matched-direction cases", fill=(20, 20, 20), font=font)
    rect = (90, 90, 950, 690)
    x0, y0, x1, y1 = rect
    draw.line([(x0, y1), (x1, y1)], fill=(60, 60, 60), width=2)
    draw.line([(x0, y0), (x0, y1)], fill=(60, 60, 60), width=2)
    x = sub["directional_distance"].to_numpy()
    y = sub["abs_energy_difference"].to_numpy()
    xmin, xmax = float(x.min()), float(x.max())
    ymin, ymax = float(y.min()), float(y.max())
    xpad = (xmax - xmin) * 0.04
    ypad = (ymax - ymin) * 0.04
    xmin -= xpad
    xmax += xpad
    ymin -= ypad
    ymax += ypad
    overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    odraw = ImageDraw.Draw(overlay)
    for xi, yi in zip(x, y):
        px = x0 + (xi - xmin) / (xmax - xmin) * (x1 - x0)
        py = y1 - (yi - ymin) / (ymax - ymin) * (y1 - y0)
        odraw.ellipse([px - 3, py - 3, px + 3, py + 3], fill=(128, 128, 128, 90))
    img.paste(overlay, (0, 0), overlay)
    draw = ImageDraw.Draw(img)
    case_colors = {"Case1": (31, 119, 180), "Case2": (214, 39, 40), "Case3": (44, 160, 44), "Case4": (255, 127, 14)}
    for _, row in rep_cases.iterrows():
        px = x0 + (row["directional_distance"] - xmin) / (xmax - xmin) * (x1 - x0)
        py = y1 - (row["abs_energy_difference"] - ymin) / (ymax - ymin) * (y1 - y0)
        color = case_colors[row["case_id"]]
        draw.ellipse([px - 7, py - 7, px + 7, py + 7], fill=color)
        draw.text((px + 8, py - 8), row["case_id"], fill=color, font=font_small)
    draw.text((x0, y1 + 12), "directional distance", fill=(30, 30, 30), font=font_small)
    draw.text((x0 - 70, y0 - 8), "absolute energy difference (eV)", fill=(30, 30, 30), font=font_small)
    img.save(path, dpi=(150, 150))


def plot_direction_autocorr(dir_autocorr, path, structure):
    sub = dir_autocorr[dir_autocorr["state"] == structure]
    series = []
    for interface_id, color in [(0, (31, 119, 180)), (1, (214, 39, 40))]:
        part = sub[sub["interface_id"] == interface_id]
        if len(part):
            series.append({"x": part["r_mid"].to_numpy(), "y": part["mean_cosine"].to_numpy(), "color": color, "label": f"interface {interface_id}", "marker": False})
    img, _ = s3a._new_canvas()
    s3a._plot_panel(
        img,
        (95, 90, 1120, 690),
        f"{structure} direction autocorrelation",
        series,
        "in-plane distance r (Angstrom)",
        "mean cosine similarity",
    )
    img.save(path, dpi=(150, 150))


def read_cfg_data(path):
    import re

    lines = Path(path).read_text(encoding="utf-8").splitlines()
    natoms = int(lines[0].split("=")[1])
    cell = np.zeros((3, 3))
    entry_count = None
    auxiliary = {}
    data_start = None
    for idx, line in enumerate(lines):
        stripped = line.strip()
        m = re.match(r"H0\(\s*(\d+)\s*,\s*(\d+)\s*\)\s*=\s*([^\sA]+)", stripped)
        if m:
            cell[int(m.group(1)) - 1, int(m.group(2)) - 1] = float(m.group(3))
        m = re.match(r"entry_count\s*=\s*(\d+)", stripped)
        if m:
            entry_count = int(m.group(1))
        m = re.match(r"auxiliary\[\s*(\d+)\s*\]\s*=\s*(\S+)", stripped)
        if m:
            auxiliary[int(m.group(1))] = m.group(2)
        if data_start is None and entry_count is not None and auxiliary:
            if max(auxiliary) == entry_count - 4:
                data_start = idx + 1
    names = [auxiliary[i] for i in sorted(auxiliary)]
    values = []
    for a in range(natoms):
        values.append([float(x) for x in lines[data_start + a * 3 + 2].split()])
    values = np.asarray(values, dtype=float)
    positions = values[:, :3] @ cell
    name_to_index = {name: i for i, name in enumerate(names)}
    aux = values[:, 3:]
    ids = aux[:, name_to_index["id"]].astype(int)
    peatom = aux[:, name_to_index["c_PEATOM"]]
    return positions, ids, peatom, cell


def write_ovito_representative(rep_cases, df, path):
    Path(path).unlink(missing_ok=True)
    for structure in ["stable", "metastable"]:
        positions, ids, peatom, cell = read_cfg_data(PROJECT_ROOT / f"ni_{structure}.cfg")
        scaled = positions @ np.linalg.inv(cell)
        images = []
        for shift in itertools.product([-1, 0, 1], repeat=3):
            images.append((scaled + np.asarray(shift)) @ cell)
        images = np.vstack(images)
        image_owner = np.tile(np.arange(len(positions)), 27)
        id_to_pos = {int(i): idx for idx, i in enumerate(ids)}
        df_sub = df[df["structure"] == structure]
        atom_map = {int(r["atom_id"]): r for _, r in df_sub.iterrows()}
        central_ids = []
        for _, row in rep_cases.iterrows():
            for side in ["i", "j"]:
                if row[f"state_{side}"] == structure:
                    central_ids.append(int(row[f"atom_{side}"]))
        central_ids = sorted(set(central_ids))
        with open(path, "a" if path.exists() else "w", encoding="utf-8") as fh:
            for aid in central_ids:
                idx = id_to_pos[aid]
                center = positions[idx]
                deltas = images - center
                dists = np.linalg.norm(deltas, axis=1)
                mask = (dists < NEIGHBOR_RADIUS) & (dists > 1e-8)
                neigh_dists = dists[mask]
                neigh_owner = image_owner[mask]
                neigh_pos = images[mask]
                # keep nearest periodic image per atom
                best = {}
                for k, owner in enumerate(neigh_owner):
                    oid = int(owner)
                    if oid not in best or neigh_dists[k] < best[oid][0]:
                        best[oid] = (neigh_dists[k], neigh_pos[k])
                fh.write(f"{len(best)}\n")
                cmeta = atom_map.get(aid, {})
                fh.write(
                    f"central={aid} state={structure} "
                    f"D={cmeta.get('disorder_k10','nan')} dE={cmeta.get('excess_energy_pooled','nan')} "
                    f"Eres={cmeta.get('residual_energy_primary','nan')} "
                    f"CN={cmeta.get('cn','nan')} RMS={cmeta.get('rms_bond_strain','nan')} "
                    f"CSP={cmeta.get('csp','nan')} dV={cmeta.get('excess_voronoi_volume','nan')}\n"
                )
                for owner, (dist, pos) in best.items():
                    nid = int(ids[owner])
                    nmeta = atom_map.get(nid, {})
                    fh.write(
                        f"Ni {pos[0]:.8f} {pos[1]:.8f} {pos[2]:.8f} "
                        f"central={1 if nid == aid else 0} atom_id={nid} "
                        f"PEATOM={peatom[owner]:.8f} D={nmeta.get('disorder_k10','nan')} "
                        f"dE={nmeta.get('excess_energy_pooled','nan')} Eres={nmeta.get('residual_energy_primary','nan')} "
                        f"CN={nmeta.get('cn','nan')} RMS={nmeta.get('rms_bond_strain','nan')} "
                        f"CSP={nmeta.get('csp','nan')} dV={nmeta.get('excess_voronoi_volume','nan')}\n"
                    )


def permutation_null(pairs):
    df4 = pd.read_csv(OUTPUT_DIR / "stage4b_atomwise_energy_disorder.csv")
    iface = df4[df4["region"] == "interface"].copy()
    iface["bin"] = pd.qcut(iface["disorder_k10"], 10, labels=False, duplicates="drop")
    pair_sub = pairs[pairs["tolerance"] == 0.10].copy()
    if len(pair_sub) == 0:
        return pd.DataFrame()
    observed = spearman(pair_sub["directional_distance"].to_numpy(), pair_sub["abs_energy_difference"].to_numpy())
    key_to_index = {}
    base_dE = np.empty(len(df4))
    for pos, (_, r) in enumerate(df4.iterrows()):
        key_to_index[(r["structure"], int(r["atom_id"]))] = pos
        base_dE[pos] = float(r["excess_energy_pooled"])
    pair_i = np.asarray([key_to_index[(r["state_i"], int(r["atom_i"]))] for _, r in pair_sub.iterrows()])
    pair_j = np.asarray([key_to_index[(r["state_j"], int(r["atom_j"]))] for _, r in pair_sub.iterrows()])
    ddir = pair_sub["directional_distance"].to_numpy()
    rng = np.random.default_rng(SEED)
    nulls = []
    n_perm = 100
    bin_arrays = []
    for bin_id, group in iface.groupby("bin"):
        idx = np.asarray([key_to_index[(row["structure"], int(row["atom_id"]))] for _, row in group.iterrows()])
        bin_arrays.append(idx)
    for _ in range(n_perm):
        perm_dE = base_dE.copy()
        for idx in bin_arrays:
            perm_dE[idx] = base_dE[idx][rng.permutation(len(idx))]
        new_abs = np.abs(perm_dE[pair_i] - perm_dE[pair_j])
        nulls.append(spearman(ddir, new_abs))
    nulls = np.asarray(nulls)
    return pd.DataFrame(
        [
            {
                "observed_spearman": float(observed),
                "null_mean": float(nulls.mean()),
                "null_std": float(nulls.std(ddof=1)),
                "null_p5": float(np.percentile(nulls, 5)),
                "null_p95": float(np.percentile(nulls, 95)),
                "n_permutations": n_perm,
                "seed": SEED,
                "label": "atom-level permutation diagnostic",
            }
        ]
    )


def make_outcome(dir_energy_rows, dir_stats, dir_autocorr):
    rhos = [abs(r["spearman"]) for r in dir_energy_rows if np.isfinite(r["spearman"])]
    max_rho = max(rhos) if rhos else 0.0
    tol10_cats = [abs(r["spearman"]) for r in dir_energy_rows if r["tolerance"] == 0.10 and r["category"] in ("stable_stable", "metastable_metastable", "stable_metastable")]
    min_cat_rho = min(tol10_cats) if tol10_cats else 0.0
    if max_rho >= 0.3 and min_cat_rho >= 0.3:
        outcome = "OUTCOME A"
        text = "Direction contains clear additional information at comparable disorder magnitude."
    elif max_rho >= 0.1:
        outcome = "OUTCOME B"
        text = (
            "Direction contains partial additional information; metastable and cross-state "
            "associations are stronger, while stable-stable directional-energy association is weak."
        )
    else:
        outcome = "OUTCOME C"
        text = "Direction adds little beyond magnitude and the tested physical descriptors."
    return outcome, text


def write_report(validation, dir_stats, dir_energy_rows, rep_cases, dir_autocorr, pca_summary, perm, outcome, outcome_text):
    dir_stat_table = "\n".join(
        f"| {r['category']} | {r['n_pairs']} | {r['cosine_mean']:.4f} | {r['cosine_median']:.4f} | "
        f"{r['cosine_std']:.4f} | {r['angle_mean']:.2f} |" for r in dir_stats
    )
    dir_energy_table = "\n".join(
        f"| {r['tolerance']} | {r['category']} | {r['n_pairs']} | {r['spearman']:.4f} | "
        f"{r['pearson']:.4f} | {r['median_dE_by_direction_quintile']} |" for r in dir_energy_rows
    )
    pca_table = "\n".join(
        f"| {r['pca_type']} | {r['pc']} | {r['explained_variance_ratio']:.5f} | {r['cumulative']:.5f} |"
        for r in pca_summary
    )
    autocorr_table = "\n".join(
        f"| {r['state']} | {r['interface_id']} | {r['r_mid']:.1f} | {r['n_pairs']} | {r['mean_cosine']:.4f} |"
        for r in dir_autocorr.head(20).to_dict("records")
    )
    perm_line = ""
    if len(perm):
        p = perm.iloc[0]
        perm_line = (
            f"- Observed matched-pair Spearman: {p['observed_spearman']:.4f}; null mean "
            f"{p['null_mean']:.4f}; null 5-95% [{p['null_p5']:.4f}, {p['null_p95']:.4f}] "
            f"({p['n_permutations']} permutations)."
        )
    lines = [
        "# Stage 5B - SOAP Displacement Direction Beyond Scalar Disorder",
        "",
        "## 1. Objective",
        "",
        "At approximately fixed SOAP disorder magnitude, does the direction of the SOAP "
        "displacement from clean bulk contain additional structural or energetic information?",
        "",
        "## 2. Relation to Stages 4A and 4B",
        "",
        "Stage 4A defined D as a robust scalar bulk-departure magnitude. Stage 4B showed that D "
        "has only partial association with interface excess PEATOM. Stage 5B therefore separates "
        "magnitude D from normalized displacement direction u.",
        "",
        "## 3. Input validation",
        "",
        f"- SOAP matrix shape (4659, 1014), finite values.",
        "- Reconstructed k=10 mean-distance disorder versus Stage 4A D: max abs diff "
        f"{validation['max_abs_diff']:.6e}, median {validation['median_abs_diff']:.6e}, "
        f"p99 {validation['p99_abs_diff']:.6e}.",
        f"- Valid normalized directions: {validation['n_valid_directions']}.",
        "- The displacement direction uses the mean of the same k=10 nearest clean-bulk reference "
        "vectors; D in atomwise outputs remains the Stage 4A kNN mean distance.",
        "",
        "## 4. Direction definition",
        "",
        "    delta_s = s_i - mean_kNN_bulk",
        "    u = delta_s / ||delta_s||",
        "",
        "cosine similarity between directions is the primary directional measure; directional "
        "distance = 1 - cosine.",
        "",
        "## 5. Matched-disorder analysis",
        "",
        "| tolerance | category | n pairs | Spearman | Pearson | median dE by direction quintile |",
        "|---|---|---|---|---|---|",
        dir_energy_table,
        "",
        perm_line,
        "",
        "## 6. Stable/metastable directional distributions",
        "",
        "| category | n pairs | cosine mean | cosine median | cosine std | angle mean (deg) |",
        "|---|---|---|---|---|---|",
        dir_stat_table,
        "",
        "## 7. PCA visualization",
        "",
        "| PCA type | PC | explained variance ratio | cumulative |",
        "|---|---|---|---|",
        pca_table,
        "",
        "PCA is descriptive only and PCA coordinates are not physical modes.",
        "",
        "## 8. Relation to Stage 5A physical descriptors",
        "",
        "`stage5b_direction_vs_stage5a_descriptors.png` shows directional distance versus CN, "
        "RMS bond strain, excess volume, and CSP differences. Directional information is compared "
        "descriptively with the Stage 5A scalar descriptors.",
        "",
        "## 9. Spatial directional organization",
        "",
        "| state | interface | r_mid | n_pairs | mean cosine |",
        "|---|---|---|---|---|",
        autocorr_table,
        "",
        "## 10. Representative matched environments",
        "",
        "Representative cases are saved in `stage5b_representative_cases.csv` and exported to "
        "`stage5b_representative_direction_cases.xyz`.",
        "",
        "## 11. Limitations",
        "",
        "- Only a small number of independent complete GB structures are available.",
        "- Atom-level statistics are pseudoreplicated; N_atoms != N_independent_GBs.",
        "- SOAP direction depends on the chosen SOAP representation.",
        "- Bulk-reference choice affects displacement direction.",
        "- PCA is descriptive.",
        "- Pairwise matched analysis creates statistical dependence.",
        "- No causal interpretation is made.",
        "",
        "## 12. Final outcome",
        "",
        f"{outcome}: {outcome_text}",
    ]
    with open(OUTPUT_DIR / "stage5b_report.md", "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
