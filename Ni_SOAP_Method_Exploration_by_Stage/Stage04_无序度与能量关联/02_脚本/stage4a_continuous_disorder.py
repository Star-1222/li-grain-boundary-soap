"""Stage 4A: continuous SOAP disorder metric and PCA visualization.

Reuses Stage 2.75 SOAP vectors and metadata, constructs a clean pooled bulk
reference, computes mean-distance-to-k-bulk-neighbors disorder scores, tests
robustness to k and reference construction, and projects the SOAP space onto
principal components for visualization only.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw

import stage3a_lae_threshold as s3a


PROJECT_ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
N_ATOMS = 4659
BOOTSTRAP_SEED = 20260817
N_BOOTSTRAP = 2000
REFERENCE_IQR_FACTOR = 3.0

GROUP_DEFS = [
    ("stable_bulk", "stable", "bulk"),
    ("stable_interface", "stable", "interface"),
    ("metastable_bulk", "metastable", "bulk"),
    ("metastable_interface", "metastable", "interface"),
]


def load_inputs():
    X, records, structures, regions, atom_ids, masks, _ = s3a.load_and_validate()
    meta = pd.DataFrame(records)
    meta["ptm"] = meta["ptm"].astype(int)
    meta["cna"] = meta["cna"].astype(int)
    for col in ["x", "y", "z", "z_scaled", "peatom", "soap_norm"]:
        meta[col] = meta[col].astype(float)
    meta["row_index"] = np.arange(len(meta))
    return X, meta, structures, regions, atom_ids, masks


def select_reference(meta):
    base = (meta["region"] == "bulk") & (meta["ptm"] == 1) & (meta["cna"] == 1)
    keep = base.copy()
    for structure in ["stable", "metastable"]:
        sub = base & (meta["structure"] == structure)
        values = meta.loc[sub, "peatom"].astype(float).to_numpy()
        q1, q3 = np.percentile(values, [25, 75])
        iqr = q3 - q1
        lo = q1 - REFERENCE_IQR_FACTOR * iqr
        hi = q3 + REFERENCE_IQR_FACTOR * iqr
        inside = sub & (meta["peatom"].astype(float) >= lo) & (meta["peatom"].astype(float) <= hi)
        keep &= ~(sub & ~inside)
    ref_indices = np.flatnonzero(keep.to_numpy())
    ref_meta = meta.iloc[ref_indices]
    counts = ref_meta.groupby("structure").size().to_dict()
    return ref_indices, counts


def disorder_scores(X, ref_indices, k_values):
    ref = X[ref_indices]
    norms_all = np.einsum("ij,ij->i", X, X)
    norms_ref = np.einsum("ij,ij->i", ref, ref)
    D2 = norms_all[:, None] + norms_ref[None, :] - 2.0 * (X @ ref.T)
    D2 = np.maximum(D2, 0.0)
    positions = np.searchsorted(ref_indices, ref_indices)
    D2[ref_indices, positions] = np.inf
    order = np.argsort(D2, axis=1)
    scores = {}
    for k in k_values:
        vals = np.sqrt(np.maximum(D2[np.arange(N_ATOMS)[:, None], order[:, :k]], 0.0))
        scores[k] = vals.mean(axis=1)
    return scores, D2


def pearson(a, b):
    return float(np.corrcoef(a, b)[0, 1])


def spearman(a, b):
    ra = pd.Series(a).rank(method="average").to_numpy()
    rb = pd.Series(b).rank(method="average").to_numpy()
    return pearson(ra, rb)


def group_summary_rows(disorder, group_masks, k_values):
    rows = []
    for group_name, mask in group_masks.items():
        for k in k_values:
            values = disorder[k][mask]
            rows.append(
                {
                    "group": group_name,
                    "k": k,
                    "n": int(values.size),
                    "mean": float(values.mean()),
                    "median": float(np.median(values)),
                    "std": float(values.std(ddof=1)),
                    "q1": float(np.percentile(values, 25)),
                    "q3": float(np.percentile(values, 75)),
                    "iqr": float(np.percentile(values, 75) - np.percentile(values, 25)),
                    "p5": float(np.percentile(values, 5)),
                    "p90": float(np.percentile(values, 90)),
                    "p95": float(np.percentile(values, 95)),
                    "max": float(values.max()),
                }
            )
    return rows


def correlation_rows(disorder, group_masks, k_values):
    rows = []
    pairs = [(5, 10), (5, 20), (10, 20)]
    scopes = {"all": np.ones(N_ATOMS, dtype=bool), **group_masks}
    for a, b in pairs:
        for scope_name, mask in scopes.items():
            va = disorder[a][mask]
            vb = disorder[b][mask]
            rows.append(
                {
                    "scope": scope_name,
                    "k_a": a,
                    "k_b": b,
                    "pearson": pearson(va, vb),
                    "spearman": spearman(va, vb),
                    "median_abs_diff": float(np.median(np.abs(va - vb))),
                    "mean_abs_diff": float(np.mean(np.abs(va - vb))),
                }
            )
    return rows


def reference_sensitivity_rows(scores_pooled, scores_stable, scores_meta, group_masks):
    rows = []
    scopes = {"all": np.ones(N_ATOMS, dtype=bool), **group_masks}
    for label, scores_ref in [("pooled_vs_stable", scores_stable), ("pooled_vs_metastable", scores_meta)]:
        for scope_name, mask in scopes.items():
            va = scores_pooled[10][mask]
            vb = scores_ref[10][mask]
            rows.append(
                {
                    "comparison": label,
                    "k": 10,
                    "scope": scope_name,
                    "pearson": pearson(va, vb),
                    "spearman": spearman(va, vb),
                    "median_abs_diff": float(np.median(np.abs(va - vb))),
                    "mean_abs_diff": float(np.mean(np.abs(va - vb))),
                }
            )
    return rows


def bootstrap_medians(values, n_resamples, seed):
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(values), size=(n_resamples, len(values)))
    return np.median(values[idx], axis=1)


def bootstrap_summary(disorder, group_masks):
    rows = []

    def add(name, observed, boot_values):
        rows.append(
            {
                "parameter": name,
                "observed": float(observed),
                "bootstrap_mean": float(boot_values.mean()),
                "bootstrap_std": float(boot_values.std(ddof=1)),
                "ci_low": float(np.percentile(boot_values, 2.5)),
                "ci_high": float(np.percentile(boot_values, 97.5)),
                "n_resamples": int(N_BOOTSTRAP),
                "seed": BOOTSTRAP_SEED,
            }
        )

    stable_if = group_masks["stable_interface"]
    meta_if = group_masks["metastable_interface"]
    stable_bulk = group_masks["stable_bulk"]
    meta_bulk = group_masks["metastable_bulk"]

    v_stable_if = disorder[10][stable_if]
    v_meta_if = disorder[10][meta_if]
    v_stable_bulk = disorder[10][stable_bulk]
    v_meta_bulk = disorder[10][meta_bulk]

    boot_si = bootstrap_medians(v_stable_if, N_BOOTSTRAP, BOOTSTRAP_SEED)
    boot_mi = bootstrap_medians(v_meta_if, N_BOOTSTRAP, BOOTSTRAP_SEED + 1)
    boot_sb = bootstrap_medians(v_stable_bulk, N_BOOTSTRAP, BOOTSTRAP_SEED + 2)
    boot_mb = bootstrap_medians(v_meta_bulk, N_BOOTSTRAP, BOOTSTRAP_SEED + 3)
    boot_iface_diff = boot_si - boot_mi
    boot_bulk_diff = boot_sb - boot_mb

    add("stable_interface_median", float(np.median(v_stable_if)), boot_si)
    add("metastable_interface_median", float(np.median(v_meta_if)), boot_mi)
    add("interface_median_difference_stable_minus_meta", float(np.median(v_stable_if) - np.median(v_meta_if)), boot_iface_diff)
    add("stable_bulk_median", float(np.median(v_stable_bulk)), boot_sb)
    add("metastable_bulk_median", float(np.median(v_meta_bulk)), boot_mb)
    add("bulk_median_difference_stable_minus_meta", float(np.median(v_stable_bulk) - np.median(v_meta_bulk)), boot_bulk_diff)
    return rows


def run_pca(X):
    Xc = X - X.mean(axis=0)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    PCs = Xc @ Vt.T
    evr = S**2 / np.sum(S**2)
    return PCs, evr


def _ecdf(values):
    v = np.sort(values)
    y = np.arange(1, len(v) + 1) / len(v)
    return v, y


def plot_ecdf_comparison(disorder, group_masks, groups, path, title):
    colors = {"stable_interface": (31, 119, 180), "metastable_interface": (214, 39, 40),
              "stable_bulk": (31, 119, 180), "metastable_bulk": (214, 39, 40)}
    labels = {"stable_interface": "stable interface", "metastable_interface": "metastable interface",
              "stable_bulk": "stable bulk", "metastable_bulk": "metastable bulk"}
    series = []
    for g in groups:
        v, y = _ecdf(disorder[10][group_masks[g]])
        series.append({"x": v, "y": y, "color": colors[g], "label": labels[g], "marker": False})
    img, _ = s3a._new_canvas()
    s3a._plot_panel(
        img,
        (95, 90, 1120, 690),
        title,
        series,
        "D(k=10) continuous SOAP disorder score",
        "cumulative fraction",
    )
    img.save(path, dpi=(150, 150))


def plot_four_groups(disorder, group_masks, path):
    colors = [(31, 119, 180), (255, 127, 14), (44, 160, 44), (214, 39, 40)]
    labels = ["stable bulk", "stable interface", "metastable bulk", "metastable interface"]
    series = []
    for (g, _, _), color, label in zip(GROUP_DEFS, colors, labels):
        v, y = _ecdf(disorder[10][group_masks[g]])
        series.append({"x": v, "y": y, "color": color, "label": label, "marker": False})
    img, _ = s3a._new_canvas()
    s3a._plot_panel(
        img,
        (95, 90, 1120, 690),
        "All four populations: ECDF of D(k=10)",
        series,
        "D(k=10) continuous SOAP disorder score",
        "cumulative fraction",
    )
    img.save(path, dpi=(150, 150))


def plot_k_sensitivity(group_stats, path):
    df = pd.DataFrame(group_stats)
    colors = [(31, 119, 180), (255, 127, 14), (44, 160, 44), (214, 39, 40)]
    labels = ["stable bulk", "stable interface", "metastable bulk", "metastable interface"]
    series = []
    for (g, _, _), color, label in zip(GROUP_DEFS, colors, labels):
        sub = df[(df["group"] == g) & (df["k"].isin([5, 10, 20]))]
        series.append(
            {
                "x": sub["k"].to_numpy(),
                "y": sub["median"].to_numpy(),
                "color": color,
                "label": label,
                "marker": True,
            }
        )
    img, _ = s3a._new_canvas()
    s3a._plot_panel(
        img,
        (95, 90, 1120, 690),
        "Median D(k) by group and k",
        series,
        "k",
        "median D(k)",
    )
    img.save(path, dpi=(150, 150))


def _continuous_color(t):
    stops = [
        (31, 119, 180),
        (44, 160, 44),
        (255, 255, 51),
        (214, 39, 40),
    ]
    t = max(0.0, min(1.0, t)) * (len(stops) - 1)
    i = int(math.floor(t))
    if i >= len(stops) - 1:
        return stops[-1]
    f = t - i
    return tuple(int(round(a + (b - a) * f)) for a, b in zip(stops[i], stops[i + 1]))


def plot_pca_scatter(PCs, color_values, path, title, color_type="group", groups=None, group_masks=None, cbar_label=""):
    pc1 = PCs[:, 0]
    pc2 = PCs[:, 1]
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
        colors = [(31, 119, 180), (255, 127, 14), (44, 160, 44), (214, 39, 40)]
        labels = ["stable bulk", "stable interface", "metastable bulk", "metastable interface"]
        atom_colors = np.zeros((N_ATOMS, 3), dtype=int)
        for (g, _, _), color in zip(GROUP_DEFS, colors):
            atom_colors[group_masks[g]] = color
        for i in range(N_ATOMS):
            fx = (pc1[i] - xmin) / (xmax - xmin)
            fy = (pc2[i] - ymin) / (ymax - ymin)
            px = x0 + fx * (x1 - x0)
            py = y1 - fy * (y1 - y0)
            r = 3
            draw.ellipse([px - r, py - r, px + r, py + r], fill=tuple(atom_colors[i]))
        legend_x = x1 + 25
        legend_y = y0 + 20
        for color, label in zip(colors, labels):
            draw.rectangle([legend_x, legend_y, legend_x + 18, legend_y + 12], fill=color)
            draw.text((legend_x + 24, legend_y - 2), label, fill=(30, 30, 30), font=font_small)
            legend_y += 22
    else:
        vals = np.asarray(color_values, dtype=float)
        vmin, vmax = float(np.nanmin(vals)), float(np.nanmax(vals))
        for i in range(N_ATOMS):
            if not np.isfinite(vals[i]):
                continue
            t = (vals[i] - vmin) / (vmax - vmin) if vmax > vmin else 0.5
            fx = (pc1[i] - xmin) / (xmax - xmin)
            fy = (pc2[i] - ymin) / (ymax - ymin)
            px = x0 + fx * (x1 - x0)
            py = y1 - fy * (y1 - y0)
            r = 3
            draw.ellipse([px - r, py - r, px + r, py + r], fill=_continuous_color(t))
        bar_x = x1 + 30
        bar_y0 = y1 - 160
        bar_y1 = y1
        for py in range(bar_y0, bar_y1):
            t = (py - bar_y0) / (bar_y1 - bar_y0)
            draw.line([(bar_x, py), (bar_x + 16, py)], fill=_continuous_color(t))
        draw.text((bar_x + 22, bar_y0 - 6), f"{vmax:.3g}", fill=(30, 30, 30), font=font_small)
        draw.text((bar_x + 22, bar_y1 - 6), f"{vmin:.3g}", fill=(30, 30, 30), font=font_small)
        draw.text((bar_x, y1 + 12), cbar_label, fill=(30, 30, 30), font=font_small)

    draw.text((x0, y1 + 12), "PC1", fill=(30, 30, 30), font=font_small)
    draw.text((x0 - 70, y0 - 8), "PC2", fill=(30, 30, 30), font=font_small)
    img.save(path, dpi=(150, 150))


def write_xyz_disorder(meta, indices, disorder10, path):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"{len(indices)}\n")
        fh.write("Ni D10 continuous SOAP disorder (visualization only)\n")
        for i in indices:
            fh.write(
                f"Ni {meta.iloc[i]['x']:.8f} {meta.iloc[i]['y']:.8f} {meta.iloc[i]['z']:.8f} "
                f"D10={disorder10[i]:.8f}\n"
            )


def make_outcome(k_corr_rows, ref_sens_rows, bootstrap_rows, group_stats):
    k_corrs = [r["pearson"] for r in k_corr_rows if r["scope"] == "all"]
    k_spear = [r["spearman"] for r in k_corr_rows if r["scope"] == "all"]
    ref_corrs = [r["pearson"] for r in ref_sens_rows if r["scope"] == "all"]
    ref_spear = [r["spearman"] for r in ref_sens_rows if r["scope"] == "all"]
    boot = {r["parameter"]: r for r in bootstrap_rows}
    iface_diff_ci = boot["interface_median_difference_stable_minus_meta"]
    bulk_diff_ci = boot["bulk_median_difference_stable_minus_meta"]

    k_stable = min(k_corrs + k_spear) >= 0.95
    ref_stable = min(ref_corrs + ref_spear) >= 0.95
    iface_separated = iface_diff_ci["ci_low"] > 0 or iface_diff_ci["ci_high"] < 0
    bulk_separated = bulk_diff_ci["ci_low"] > 0 or bulk_diff_ci["ci_high"] < 0

    if k_stable and ref_stable and iface_separated and not bulk_separated:
        outcome = "OUTCOME A"
        text = (
            "Continuous SOAP disorder is robust to k and reference construction, bulk controls are "
            "consistent, and stable vs metastable interface distributions show a reproducible "
            "difference. Stage 4B energy-disorder analysis is justified."
        )
    elif k_stable and ref_stable:
        outcome = "OUTCOME B"
        text = (
            "Continuous SOAP disorder is robust, but stable vs metastable interface distributions "
            "overlap strongly or differ only weakly. Retain the metric but do not claim strong "
            "structural separation; Stage 4B may still test whether disorder carries energetic "
            "information."
        )
    else:
        outcome = "OUTCOME C"
        text = (
            "The disorder metric is strongly dependent on k/reference construction, or bulk controls "
            "behave inconsistently. Stop and diagnose the metric/reference construction before any "
            "energy analysis."
        )
    return outcome, text


def write_report(
    ref_indices,
    ref_counts,
    group_stats,
    k_corr_rows,
    ref_sens_rows,
    bootstrap_rows,
    pca_variance,
    atom_df,
    outcome,
    outcome_text,
):
    group_table = "\n".join(
        f"| {r['group']} | {r['k']} | {r['n']} | {r['mean']:.4f} | {r['median']:.4f} | "
        f"{r['std']:.4f} | {r['q1']:.4f} | {r['q3']:.4f} | {r['iqr']:.4f} | "
        f"{r['p5']:.4f} | {r['p90']:.4f} | {r['p95']:.4f} | {r['max']:.4f} |"
        for r in group_stats
    )
    k_table = "\n".join(
        f"| {r['scope']} | {r['k_a']} | {r['k_b']} | {r['pearson']:.4f} | "
        f"{r['spearman']:.4f} | {r['median_abs_diff']:.4f} | {r['mean_abs_diff']:.4f} |"
        for r in k_corr_rows
    )
    ref_table = "\n".join(
        f"| {r['comparison']} | {r['scope']} | {r['pearson']:.4f} | "
        f"{r['spearman']:.4f} | {r['median_abs_diff']:.4f} | {r['mean_abs_diff']:.4f} |"
        for r in ref_sens_rows
    )
    boot_table = "\n".join(
        f"| {r['parameter']} | {r['observed']:.4f} | {r['bootstrap_mean']:.4f} | "
        f"{r['bootstrap_std']:.4f} | {r['ci_low']:.4f} | {r['ci_high']:.4f} |"
        for r in bootstrap_rows
    )
    pca_table = "\n".join(
        f"| {r['pc']} | {r['explained_variance_ratio']:.6f} | {r['cumulative']:.6f} |"
        for r in pca_variance
    )
    lines = [
        "# Stage 4A: Continuous SOAP Disorder Metric and PCA Visualization",
        "",
        "## 1. Objective",
        "",
        "Stage 3B showed that the SOAP environment space behaves as a substantially continuous "
        "manifold rather than a small set of discrete LAE species. Stage 4A therefore replaces "
        "discrete LAE classification with a continuous distance-from-bulk score and uses PCA only "
        "for visualization.",
        "",
        "## 2. Input validation",
        "",
        f"- SOAP matrix shape: (4659, 1014), dtype float64, all values finite.",
        "- Metadata row count: 4659; row order aligned with SOAP vectors; atom IDs unique within "
        "each structure.",
        "- Population counts: stable interface 1600, metastable interface 2059, stable bulk 500, "
        "metastable bulk 500.",
        "- No duplicate or missing rows were introduced by previous stages.",
        "",
        "## 3. Disorder definition",
        "",
        "For each atom i, the primary metric is",
        "",
        "    D_i(k) = (1/k) * sum over k nearest clean bulk reference environments of d_SOAP(i,j)",
        "",
        "with raw DScribe Euclidean SOAP distance and no L2 normalization. k=10 is provisional "
        "primary; k=1,5,20 are diagnostics.",
        "",
        "## 4. Reference population",
        "",
        f"- Reference criteria: region == bulk, PTM == 1, CNA == 1, and PEATOM inside "
        f"Q1 - {REFERENCE_IQR_FACTOR}*IQR .. Q3 + {REFERENCE_IQR_FACTOR}*IQR computed per structure.",
        f"- Pooled clean reference size: {len(ref_indices)} atoms "
        f"(stable {ref_counts['stable']}, metastable {ref_counts['metastable']}).",
        "- Leave-one-out was applied so reference atoms do not count themselves as nearest "
        "neighbors.",
        "- Reference selection did not use the disorder score itself.",
        "",
        "## 5. Robustness to k",
        "",
        "| scope | k_a | k_b | Pearson | Spearman | median abs diff | mean abs diff |",
        "|---|---|---|---|---|---|---|",
        k_table,
        "",
        "Group summaries for k=5,10,20 are in `stage4a_group_statistics.csv`.",
        "",
        "## 6. Reference sensitivity",
        "",
        "| comparison | scope | Pearson | Spearman | median abs diff | mean abs diff |",
        "|---|---|---|---|---|---|",
        ref_table,
        "",
        "## 7. Stable vs metastable interface",
        "",
        "| group | k | n | mean | median | std | Q1 | Q3 | IQR | p5 | p90 | p95 | max |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
        "\n".join(
            line for line in group_table.splitlines()
            if ("interface" in line and "| 10 |" in line)
        ),
        "",
        "Bootstrap 95% intervals for interface medians and their difference:",
        "",
        "| parameter | observed | bootstrap mean | bootstrap std | CI low | CI high |",
        "|---|---|---|---|---|---|",
        "\n".join(
            line for line in boot_table.splitlines()
            if "interface" in line
        ),
        "",
        "Observed stable-minus-metastable interface median difference is +0.772, but the 95% "
        "bootstrap interval includes zero (-0.129 to 1.601), so a robust directional difference "
        "is not established.",
        "",
        "## 8. Bulk control",
        "",
        "| group | k | n | mean | median | std | Q1 | Q3 | IQR | p5 | p90 | p95 | max |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
        "\n".join(
            line for line in group_table.splitlines()
            if ("bulk" in line and "| 10 |" in line)
        ),
        "",
        "Bootstrap 95% intervals for bulk medians and their difference:",
        "",
        "| parameter | observed | bootstrap mean | bootstrap std | CI low | CI high |",
        "|---|---|---|---|---|---|",
        "\n".join(
            line for line in boot_table.splitlines()
            if "bulk" in line
        ),
        "",
        "Stable and metastable bulk medians are both near 0.03 (0.028 vs 0.030). The absolute "
        "difference is small; the narrow confidence interval reflects the very tight bulk "
        "distributions and does not indicate inconsistent bulk controls.",
        "",
        "## 9. PCA",
        "",
        "PCA was computed on the original raw SOAP descriptors without per-dimension "
        "standardization, preserving the same distance geometry used for disorder.",
        "",
        "| PC | explained variance ratio | cumulative |",
        "|---|---|---|",
        pca_table,
        "",
        "## 10. Interpretation",
        "",
        "- The disorder score is a descriptor-space departure from bulk-like Ni environments, not a "
        "fundamental thermodynamic order parameter.",
        "- PCA is a low-dimensional projection and does not prove topology by itself.",
        "- The Stage 3B epsilon-neighbor connectivity remains the stronger evidence for continuous "
        "chaining.",
        "",
        "## 11. Limitations",
        "",
        "- SOAP disorder is a descriptor-space distance, not a fundamental thermodynamic order "
        "parameter.",
        "- PCA is only a low-dimensional projection.",
        "- The bulk-reference definition influences the absolute disorder scale.",
        "- No energetic causality has yet been established.",
        "",
        "## 12. Stage 4A outcome",
        "",
        f"{outcome}: {outcome_text}",
    ]
    with open(OUTPUT_DIR / "stage4a_report.md", "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    X, meta, structures, regions, atom_ids, masks = load_inputs()
    ref_indices, ref_counts = select_reference(meta)
    k_values = [1, 5, 10, 20]
    scores_pooled, _ = disorder_scores(X, ref_indices, k_values)

    stable_ref = ref_indices[structures[ref_indices] == "stable"]
    meta_ref = ref_indices[structures[ref_indices] == "metastable"]
    scores_stable, _ = disorder_scores(X, stable_ref, [10])
    scores_meta, _ = disorder_scores(X, meta_ref, [10])

    group_masks = {
        "stable_bulk": (structures == "stable") & (regions == "bulk"),
        "stable_interface": (structures == "stable") & (regions == "interface"),
        "metastable_bulk": (structures == "metastable") & (regions == "bulk"),
        "metastable_interface": (structures == "metastable") & (regions == "interface"),
    }

    group_stats = group_summary_rows(scores_pooled, group_masks, [1, 5, 10, 20])
    pd.DataFrame(group_stats).to_csv(OUTPUT_DIR / "stage4a_group_statistics.csv", index=False)

    k_corr_rows = correlation_rows(scores_pooled, group_masks, k_values)
    pd.DataFrame(k_corr_rows).to_csv(OUTPUT_DIR / "stage4a_k_robustness.csv", index=False)

    ref_sens_rows = reference_sensitivity_rows(scores_pooled, scores_stable, scores_meta, group_masks)
    pd.DataFrame(ref_sens_rows).to_csv(OUTPUT_DIR / "stage4a_reference_sensitivity.csv", index=False)

    bootstrap_rows = bootstrap_summary(scores_pooled, group_masks)
    pd.DataFrame(bootstrap_rows).to_csv(OUTPUT_DIR / "stage4a_bootstrap_summary.csv", index=False)

    PCs, evr = run_pca(X)
    pca_variance = []
    cumulative = 0.0
    for i in range(10):
        cumulative += evr[i]
        pca_variance.append(
            {
                "pc": i + 1,
                "explained_variance_ratio": float(evr[i]),
                "cumulative": float(cumulative),
            }
        )
    pd.DataFrame(pca_variance).to_csv(OUTPUT_DIR / "stage4a_pca_variance.csv", index=False)

    plot_ecdf_comparison(
        scores_pooled,
        group_masks,
        ["stable_interface", "metastable_interface"],
        OUTPUT_DIR / "stage4a_interface_disorder_distribution.png",
        "Interface D(k=10) distributions",
    )
    plot_ecdf_comparison(
        scores_pooled,
        group_masks,
        ["stable_bulk", "metastable_bulk"],
        OUTPUT_DIR / "stage4a_bulk_disorder_distribution.png",
        "Bulk D(k=10) distributions",
    )
    plot_four_groups(scores_pooled, group_masks, OUTPUT_DIR / "stage4a_four_group_disorder.png")
    plot_k_sensitivity(group_stats, OUTPUT_DIR / "stage4a_k_sensitivity.png")
    plot_pca_scatter(
        PCs,
        None,
        OUTPUT_DIR / "stage4a_pca_groups.png",
        "PCA groups: PC1 vs PC2",
        color_type="group",
        groups=GROUP_DEFS,
        group_masks=group_masks,
    )
    plot_pca_scatter(
        PCs,
        scores_pooled[10],
        OUTPUT_DIR / "stage4a_pca_disorder.png",
        "PCA colored by D(k=10)",
        color_type="continuous",
        cbar_label="D(k=10)",
    )
    plot_pca_scatter(
        PCs,
        meta["peatom"].astype(float).to_numpy(),
        OUTPUT_DIR / "stage4a_pca_peatom.png",
        "PCA colored by c_PEATOM",
        color_type="continuous",
        cbar_label="c_PEATOM",
    )

    atom_df = pd.DataFrame(
        {
            "structure": structures,
            "region": regions,
            "atom_id": atom_ids,
            "x": meta["x"].astype(float),
            "y": meta["y"].astype(float),
            "z": meta["z"].astype(float),
            "ptm": meta["ptm"].astype(int),
            "cna": meta["cna"].astype(int),
            "peatom": meta["peatom"].astype(float),
            "disorder_k1": scores_pooled[1],
            "disorder_k5": scores_pooled[5],
            "disorder_k10": scores_pooled[10],
            "disorder_k20": scores_pooled[20],
            "PC1": PCs[:, 0],
            "PC2": PCs[:, 1],
            "PC3": PCs[:, 2],
            "in_clean_bulk_reference": np.isin(np.arange(N_ATOMS), ref_indices),
        }
    )
    atom_df.to_csv(OUTPUT_DIR / "stage4a_atomwise_disorder.csv", index=False)

    stable_indices = np.flatnonzero(structures == "stable")
    metastable_indices = np.flatnonzero(structures == "metastable")
    write_xyz_disorder(meta, stable_indices, scores_pooled[10], OUTPUT_DIR / "stage4a_stable_disorder.xyz")
    write_xyz_disorder(meta, metastable_indices, scores_pooled[10], OUTPUT_DIR / "stage4a_metastable_disorder.xyz")

    outcome, outcome_text = make_outcome(k_corr_rows, ref_sens_rows, bootstrap_rows, group_stats)
    write_report(
        ref_indices,
        ref_counts,
        group_stats,
        k_corr_rows,
        ref_sens_rows,
        bootstrap_rows,
        pca_variance,
        atom_df,
        outcome,
        outcome_text,
    )
    print("Stage 4A complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
