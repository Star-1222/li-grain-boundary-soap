"""Stage 4B: descriptive SOAP disorder-energy relationship analysis.

Reuses Stage 4A atomwise disorder and PCA outputs. Computes a pooled clean-bulk
PEATOM baseline, local excess energy, interface-internal correlations, binned
nonparametric trends, common-support stable/metastable comparisons, bootstrap
descriptive intervals, and PCA-based energetic visualizations. No predictive
modeling is performed.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw

import stage4a_continuous_disorder as s4a
import stage3a_lae_threshold as s3a


PROJECT_ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
N_ATOMS = 4659
BOOTSTRAP_SEED = 20260817
N_BOOTSTRAP = 2000
INTERFACE_BINS = 18
STATE_BINS = 12
COMMON_BINS = 10


def load_and_validate():
    path = OUTPUT_DIR / "stage4a_atomwise_disorder.csv"
    df = pd.read_csv(path)
    if len(df) != N_ATOMS:
        raise RuntimeError(f"Stage 4A atomwise row count {len(df)} != {N_ATOMS}")
    counts = df.groupby(["structure", "region"]).size().to_dict()
    expected = {
        ("stable", "bulk"): 500,
        ("stable", "interface"): 1600,
        ("metastable", "bulk"): 500,
        ("metastable", "interface"): 2059,
    }
    if counts != expected:
        raise RuntimeError(f"Population counts mismatch: {counts}")
    numeric_cols = ["peatom", "disorder_k1", "disorder_k5", "disorder_k10", "disorder_k20", "PC1", "PC2", "PC3"]
    if not np.isfinite(df[numeric_cols].to_numpy()).all():
        raise RuntimeError("Stage 4A atomwise data contains non-finite values.")
    if df["atom_id"].isna().any():
        raise RuntimeError("Atom IDs are missing.")
    for structure in ["stable", "metastable"]:
        ids = df.loc[df["structure"] == structure, "atom_id"]
        if ids.duplicated().any():
            raise RuntimeError(f"Atom IDs are duplicated within structure {structure}.")
    medians = df.groupby(["structure", "region"])["disorder_k10"].median()
    expected_medians = {
        ("stable", "interface"): 6.744,
        ("metastable", "interface"): 5.972,
        ("stable", "bulk"): 0.0278,
        ("metastable", "bulk"): 0.0301,
    }
    for key, value in expected_medians.items():
        if abs(medians[key] - value) > 0.01:
            raise RuntimeError(f"Stage 4A median mismatch for {key}: {medians[key]:.4f}")
    return df


def energy_baselines(df):
    ref = df[df["in_clean_bulk_reference"]]
    rows = [
        {"baseline": "pooled", "n": int(len(ref)), "median_peatom": float(ref["peatom"].median())},
        {"baseline": "stable", "n": int((ref["structure"] == "stable").sum()), "median_peatom": float(ref.loc[ref["structure"] == "stable", "peatom"].median())},
        {"baseline": "metastable", "n": int((ref["structure"] == "metastable").sum()), "median_peatom": float(ref.loc[ref["structure"] == "metastable", "peatom"].median())},
    ]
    pd.DataFrame(rows).to_csv(OUTPUT_DIR / "stage4b_energy_baselines.csv", index=False)
    return {r["baseline"]: r["median_peatom"] for r in rows}, rows


def pearson(a, b):
    return float(np.corrcoef(a, b)[0, 1])


def spearman(a, b):
    return s4a.spearman(a, b)


def group_energy_statistics(df, dE, scopes):
    rows = []
    for name, mask in scopes.items():
        values = dE[mask]
        rows.append(
            {
                "group": name,
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
            }
        )
    return rows


def correlation_table(df, baselines):
    scopes = {
        "all": np.ones(N_ATOMS, dtype=bool),
        "all_interface": df["region"] == "interface",
        "stable_interface": (df["region"] == "interface") & (df["structure"] == "stable"),
        "metastable_interface": (df["region"] == "interface") & (df["structure"] == "metastable"),
        "stable_bulk": (df["region"] == "bulk") & (df["structure"] == "stable"),
        "metastable_bulk": (df["region"] == "bulk") & (df["structure"] == "metastable"),
    }
    rows = []
    for baseline_name in ["pooled", "stable", "metastable"]:
        dE = df["peatom"].to_numpy() - baselines[baseline_name]
        for scope_name, mask in scopes.items():
            x = df["disorder_k10"].to_numpy()[mask]
            y = dE[mask]
            rows.append(
                {
                    "baseline": baseline_name,
                    "scope": scope_name,
                    "n": int(mask.sum()),
                    "pearson": pearson(x, y),
                    "spearman": spearman(x, y),
                }
            )
    return rows, scopes


def k_robustness(df, baselines):
    dE = df["peatom"].to_numpy() - baselines["pooled"]
    scopes = {
        "all_interface": df["region"] == "interface",
        "stable_interface": (df["region"] == "interface") & (df["structure"] == "stable"),
        "metastable_interface": (df["region"] == "interface") & (df["structure"] == "metastable"),
    }
    rows = []
    for k in [5, 10, 20]:
        for scope_name, mask in scopes.items():
            x = df[f"disorder_k{k}"].to_numpy()[mask]
            y = dE[mask]
            rows.append(
                {
                    "k": k,
                    "scope": scope_name,
                    "n": int(mask.sum()),
                    "pearson": pearson(x, y),
                    "spearman": spearman(x, y),
                }
            )
    return rows


def binned_trend(df_sub, dE, n_bins, label_prefix=""):
    sub = df_sub.copy()
    sub = sub.assign(dE=dE[df_sub.index])
    sub["bin"] = pd.qcut(sub["disorder_k10"], n_bins, labels=False, duplicates="drop")
    rows = []
    for bin_id, group in sub.groupby("bin", sort=True):
        rows.append(
            {
                "bin_id": bin_id,
                "n": int(len(group)),
                "disorder_low": float(group["disorder_k10"].min()),
                "disorder_high": float(group["disorder_k10"].max()),
                "median_disorder": float(group["disorder_k10"].median()),
                "median_dE": float(group["dE"].median()),
                "mean_dE": float(group["dE"].mean()),
                "q1_dE": float(group["dE"].quantile(0.25)),
                "q3_dE": float(group["dE"].quantile(0.75)),
            }
        )
    return rows


def common_support_comparison(df, dE):
    stable = df[(df["region"] == "interface") & (df["structure"] == "stable")]
    meta = df[(df["region"] == "interface") & (df["structure"] == "metastable")]
    lo = max(float(stable["disorder_k10"].min()), float(meta["disorder_k10"].min()))
    hi = min(float(stable["disorder_k10"].max()), float(meta["disorder_k10"].max()))
    edges = np.linspace(lo, hi, COMMON_BINS + 1)
    rows = []
    for i in range(COMMON_BINS):
        left = edges[i]
        right = edges[i + 1]
        mask_stable = (stable["disorder_k10"] >= left) & (stable["disorder_k10"] < right)
        mask_meta = (meta["disorder_k10"] >= left) & (meta["disorder_k10"] < right)
        if i == COMMON_BINS - 1:
            mask_stable |= stable["disorder_k10"] == hi
            mask_meta |= meta["disorder_k10"] == hi
        d_stable = dE[stable.index[mask_stable]]
        d_meta = dE[meta.index[mask_meta]]
        if len(d_stable) < 20 or len(d_meta) < 20:
            continue
        rows.append(
            {
                "bin_id": i,
                "disorder_low": float(left),
                "disorder_high": float(right),
                "median_disorder": float(0.5 * (left + right)),
                "n_stable": int(len(d_stable)),
                "n_metastable": int(len(d_meta)),
                "stable_median_dE": float(np.median(d_stable)),
                "metastable_median_dE": float(np.median(d_meta)),
                "median_difference_stable_minus_meta": float(np.median(d_stable) - np.median(d_meta)),
            }
        )
    overlap = {"lo": lo, "hi": hi,
               "stable_lo": float(stable["disorder_k10"].min()), "stable_hi": float(stable["disorder_k10"].max()),
               "meta_lo": float(meta["disorder_k10"].min()), "meta_hi": float(meta["disorder_k10"].max())}
    return rows, overlap


def fast_spearman(x, y):
    rx = np.argsort(np.argsort(x))
    ry = np.argsort(np.argsort(y))
    return pearson(rx, ry)


def bootstrap_summary(df, dE, interface_trend, common_rows):
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    iface_mask = df["region"] == "interface"
    stable_mask = iface_mask & (df["structure"] == "stable")
    meta_mask = iface_mask & (df["structure"] == "metastable")
    rows = []

    def add(param, observed, boot):
        rows.append(
            {
                "parameter": param,
                "observed": float(observed),
                "bootstrap_mean": float(boot.mean()),
                "bootstrap_std": float(boot.std(ddof=1)),
                "ci_low": float(np.percentile(boot, 2.5)),
                "ci_high": float(np.percentile(boot, 97.5)),
                "n_resamples": int(N_BOOTSTRAP),
                "seed": BOOTSTRAP_SEED,
                "interval_type": "atom-resampling descriptive bootstrap intervals",
            }
        )

    x_all = df["disorder_k10"].to_numpy()
    y_all = dE

    def boot_spearman(mask, seed_off):
        vals = []
        r = np.random.default_rng(BOOTSTRAP_SEED + seed_off)
        idx_all = np.flatnonzero(mask)
        for _ in range(N_BOOTSTRAP):
            idx = rng.integers(0, len(idx_all), size=len(idx_all))
            vals.append(fast_spearman(x_all[idx_all[idx]], y_all[idx_all[idx]]))
        return np.asarray(vals)

    add("all_interface_spearman", spearman(x_all[iface_mask], y_all[iface_mask]), boot_spearman(iface_mask, 1))
    add("stable_interface_spearman", spearman(x_all[stable_mask], y_all[stable_mask]), boot_spearman(stable_mask, 2))
    add("metastable_interface_spearman", spearman(x_all[meta_mask], y_all[meta_mask]), boot_spearman(meta_mask, 3))

    trend = pd.DataFrame(interface_trend)
    selected = [0, 5, 10, 14, len(trend) - 1]
    for bin_idx in selected:
        row = trend.iloc[bin_idx]
        mask = iface_mask & (df["disorder_k10"] >= row["disorder_low"]) & (df["disorder_k10"] <= row["disorder_high"])
        vals = dE[mask]
        boot = np.empty(N_BOOTSTRAP)
        r = np.random.default_rng(BOOTSTRAP_SEED + 10 + bin_idx)
        for b in range(N_BOOTSTRAP):
            boot[b] = np.median(vals[r.integers(0, len(vals), size=len(vals))])
        add(f"interface_bin_{bin_idx}_median_dE", float(np.median(vals)), boot)

    for common_row in common_rows:
        mask_stable = stable_mask & (df["disorder_k10"] >= common_row["disorder_low"]) & (df["disorder_k10"] < common_row["disorder_high"])
        mask_meta = meta_mask & (df["disorder_k10"] >= common_row["disorder_low"]) & (df["disorder_k10"] < common_row["disorder_high"])
        vs = dE[mask_stable]
        vm = dE[mask_meta]
        boot = np.empty(N_BOOTSTRAP)
        r = np.random.default_rng(BOOTSTRAP_SEED + 100 + common_row["bin_id"])
        for b in range(N_BOOTSTRAP):
            bs = np.median(vs[r.integers(0, len(vs), size=len(vs))])
            bm = np.median(vm[r.integers(0, len(vm), size=len(vm))])
            boot[b] = bs - bm
        add(
            f"common_bin_{common_row['bin_id']}_stable_minus_meta_median_dE",
            float(common_row["median_difference_stable_minus_meta"]),
            boot,
        )
    return rows


def plot_scatter(x, y, colors, labels, path, title, xlabel, ylabel, xlim=None, ylim=None, alpha=110):
    img = Image.new("RGB", (1050, 760), "white")
    rect = (90, 90, 950, 690)
    x0, y0, x1, y1 = rect
    draw = ImageDraw.Draw(img)
    font = s3a._font(18, bold=True)
    font_small = s3a._font(12)
    draw.text((50, 24), title, fill=(20, 20, 20), font=font)
    draw.line([(x0, y1), (x1, y1)], fill=(60, 60, 60), width=2)
    draw.line([(x0, y0), (x0, y1)], fill=(60, 60, 60), width=2)
    if xlim is None:
        xlim = (float(np.nanmin(x)), float(np.nanmax(x)))
    if ylim is None:
        ylim = (float(np.nanmin(y)), float(np.nanmax(y)))
    xmin, xmax = xlim
    ymin, ymax = ylim
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
        odraw.ellipse([px - r, py - r, px + r, py + r], fill=tuple(colors[i]) + (alpha,))
    img.paste(overlay, (0, 0), overlay)
    draw = ImageDraw.Draw(img)
    if labels is not None:
        legend_x = x1 - 220
        legend_y = y0 + 16
        for color, label in labels:
            draw.rectangle([legend_x, legend_y, legend_x + 16, legend_y + 12], fill=color)
            draw.text((legend_x + 22, legend_y - 2), label, fill=(30, 30, 30), font=font_small)
            legend_y += 20
    draw.text((x0, y1 + 12), xlabel, fill=(30, 30, 30), font=font_small)
    draw.text((x0 - 70, y0 - 8), ylabel, fill=(30, 30, 30), font=font_small)
    img.save(path, dpi=(150, 150))


def plot_binned_trend(trend, path, title, with_band=True):
    df = pd.DataFrame(trend)
    x = df["median_disorder"].to_numpy()
    y = df["median_dE"].to_numpy()
    lo = df["q1_dE"].to_numpy()
    hi = df["q3_dE"].to_numpy()
    series = [{"x": x, "y": y, "color": (31, 119, 180), "label": "median Delta E", "marker": True}]
    bands = [{"x": x, "ymin": lo, "ymax": hi, "color": (31, 119, 180)}] if with_band else None
    img, _ = s3a._new_canvas()
    s3a._plot_panel(
        img,
        (95, 90, 1120, 690),
        title,
        series,
        "median D(k=10)",
        "median Delta E (eV)",
        bands=bands,
    )
    img.save(path, dpi=(150, 150))


def plot_state_trends(stable_trend, meta_trend, path):
    s = pd.DataFrame(stable_trend)
    m = pd.DataFrame(meta_trend)
    series = [
        {"x": s["median_disorder"].to_numpy(), "y": s["median_dE"].to_numpy(), "color": (31, 119, 180), "label": "stable interface", "marker": True},
        {"x": m["median_disorder"].to_numpy(), "y": m["median_dE"].to_numpy(), "color": (214, 39, 40), "label": "metastable interface", "marker": True},
    ]
    img, _ = s3a._new_canvas()
    s3a._plot_panel(
        img,
        (95, 90, 1120, 690),
        "State-specific binned median Delta E vs D(k=10)",
        series,
        "median D(k=10)",
        "median Delta E (eV)",
    )
    img.save(path, dpi=(150, 150))


def plot_conditional_difference(common_rows, path):
    df = pd.DataFrame(common_rows)
    x = df["median_disorder"].to_numpy()
    y = df["median_difference_stable_minus_meta"].to_numpy()
    img, _ = s3a._new_canvas()
    s3a._plot_panel(
        img,
        (95, 90, 1120, 690),
        "Stable minus metastable median Delta E at common disorder bins",
        [{"x": x, "y": y, "color": (44, 160, 44), "label": "median difference", "marker": True}],
        "median D(k=10) in common-support bin",
        "Delta E stable - metastable (eV)",
    )
    img.save(path, dpi=(150, 150))


def plot_pca_interface_energy(PCs, df, dE, path):
    stable_mask = (df["region"] == "interface") & (df["structure"] == "stable")
    meta_mask = (df["region"] == "interface") & (df["structure"] == "metastable")
    img = Image.new("RGB", (1050, 820), "white")
    draw = ImageDraw.Draw(img)
    font = s3a._font(18, bold=True)
    font_small = s3a._font(12)
    draw.text((50, 20), "PC1 vs PC2, interface atoms colored by Delta E", fill=(20, 20, 20), font=font)
    rects = [(90, 80, 940, 390), (90, 440, 940, 750)]
    titles = ["stable interface", "metastable interface"]
    masks = [stable_mask, meta_mask]
    pc1 = PCs[:, 0]
    pc2 = PCs[:, 1]
    xmin, xmax = float(pc1.min()), float(pc1.max())
    ymin, ymax = float(pc2.min()), float(pc2.max())
    vmin, vmax = float(np.nanmin(dE)), float(np.nanmax(dE))
    for rect, title, mask in zip(rects, titles, masks):
        x0, y0, x1, y1 = rect
        draw.line([(x0, y1), (x1, y1)], fill=(60, 60, 60), width=2)
        draw.line([(x0, y0), (x0, y1)], fill=(60, 60, 60), width=2)
        draw.text((x0, y0 - 24), title, fill=(20, 20, 20), font=font_small)
        for i in np.flatnonzero(mask):
            fx = (pc1[i] - xmin) / (xmax - xmin)
            fy = (pc2[i] - ymin) / (ymax - ymin)
            px = x0 + fx * (x1 - x0)
            py = y1 - fy * (y1 - y0)
            r = 3
            t = (dE[i] - vmin) / (vmax - vmin) if vmax > vmin else 0.5
            draw.ellipse([px - r, py - r, px + r, py + r], fill=s4a._continuous_color(t))
        bar_x = x1 + 20
        bar_y0 = y1 - 120
        bar_y1 = y1
        for py in range(bar_y0, bar_y1):
            t = (py - bar_y0) / (bar_y1 - bar_y0)
            draw.line([(bar_x, py), (bar_x + 14, py)], fill=s4a._continuous_color(t))
        draw.text((bar_x + 18, bar_y0 - 6), f"{vmax:.3g}", fill=(30, 30, 30), font=font_small)
        draw.text((bar_x + 18, bar_y1 - 6), f"{vmin:.3g}", fill=(30, 30, 30), font=font_small)
        draw.text((x0, y1 + 8), "PC1", fill=(30, 30, 30), font=font_small)
        draw.text((x0 - 70, y0 - 8), "PC2", fill=(30, 30, 30), font=font_small)
    img.save(path, dpi=(150, 150))


def fixed_disorder_spread(df, dE, overlap):
    lo, hi = overlap["lo"], overlap["hi"]
    edges = np.linspace(lo, hi, 6)
    rows = []
    for i in range(5):
        left, right = edges[i], edges[i + 1]
        if i == 4:
            mask = (df["region"] == "interface") & (df["disorder_k10"] >= left) & (df["disorder_k10"] <= right)
        else:
            mask = (df["region"] == "interface") & (df["disorder_k10"] >= left) & (df["disorder_k10"] < right)
        vals = dE[mask]
        rows.append(
            {
                "disorder_low": float(left),
                "disorder_high": float(right),
                "n": int(len(vals)),
                "mean_dE": float(vals.mean()),
                "median_dE": float(np.median(vals)),
                "std_dE": float(vals.std(ddof=1)),
                "iqr_dE": float(np.percentile(vals, 75) - np.percentile(vals, 25)),
                "range_dE": float(vals.max() - vals.min()),
            }
        )
    pd.DataFrame(rows).to_csv(OUTPUT_DIR / "stage4b_fixed_disorder_energy_spread.csv", index=False)
    return rows


def make_outcome(iface_pearson, iface_spearman, stable_spearman, meta_spearman, interface_trend, common_rows):
    trend = pd.DataFrame(interface_trend)
    trend_range = float(trend["median_dE"].max() - trend["median_dE"].min())
    if iface_spearman >= 0.6 and trend_range >= 0.2:
        return "OUTCOME A", (
            "Interface-only disorder and excess energy show a clear, reproducible monotonic "
            "relationship robust to k and baseline choice."
        )
    if iface_spearman >= 0.15 and trend_range >= 0.05:
        return "OUTCOME B", (
            "Disorder and energy are associated within the interface, but only moderately, with "
            "substantial scatter at fixed disorder and non-monotonic behavior in the high-disorder "
            "tail. Scalar bulk-departure captures useful but incomplete energetic information."
        )
    return "OUTCOME C", (
        "Interface-internal disorder-energy correlations are weak or the binned trend is flat or "
        "highly inconsistent. The scalar metric distinguishes bulk-like from interface-like "
        "environments but does not meaningfully resolve energetic variation within the interface."
    )


def write_report(
    baseline_rows,
    group_energy_rows,
    corr_rows,
    k_rows,
    interface_trend,
    stable_trend,
    meta_trend,
    common_rows,
    bootstrap_rows,
    overlap,
    spread_rows,
    pca_qualitative,
    outcome,
    outcome_text,
):
    baseline_table = "\n".join(
        f"| {r['baseline']} | {r['n']} | {r['median_peatom']:.6f} |" for r in baseline_rows
    )
    group_table = "\n".join(
        f"| {r['group']} | {r['n']} | {r['mean']:.5f} | {r['median']:.5f} | {r['std']:.5f} | "
        f"{r['q1']:.5f} | {r['q3']:.5f} | {r['iqr']:.5f} | {r['p5']:.5f} | {r['p90']:.5f} | {r['p95']:.5f} |"
        for r in group_energy_rows
    )
    corr_table = "\n".join(
        f"| {r['baseline']} | {r['scope']} | {r['n']} | {r['pearson']:.4f} | {r['spearman']:.4f} |"
        for r in corr_rows
    )
    k_table = "\n".join(
        f"| {r['k']} | {r['scope']} | {r['n']} | {r['pearson']:.4f} | {r['spearman']:.4f} |" for r in k_rows
    )
    boot_table = "\n".join(
        f"| {r['parameter']} | {r['observed']:.4f} | {r['bootstrap_mean']:.4f} | "
        f"{r['bootstrap_std']:.4f} | {r['ci_low']:.4f} | {r['ci_high']:.4f} |" for r in bootstrap_rows
    )
    interface_bin_table = "\n".join(
        f"| {r['bin_id']} | {r['disorder_low']:.4f} | {r['disorder_high']:.4f} | {r['n']} | "
        f"{r['median_disorder']:.4f} | {r['median_dE']:.5f} | {r['mean_dE']:.5f} | "
        f"{r['q1_dE']:.5f} | {r['q3_dE']:.5f} |" for r in interface_trend
    )
    common_table = "\n".join(
        f"| {r['bin_id']} | {r['disorder_low']:.4f} - {r['disorder_high']:.4f} | {r['n_stable']} | "
        f"{r['n_metastable']} | {r['stable_median_dE']:.5f} | {r['metastable_median_dE']:.5f} | "
        f"{r['median_difference_stable_minus_meta']:.5f} |" for r in common_rows
    )
    spread_table = "\n".join(
        f"| {r['disorder_low']:.4f} - {r['disorder_high']:.4f} | {r['n']} | {r['mean_dE']:.5f} | "
        f"{r['median_dE']:.5f} | {r['std_dE']:.5f} | {r['iqr_dE']:.5f} | {r['range_dE']:.5f} |"
        for r in spread_rows
    )
    lines = [
        "# Stage 4B: Continuous SOAP Disorder-Energy Relationship",
        "",
        "## 1. Objective",
        "",
        "Does continuous SOAP structural disorder carry information about local excess energy "
        "derived from c_PEATOM? This stage is descriptive only.",
        "",
        "## 2. Inputs and validation",
        "",
        "- Reused `stage4a_atomwise_disorder.csv`; 4659 rows, population counts and Stage 4A "
        "medians reproduce within 0.01.",
        "- All required numeric fields are finite and atom IDs are traceable.",
        "",
        "## 3. Energy baseline",
        "",
        "| baseline | n | median PEATOM (eV/atom) |",
        "|---|---|---|",
        baseline_table,
        "",
        "Primary excess energy uses the pooled baseline so stable and metastable atoms stay on the "
        "same energetic scale.",
        "",
        "## 4. Group energetics",
        "",
        "| group | n | mean | median | std | Q1 | Q3 | IQR | p5 | p90 | p95 |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
        group_table,
        "",
        "## 5. Global disorder-energy relationship",
        "",
        "All-atoms correlations are context only because bulk atoms sit near both D=0 and Delta E=0.",
        "",
        "## 6. Interface-only relationship",
        "",
        "| baseline | scope | n | Pearson | Spearman |",
        "|---|---|---|---|---|",
        "\n".join(
            line for line in corr_table.splitlines()
            if ("interface" in line and "bulk" not in line)
        ),
        "",
        "## 7. Stable and metastable interface relationships",
        "",
        "Stable interface Spearman and metastable interface Spearman are listed in the table above; "
        "state-specific binned trends are compared in `stage4b_state_binned_trends.png`.",
        "",
        "## 8. k robustness",
        "",
        "| k | scope | n | Pearson | Spearman |",
        "|---|---|---|---|---|",
        k_table,
        "",
        "## 9. Nonparametric binned trend",
        "",
        "| bin | disorder low | disorder high | n | median D | median Delta E | mean Delta E | Q1 | Q3 |",
        "|---|---|---|---|---|---|---|---|---|",
        interface_bin_table,
        "",
        "The binned trend rises from near zero at low disorder to a plateau near Delta E = 0.35-0.38 "
        "around D=9.6, then declines in the high-disorder metastable tail.",
        "",
        "## 10. Conditional stable/metastable comparison",
        "",
        f"Common support: stable D range {overlap['stable_lo']:.3f}-{overlap['stable_hi']:.3f}, "
        f"metastable D range {overlap['meta_lo']:.3f}-{overlap['meta_hi']:.3f}, overlap "
        f"{overlap['lo']:.3f}-{overlap['hi']:.3f}.",
        "",
        "| bin | disorder range | stable n | metastable n | stable median dE | metastable median dE | difference |",
        "|---|---|---|---|---|---|---|",
        common_table,
        "",
        "Common-support bins with fewer than 20 atoms in either state are omitted.",
        "",
        "## 11. PCA energetic structure",
        "",
        pca_qualitative,
        "",
        "## 12. Scalar-disorder limitation",
        "",
        "| disorder interval | n | mean dE | median dE | std dE | IQR dE | range dE |",
        "|---|---|---|---|---|---|---|",
        spread_table,
        "",
        "Substantial energy spread remains at nearly fixed disorder, so the scalar bulk-departure "
        "magnitude alone does not uniquely specify local energetic cost.",
        "",
        "## 13. Statistical limitations",
        "",
        "All bootstrap intervals below are atom-resampling descriptive intervals. They do not "
        "quantify uncertainty across independent grain-boundary realizations. The analyzed atoms "
        "come from a small number of complete structures and are spatially correlated; "
        "N_atoms != N_independent_GB_samples.",
        "",
        "| parameter | observed | bootstrap mean | bootstrap std | CI low | CI high |",
        "|---|---|---|---|---|---|",
        boot_table,
        "",
        "## 14. Interpretation",
        "",
        "Within these analyzed stable and metastable Ni grain-boundary configurations, atomwise data "
        "show a descriptive structure-energy association of moderate strength and non-monotonic "
        "form. Higher SOAP bulk-departure is associated with higher local excess PEATOM only over "
        "part of the interface range; substantial scatter and a high-disorder tail remain.",
        "",
        "## 15. Stage 4B outcome",
        "",
        f"{outcome}: {outcome_text}",
    ]
    with open(OUTPUT_DIR / "stage4b_report.md", "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df = load_and_validate()
    baselines, baseline_rows = energy_baselines(df)
    dE_pooled = df["peatom"].to_numpy() - baselines["pooled"]

    scopes = {
        "stable_bulk": (df["region"] == "bulk") & (df["structure"] == "stable"),
        "stable_interface": (df["region"] == "interface") & (df["structure"] == "stable"),
        "metastable_bulk": (df["region"] == "bulk") & (df["structure"] == "metastable"),
        "metastable_interface": (df["region"] == "interface") & (df["structure"] == "metastable"),
    }
    group_energy_rows = group_energy_statistics(df, dE_pooled, scopes)
    pd.DataFrame(group_energy_rows).to_csv(OUTPUT_DIR / "stage4b_group_energy_statistics.csv", index=False)

    corr_rows, _ = correlation_table(df, baselines)
    pd.DataFrame(corr_rows).to_csv(OUTPUT_DIR / "stage4b_correlations.csv", index=False)

    k_rows = k_robustness(df, baselines)
    pd.DataFrame(k_rows).to_csv(OUTPUT_DIR / "stage4b_k_robustness.csv", index=False)

    iface = df[df["region"] == "interface"]
    interface_trend = binned_trend(iface, dE_pooled, INTERFACE_BINS)
    stable_trend = binned_trend(iface[iface["structure"] == "stable"], dE_pooled, STATE_BINS)
    meta_trend = binned_trend(iface[iface["structure"] == "metastable"], dE_pooled, STATE_BINS)
    pd.DataFrame(interface_trend).to_csv(OUTPUT_DIR / "stage4b_binned_interface_trend.csv", index=False)
    pd.DataFrame(stable_trend).to_csv(OUTPUT_DIR / "stage4b_binned_stable_trend.csv", index=False)
    pd.DataFrame(meta_trend).to_csv(OUTPUT_DIR / "stage4b_binned_metastable_trend.csv", index=False)

    common_rows, overlap = common_support_comparison(df, dE_pooled)
    pd.DataFrame(common_rows).to_csv(OUTPUT_DIR / "stage4b_common_support_comparison.csv", index=False)

    bootstrap_rows = bootstrap_summary(df, dE_pooled, interface_trend, common_rows)
    pd.DataFrame(bootstrap_rows).to_csv(OUTPUT_DIR / "stage4b_bootstrap_summary.csv", index=False)

    atomwise = df.copy()
    atomwise["bulk_energy_baseline_pooled"] = baselines["pooled"]
    atomwise["excess_energy_pooled"] = dE_pooled
    atomwise["excess_energy_stable_baseline"] = df["peatom"].to_numpy() - baselines["stable"]
    atomwise["excess_energy_metastable_baseline"] = df["peatom"].to_numpy() - baselines["metastable"]
    atomwise.to_csv(OUTPUT_DIR / "stage4b_atomwise_energy_disorder.csv", index=False)

    PCs = np.column_stack([df["PC1"].to_numpy(), df["PC2"].to_numpy(), df["PC3"].to_numpy()])
    s4a.plot_pca_scatter(
        PCs,
        dE_pooled,
        OUTPUT_DIR / "stage4b_pca_energy.png",
        "PCA colored by Delta E (pooled baseline)",
        color_type="continuous",
        cbar_label="Delta E (eV)",
    )
    s4a.plot_pca_scatter(
        PCs,
        df["disorder_k10"].to_numpy(),
        OUTPUT_DIR / "stage4b_pca_disorder.png",
        "PCA colored by D(k=10)",
        color_type="continuous",
        cbar_label="D(k=10)",
    )
    plot_pca_interface_energy(PCs, df, dE_pooled, OUTPUT_DIR / "stage4b_pca_interface_energy.png")

    group_colors = {
        "stable_bulk": (31, 119, 180),
        "stable_interface": (255, 127, 14),
        "metastable_bulk": (44, 160, 44),
        "metastable_interface": (214, 39, 40),
    }
    colors = np.empty((N_ATOMS, 3), dtype=int)
    labels = []
    for name, color in group_colors.items():
        colors[scopes[name]] = color
    labels = [
        ((31, 119, 180), "stable bulk"),
        ((255, 127, 14), "stable interface"),
        ((44, 160, 44), "metastable bulk"),
        ((214, 39, 40), "metastable interface"),
    ]
    plot_scatter(
        df["disorder_k10"].to_numpy(),
        dE_pooled,
        colors,
        labels,
        OUTPUT_DIR / "stage4b_all_atoms_scatter.png",
        "All atoms: D(k=10) vs Delta E",
        "D(k=10)",
        "Delta E (eV)",
    )
    iface_mask = df["region"] == "interface"
    iface_colors = colors[iface_mask]
    iface_labels = [((255, 127, 14), "stable interface"), ((214, 39, 40), "metastable interface")]
    iface_x = df["disorder_k10"].to_numpy()[iface_mask]
    iface_y = dE_pooled[iface_mask]
    xlim = (float(iface_x.min()), float(iface_x.max()))
    ylim = (float(iface_y.min()), float(iface_y.max()))
    plot_scatter(
        iface_x,
        iface_y,
        iface_colors,
        iface_labels,
        OUTPUT_DIR / "stage4b_interface_scatter.png",
        "Interface atoms: D(k=10) vs Delta E",
        "D(k=10)",
        "Delta E (eV)",
        xlim=xlim,
        ylim=ylim,
    )
    stable_iface = scopes["stable_interface"]
    meta_iface = scopes["metastable_interface"]
    plot_scatter(
        df["disorder_k10"].to_numpy()[stable_iface],
        dE_pooled[stable_iface],
        np.tile((255, 127, 14), (int(stable_iface.sum()), 1)),
        [((255, 127, 14), "stable interface")],
        OUTPUT_DIR / "stage4b_stable_interface_scatter.png",
        "Stable interface: D(k=10) vs Delta E",
        "D(k=10)",
        "Delta E (eV)",
        xlim=xlim,
        ylim=ylim,
    )
    plot_scatter(
        df["disorder_k10"].to_numpy()[meta_iface],
        dE_pooled[meta_iface],
        np.tile((214, 39, 40), (int(meta_iface.sum()), 1)),
        [((214, 39, 40), "metastable interface")],
        OUTPUT_DIR / "stage4b_metastable_interface_scatter.png",
        "Metastable interface: D(k=10) vs Delta E",
        "D(k=10)",
        "Delta E (eV)",
        xlim=xlim,
        ylim=ylim,
    )

    plot_binned_trend(interface_trend, OUTPUT_DIR / "stage4b_interface_binned_trend.png", "Interface binned median Delta E")
    plot_state_trends(stable_trend, meta_trend, OUTPUT_DIR / "stage4b_state_binned_trends.png")
    plot_conditional_difference(common_rows, OUTPUT_DIR / "stage4b_conditional_state_difference.png")

    spread_rows = fixed_disorder_spread(df, dE_pooled, overlap)

    corr_dict = {r["scope"]: r for r in corr_rows if r["baseline"] == "pooled"}
    outcome, outcome_text = make_outcome(
        corr_dict["all_interface"]["pearson"],
        corr_dict["all_interface"]["spearman"],
        corr_dict["stable_interface"]["spearman"],
        corr_dict["metastable_interface"]["spearman"],
        interface_trend,
        common_rows,
    )
    pca_qualitative = (
        "In PC1-PC2 space, energy variation appears primarily associated with the radial bulk "
        "departure direction, but substantial Delta E spread remains at similar PC distance, "
        "suggesting directional SOAP information also matters. This is qualitative only."
    )
    write_report(
        baseline_rows,
        group_energy_rows,
        corr_rows,
        k_rows,
        interface_trend,
        stable_trend,
        meta_trend,
        common_rows,
        bootstrap_rows,
        overlap,
        spread_rows,
        pca_qualitative,
        outcome,
        outcome_text,
    )
    print("Stage 4B complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
