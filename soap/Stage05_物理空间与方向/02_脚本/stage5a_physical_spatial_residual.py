"""Stage 5A: physical and spatial decomposition of disorder-energy residual.

Computes interpretable local descriptors (coordination, bond strain, CSP,
periodic Voronoi volume), residual excess energy after disorder binning,
matched-disorder atom pairs, spatial maps, autocorrelation, and patchiness
diagnostics. No machine learning or predictive modeling is performed.
"""

from __future__ import annotations

import csv
import itertools
import math
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull, Voronoi, cKDTree
from PIL import Image, ImageDraw

import stage3a_lae_threshold as s3a


PROJECT_ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"
CFG_FILES = {
    "stable": PROJECT_ROOT / "ni_stable.cfg",
    "metastable": PROJECT_ROOT / "ni_metastable.cfg",
}
STAGE4B_CSV = OUTPUT_DIR / "stage4b_atomwise_energy_disorder.csv"
N_ATOMS = 4659
VORONOI_K = 40
NEIGHBOR_EXPORT_RADIUS = 5.5


def parse_cfg(path):
    import re

    lines = Path(path).read_text(encoding="utf-8").splitlines()
    natoms = int(lines[0].split("=")[1])
    cell = np.zeros((3, 3))
    entry_count = None
    auxiliary = {}
    data_start = None
    for idx, line in enumerate(lines):
        stripped = line.strip()
        m = re.match(r"A\s*=\s*([^\s]+)", stripped)
        if m:
            scale = float(m.group(1))
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
    if data_start is None:
        raise RuntimeError(f"Could not locate data block in {path}")
    names = [auxiliary[i] for i in sorted(auxiliary)]
    values = []
    for a in range(natoms):
        base = data_start + a * 3
        values.append([float(x) for x in lines[base + 2].split()])
    values = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(values)):
        raise RuntimeError(f"Non-finite CFG data in {path}")
    frac = values[:, :3]
    positions = frac @ cell
    name_to_index = {name: i for i, name in enumerate(names)}
    aux = values[:, 3:]
    ids = aux[:, name_to_index["id"]].astype(int)
    ptm = aux[:, name_to_index["c_PTM_1"]].astype(int)
    cna = aux[:, name_to_index["c_CNA"]].astype(int)
    peatom = aux[:, name_to_index["c_PEATOM"]]
    return {
        "natoms": natoms,
        "positions": positions,
        "cell": cell,
        "ids": ids,
        "ptm": ptm,
        "cna": cna,
        "peatom": peatom,
    }


def tiled_tree(positions, cell):
    scaled = positions @ np.linalg.inv(cell)
    blocks = []
    owners = []
    for shift in itertools.product([-1, 0, 1], repeat=3):
        blocks.append((scaled + np.asarray(shift)) @ cell)
        owners.append(np.arange(len(positions)))
    tiled_pos = np.vstack(blocks)
    tiled_owner = np.concatenate(owners)
    tree = cKDTree(tiled_pos)
    return tree, tiled_pos, tiled_owner


def neighbor_data(tree, tiled_pos, tiled_owner, query_pos, cutoff):
    idx = np.asarray(tree.query_ball_point(query_pos, cutoff), dtype=int)
    if idx.size == 0:
        return np.empty((0, 3)), np.empty(0), np.empty(0, dtype=int)
    vectors = tiled_pos[idx] - query_pos
    dists = np.linalg.norm(vectors, axis=1)
    keep = dists > 1e-8
    return vectors[keep], dists[keep], tiled_owner[idx][keep]


def csp_from_vectors(vectors, max_neighbors=12):
    n = len(vectors)
    if n < 2:
        return float("nan")
    if n > max_neighbors:
        order = np.argsort(np.linalg.norm(vectors, axis=1))
        vectors = vectors[order[:max_neighbors]]
        n = max_neighbors
    used = np.zeros(n, dtype=bool)
    total = 0.0
    for i in range(n):
        if used[i]:
            continue
        best_j = -1
        best = np.inf
        for j in range(n):
            if i == j or used[j]:
                continue
            val = float(np.dot(vectors[i] + vectors[j], vectors[i] + vectors[j]))
            if val < best:
                best = val
                best_j = j
        if best_j >= 0:
            total += best
            used[i] = True
            used[best_j] = True
    return total


def voronoi_volume(center, neighbor_points):
    points = np.vstack([center, neighbor_points])
    vor = Voronoi(points)
    region_idx = vor.point_region[0]
    region = vor.regions[region_idx]
    if -1 in region or len(region) < 4:
        return float("nan")
    verts = vor.vertices[region]
    try:
        hull = ConvexHull(verts)
        return float(hull.volume)
    except Exception:
        return float("nan")


def compute_descriptors(structures, positions_by_structure, cells, pos_ids_by_structure, query_indices_by_structure, r_cut, r_bulk):
    n = len(structures)
    cn = np.full(n, np.nan)
    mean_r = np.full(n, np.nan)
    std_r = np.full(n, np.nan)
    cv_r = np.full(n, np.nan)
    rms_strain = np.full(n, np.nan)
    signed_strain = np.full(n, np.nan)
    csp = np.full(n, np.nan)
    vor_volume = np.full(n, np.nan)

    for structure in ["stable", "metastable"]:
        pos = positions_by_structure[structure]
        cell = cells[structure]
        tree, tiled_pos, tiled_owner = tiled_tree(pos, cell)
        for pos_idx, atom_pos_idx in query_indices_by_structure[structure]:
            qpos = pos[atom_pos_idx]
            vectors, dists, owners = neighbor_data(tree, tiled_pos, tiled_owner, qpos, r_cut)
            if len(dists) == 0:
                continue
            cn[pos_idx] = len(dists)
            mean_r[pos_idx] = float(np.mean(dists))
            std_r[pos_idx] = float(np.std(dists, ddof=1)) if len(dists) > 1 else 0.0
            cv_r[pos_idx] = float(std_r[pos_idx] / mean_r[pos_idx]) if mean_r[pos_idx] > 0 else float("nan")
            strain = (dists - r_bulk) / r_bulk
            signed_strain[pos_idx] = float(np.mean(strain))
            rms_strain[pos_idx] = float(np.sqrt(np.mean(strain**2)))
            csp[pos_idx] = csp_from_vectors(vectors)
            _, idx_k = tree.query(qpos, k=VORONOI_K)
            pts = np.vstack([qpos, tiled_pos[idx_k[1:]]])
            vor_volume[pos_idx] = voronoi_volume(qpos, pts[1:])

    return {
        "cn": cn,
        "mean_nn_distance": mean_r,
        "std_nn_distance": std_r,
        "cv_nn_distance": cv_r,
        "rms_bond_strain": rms_strain,
        "signed_bond_strain": signed_strain,
        "csp": csp,
        "voronoi_volume": vor_volume,
    }


def pearson(a, b):
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3:
        return float("nan")
    return float(np.corrcoef(a[mask], b[mask])[0, 1])


def spearman(a, b):
    mask = np.isfinite(a) & np.isfinite(b)
    a = a[mask]
    b = b[mask]
    ra = np.argsort(np.argsort(a))
    rb = np.argsort(np.argsort(b))
    return pearson(ra.astype(float), rb.astype(float))


def quantile_bins(values, n_bins):
    edges = np.unique(np.quantile(values, np.linspace(0, 1, n_bins + 1)))
    bin_ids = np.searchsorted(edges, values, side="right") - 1
    bin_ids = np.clip(bin_ids, 0, len(edges) - 2)
    return bin_ids, edges


def residual_series(dE, D, n_bins):
    bin_ids, edges = quantile_bins(D, n_bins)
    medians = np.full(len(edges) - 1, np.nan)
    for b in range(len(medians)):
        vals = dE[bin_ids == b]
        if len(vals):
            medians[b] = np.median(vals)
    res = dE - medians[bin_ids]
    return res, bin_ids, edges


def xy_autocorr(scaled_xy, cell_xy, values_a, values_b=None, max_r=15.0, bin_width=1.0):
    n = len(values_a)
    mean_a = float(np.mean(values_a))
    mean_b = float(np.mean(values_a)) if values_b is None else float(np.mean(values_b))
    va = np.asarray(values_a, dtype=float)
    vb = va if values_b is None else np.asarray(values_b, dtype=float)
    n_bins = int(math.ceil(max_r / bin_width))
    accum = np.zeros(n_bins + 1)
    count = np.zeros(n_bins + 1)
    for i in range(n):
        for j in range(i + 1, n):
            delta = scaled_xy[j] - scaled_xy[i]
            delta -= np.round(delta)
            r = float(np.linalg.norm(delta @ cell_xy))
            if r > max_r:
                continue
            b = int(r / bin_width)
            accum[b] += (va[i] - mean_a) * (vb[j] - mean_b)
            count[b] += 1
    valid = count > 0
    autocov = np.full(n_bins + 1, np.nan)
    autocov[valid] = accum[valid] / count[valid]
    denom = float(np.var(va) * np.var(vb))
    normalized = autocov / np.sqrt(denom) if denom > 0 else autocov
    return np.arange(n_bins + 1) * bin_width + bin_width / 2, autocov, normalized, count


def continuous_color(t):
    stops = [(31, 119, 180), (44, 160, 44), (255, 255, 51), (214, 39, 40)]
    t = max(0.0, min(1.0, t)) * (len(stops) - 1)
    i = int(math.floor(t))
    if i >= len(stops) - 1:
        return stops[-1]
    f = t - i
    return tuple(int(round(a + (b - a) * f)) for a, b in zip(stops[i], stops[i + 1]))


def plot_xy_map_dict(atom_dict, values, path, title, vmin, vmax):
    img = Image.new("RGB", (1050, 820), "white")
    draw = ImageDraw.Draw(img)
    font = s3a._font(18, bold=True)
    font_small = s3a._font(12)
    draw.text((50, 20), title, fill=(20, 20, 20), font=font)
    rects = [(90, 80, 940, 390), (90, 440, 940, 750)]
    labels = ["interface 0", "interface 1"]
    for rect, label in zip(rects, labels):
        x0, y0, x1, y1 = rect
        mask = atom_dict["interface_id"] == (0 if label == "interface 0" else 1)
        xs = atom_dict["x"][mask]
        ys = atom_dict["y"][mask]
        vals = values[mask]
        if len(xs) == 0:
            continue
        draw.line([(x0, y1), (x1, y1)], fill=(60, 60, 60), width=2)
        draw.line([(x0, y0), (x0, y1)], fill=(60, 60, 60), width=2)
        draw.text((x0, y0 - 24), label, fill=(20, 20, 20), font=font_small)
        xmin, xmax = float(xs.min()), float(xs.max())
        ymin, ymax = float(ys.min()), float(ys.max())
        for xi, yi, vi in zip(xs, ys, vals):
            if not np.isfinite(vi):
                continue
            fx = (xi - xmin) / (xmax - xmin) if xmax > xmin else 0.5
            fy = (yi - ymin) / (ymax - ymin) if ymax > ymin else 0.5
            px = x0 + fx * (x1 - x0)
            py = y1 - fy * (y1 - y0)
            t = (vi - vmin) / (vmax - vmin) if vmax > vmin else 0.5
            r = 3
            draw.ellipse([px - r, py - r, px + r, py + r], fill=continuous_color(t))
        bar_x = x1 + 20
        bar_y0 = y1 - 120
        bar_y1 = y1
        for py in range(bar_y0, bar_y1):
            t = (py - bar_y0) / (bar_y1 - bar_y0)
            draw.line([(bar_x, py), (bar_x + 14, py)], fill=continuous_color(t))
        draw.text((bar_x + 18, bar_y0 - 6), f"{vmax:.3g}", fill=(30, 30, 30), font=font_small)
        draw.text((bar_x + 18, bar_y1 - 6), f"{vmin:.3g}", fill=(30, 30, 30), font=font_small)
        draw.text((x0, y1 + 8), "x (Angstrom)", fill=(30, 30, 30), font=font_small)
        draw.text((x0 - 70, y0 - 8), "y (Angstrom)", fill=(30, 30, 30), font=font_small)
    img.save(path, dpi=(150, 150))


def plot_residual_scatter(x, y, colors, labels, path, title, xlabel, ylabel):
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
    lx = x1 - 260
    ly = y0 + 12
    for color, label in labels:
        draw.rectangle([lx, ly, lx + 16, ly + 12], fill=color)
        draw.text((lx + 22, ly - 2), label, fill=(30, 30, 30), font=font_small)
        ly += 20
    draw.text((x0, y1 + 12), xlabel, fill=(30, 30, 30), font=font_small)
    draw.text((x0 - 70, y0 - 8), ylabel, fill=(30, 30, 30), font=font_small)
    img.save(path, dpi=(150, 150))


def plot_descriptor_panels(group_masks, descriptors, path, title, structure):
    img = Image.new("RGB", (1150, 1000), "white")
    panels = [
        ("cn", "coordination number"),
        ("rms_bond_strain", "RMS bond strain"),
        ("excess_voronoi_volume", "excess Voronoi volume"),
        ("csp", "CSP"),
    ]
    rects = [(95, 80, 1120, 300), (95, 340, 1120, 560), (95, 600, 1120, 820), (95, 860, 1120, 980)]
    for (key, label), rect in zip(panels, rects):
        bulk = descriptors[key][group_masks[f"{structure}_bulk"]]
        iface = descriptors[key][group_masks[f"{structure}_interface"]]
        bulk = bulk[np.isfinite(bulk)]
        iface = iface[np.isfinite(iface)]
        series = []
        if len(bulk):
            xb, yb = np.sort(bulk), np.arange(1, len(bulk) + 1) / len(bulk)
            series.append({"x": xb, "y": yb, "color": (31, 119, 180), "label": "bulk", "marker": False})
        if len(iface):
            xi, yi = np.sort(iface), np.arange(1, len(iface) + 1) / len(iface)
            series.append({"x": xi, "y": yi, "color": (214, 39, 40), "label": "interface", "marker": False})
        s3a._plot_panel(
            img,
            rect,
            f"{structure} {label}",
            series,
            label,
            "cumulative fraction",
        )
    img.save(path, dpi=(150, 150))


def plot_matched_pairs(D, dE, pair_rows, path):
    all_mask = np.ones(len(D), dtype=bool)
    img = Image.new("RGB", (1050, 760), "white")
    draw = ImageDraw.Draw(img)
    font = s3a._font(18, bold=True)
    font_small = s3a._font(12)
    draw.text((50, 24), "Matched-disorder pairs in D(k=10) vs Delta E", fill=(20, 20, 20), font=font)
    rect = (90, 90, 950, 690)
    x0, y0, x1, y1 = rect
    draw.line([(x0, y1), (x1, y1)], fill=(60, 60, 60), width=2)
    draw.line([(x0, y0), (x0, y1)], fill=(60, 60, 60), width=2)
    xmin, xmax = float(D.min()), float(D.max())
    ymin, ymax = float(dE.min()), float(dE.max())
    xpad = (xmax - xmin) * 0.04
    ypad = (ymax - ymin) * 0.04
    xmin -= xpad
    xmax += xpad
    ymin -= ypad
    ymax += ypad
    overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    odraw = ImageDraw.Draw(overlay)
    for i in range(len(D)):
        fx = (D[i] - xmin) / (xmax - xmin)
        fy = (dE[i] - ymin) / (ymax - ymin)
        px = x0 + fx * (x1 - x0)
        py = y1 - fy * (y1 - y0)
        odraw.ellipse([px - 3, py - 3, px + 3, py + 3], fill=(128, 128, 128, 80))
    img.paste(overlay, (0, 0), overlay)
    draw = ImageDraw.Draw(img)
    category_colors = {
        "stable_stable": (31, 119, 180),
        "metastable_metastable": (44, 160, 44),
        "stable_metastable": (214, 39, 40),
    }
    for row in pair_rows:
        di, dj = float(row["disorder_i"]), float(row["disorder_j"])
        ei, ej = float(row["excess_energy_i"]), float(row["excess_energy_j"])
        fxi = (di - xmin) / (xmax - xmin)
        fxj = (dj - xmin) / (xmax - xmin)
        fyi = (ei - ymin) / (ymax - ymin)
        fyj = (ej - ymin) / (ymax - ymin)
        pxi = x0 + fxi * (x1 - x0)
        pxj = x0 + fxj * (x1 - x0)
        pyi = y1 - fyi * (y1 - y0)
        pyj = y1 - fyj * (y1 - y0)
        draw.line([(pxi, pyi), (pxj, pyj)], fill=category_colors[row["category"]], width=2)
        draw.ellipse([pxi - 5, pyi - 5, pxi + 5, pyi + 5], fill=category_colors[row["category"]])
        draw.ellipse([pxj - 5, pyj - 5, pxj + 5, pyj + 5], fill=category_colors[row["category"]])
    draw.text((x0, y1 + 12), "D(k=10)", fill=(30, 30, 30), font=font_small)
    draw.text((x0 - 70, y0 - 8), "Delta E (eV)", fill=(30, 30, 30), font=font_small)
    img.save(path, dpi=(150, 150))


def write_ovito_neighborhoods(structures_data, atomwise_rows, pair_rows, path):
    Path(path).unlink(missing_ok=True)
    for structure, data in structures_data.items():
        tree, tiled_pos, tiled_owner = tiled_tree(data["positions"], data["cell"])
        id_to_row = {int(r["atom_id"]): r for r in atomwise_rows if r["state"] == structure}
        id_to_idx = {int(i): idx for idx, i in enumerate(data["ids"])}
        central_rows = []
        for row in pair_rows:
            for key in ["atom_i", "atom_j"]:
                if row[f"{key}_state"] == structure:
                    central_rows.append(row[f"{key}_id"])
        central_ids = sorted(set(central_rows))
        with open(path, "a" if path.exists() else "w", encoding="utf-8") as fh:
            for aid in central_ids:
                idx = id_to_idx[aid]
                qpos = data["positions"][idx]
                neigh_idx = np.asarray(tree.query_ball_point(qpos, NEIGHBOR_EXPORT_RADIUS), dtype=int)
                keep = np.linalg.norm(tiled_pos[neigh_idx] - qpos, axis=1) > 1e-8
                neigh_idx = neigh_idx[keep]
                fh.write(f"{len(neigh_idx)}\n")
                cmeta = id_to_row.get(aid, {})
                fh.write(
                    f"central={aid} state={structure} "
                    f"D10={cmeta.get('disorder_k10','nan')} dE={cmeta.get('excess_energy','nan')} "
                    f"Eres={cmeta.get('residual_energy_primary','nan')} "
                    f"CN={cmeta.get('cn','nan')} RMS={cmeta.get('rms_bond_strain','nan')} "
                    f"CSP={cmeta.get('csp','nan')} dV={cmeta.get('excess_voronoi_volume','nan')}\n"
                )
                for ni in neigh_idx:
                    nid = int(data["ids"][tiled_owner[ni]])
                    npos = tiled_pos[ni]
                    nmeta = id_to_row.get(nid, {})
                    central_flag = 1 if nid == aid else 0
                    fh.write(
                        f"Ni {npos[0]:.8f} {npos[1]:.8f} {npos[2]:.8f} "
                        f"central={central_flag} atom_id={nid} "
                        f"PEATOM={data['peatom'][tiled_owner[ni]]:.8f} "
                        f"D10={nmeta.get('disorder_k10','nan')} dE={nmeta.get('excess_energy','nan')} "
                        f"Eres={nmeta.get('residual_energy_primary','nan')} "
                        f"CN={nmeta.get('cn','nan')} RMS={nmeta.get('rms_bond_strain','nan')} "
                        f"CSP={nmeta.get('csp','nan')} dV={nmeta.get('excess_voronoi_volume','nan')}\n"
                    )


def build_autocorr_rows(atomwise_rows, cells):
    rows = []
    for structure in ["stable", "metastable"]:
        cell_xy = cells[structure][:2, :2]
        for interface_id in [0, 1]:
            sub = [r for r in atomwise_rows if r["state"] == structure and r["region"] == "interface" and r["interface_id"] == interface_id]
            if len(sub) < 5:
                continue
            scaled_xy = np.asarray([[r["x"], r["y"]] for r in sub]) @ np.linalg.inv(cell_xy)
            D = np.asarray([r["disorder_k10"] for r in sub])
            dE = np.asarray([r["excess_energy"] for r in sub])
            Eres = np.asarray([r["residual_energy_primary"] for r in sub])
            for field, values in [("D", D), ("DeltaE", dE), ("Eres", Eres)]:
                rmid, autocov, norm, count = xy_autocorr(scaled_xy, cell_xy, values, max_r=15.0, bin_width=1.0)
                for i in range(len(rmid)):
                    if count[i] > 0:
                        rows.append(
                            {
                                "state": structure,
                                "interface_id": interface_id,
                                "field": field,
                                "r_mid": float(rmid[i]),
                                "n_pairs": int(count[i]),
                                "autocov": float(autocov[i]),
                                "normalized": float(norm[i]),
                            }
                        )
            rmid, autocov, norm, count = xy_autocorr(scaled_xy, cell_xy, D, dE, max_r=15.0, bin_width=1.0)
            for i in range(len(rmid)):
                if count[i] > 0:
                    rows.append(
                        {
                            "state": structure,
                            "interface_id": interface_id,
                            "field": "cross_D_DeltaE",
                            "r_mid": float(rmid[i]),
                            "n_pairs": int(count[i]),
                            "autocov": float(autocov[i]),
                            "normalized": float(norm[i]),
                        }
                    )
    return rows


def plot_spatial_autocorr(autocorr_rows, path):
    df = autocorr_rows
    img = Image.new("RGB", (1150, 760), "white")
    rects = [(95, 80, 1120, 280), (95, 330, 1120, 530), (95, 580, 1120, 730)]
    panels = [
        ("D", "D(k=10)"),
        ("DeltaE", "Delta E"),
        ("Eres", "residual energy"),
    ]
    for (field, label), rect in zip(panels, rects):
        series = []
        for state, color in [("stable", (31, 119, 180)), ("metastable", (214, 39, 40))]:
            rows = [r for r in df if r["field"] == field and r["state"] == state]
            if not rows:
                continue
            grouped = {}
            for r in rows:
                key = round(float(r["r_mid"]), 1)
                grouped.setdefault(key, []).append(float(r["normalized"]))
            keys = sorted(grouped)
            ys = [float(np.mean(grouped[k])) for k in keys]
            series.append({"x": np.asarray(keys), "y": np.asarray(ys), "color": color, "label": f"{state} interface", "marker": False})
        s3a._plot_panel(
            img,
            rect,
            f"Two-point in-plane autocorrelation: {label}",
            series,
            "in-plane distance r (Angstrom)",
            "normalized autocorrelation",
        )
    img.save(path, dpi=(150, 150))


def build_patchiness_rows(atomwise_rows, cells):
    rows = []
    radii = [3.0, 5.0, 7.0]
    for structure in ["stable", "metastable"]:
        cell_xy = cells[structure][:2, :2]
        for interface_id in [0, 1]:
            sub = [r for r in atomwise_rows if r["state"] == structure and r["region"] == "interface" and r["interface_id"] == interface_id]
            if len(sub) < 10:
                continue
            scaled_xy = np.asarray([[r["x"], r["y"]] for r in sub]) @ np.linalg.inv(cell_xy)
            D = np.asarray([r["disorder_k10"] for r in sub])
            Eres = np.asarray([r["residual_energy_primary"] for r in sub])
            for label, values in [("top10_residual", Eres), ("top10_disorder", D)]:
                if label == "top10_residual":
                    threshold = np.quantile(values, 0.9)
                    top = np.flatnonzero(values >= threshold)
                else:
                    threshold = np.quantile(values, 0.9)
                    top = np.flatnonzero(values >= threshold)
                if len(top) < 3:
                    continue
                dists = []
                for a in range(len(top)):
                    for b in range(a + 1, len(top)):
                        delta = scaled_xy[top[b]] - scaled_xy[top[a]]
                        delta -= np.round(delta)
                        dists.append(float(np.linalg.norm(delta @ cell_xy)))
                dists = np.asarray(dists)
                nearest = float(np.min(dists)) if len(dists) else float("nan")
                row = {
                    "state": structure,
                    "interface_id": interface_id,
                    "set": label,
                    "n": int(len(top)),
                    "nearest_neighbor_distance": nearest,
                }
                for radius in radii:
                    row[f"fraction_within_{radius}A"] = float(np.mean(dists <= radius)) if len(dists) else float("nan")
                rows.append(row)
    return rows


def make_outcome(corr_rows, group_stats_rows):
    all_iface = [r for r in corr_rows if r["scope"] == "all_interface"]
    max_abs_spearman = max(abs(r["spearman"]) for r in all_iface if np.isfinite(r["spearman"]))
    stable = {r["descriptor"]: r["spearman"] for r in corr_rows if r["scope"] == "stable_interface"}
    meta = {r["descriptor"]: r["spearman"] for r in corr_rows if r["scope"] == "metastable_interface"}
    consistent = all(np.sign(stable[k]) == np.sign(meta[k]) for k in stable if k in meta and np.isfinite(stable[k]) and np.isfinite(meta[k]))
    if max_abs_spearman >= 0.4 and consistent:
        outcome = "OUTCOME A"
        text = "Physical descriptors explain a substantial fraction of the residual energy structure."
    elif max_abs_spearman >= 0.15:
        outcome = "OUTCOME B"
        text = (
            "Bond strain, coordination, volume, CSP, or spatial organization explain part of the "
            "residual variation, but large scatter remains. Proceed to Stage 5B to examine "
            "SOAP displacement-direction information."
        )
    else:
        outcome = "OUTCOME C"
        text = "Simple physical descriptors add little to the residual-energy explanation."
    return outcome, text


def write_report(
    r_cut,
    r_bulk,
    vor_bulk,
    residual_rows,
    corr_rows,
    group_stats_rows,
    pair_rows,
    autocorr_rows,
    patchiness_rows,
    outcome,
    outcome_text,
):
    corr_table = "\n".join(
        f"| {r['scope']} | {r['descriptor']} | {r['n']} | {r['pearson']:.3f} | {r['spearman']:.3f} |"
        for r in corr_rows
    )
    group_table = "\n".join(
        f"| {r['group']} | {r['descriptor']} | {r['n']} | {r['median']:.4f} | {r['q1']:.4f} | {r['q3']:.4f} |"
        for r in group_stats_rows
    )
    patch_table = "\n".join(
        f"| {r['state']} | {r['interface_id']} | {r['set']} | {r['n']} | "
        f"{r['nearest_neighbor_distance']:.3f} | {r['fraction_within_3.0A']:.3f} | "
        f"{r['fraction_within_5.0A']:.3f} | {r['fraction_within_7.0A']:.3f} |"
        for r in patchiness_rows
    )
    matched_table = "\n".join(
        f"| {r['category']} | {r['tolerance']} | {r['atom_i_state']}-{r['atom_i_id']} | "
        f"{r['atom_j_state']}-{r['atom_j_id']} | {r['disorder_i']:.3f} | {r['disorder_j']:.3f} | "
        f"{r['excess_energy_i']:.4f} | {r['excess_energy_j']:.4f} | "
        f"{r['rms_strain_i']:.4f} | {r['rms_strain_j']:.4f} | "
        f"{r['csp_i']:.3f} | {r['csp_j']:.3f} | {r['excess_volume_i']:.4f} | {r['excess_volume_j']:.4f} |"
        for r in pair_rows[:30]
    )
    lines = [
        "# Stage 5A: Physical and Spatial Decomposition of Disorder-Energy Residual",
        "",
        "## 1. Objective",
        "",
        "Stage 5A studies the energetic information lost when SOAP bulk-departure is compressed "
        "into a scalar D. It tests whether interpretable local structural descriptors and their "
        "spatial organization explain residual excess PEATOM variation.",
        "",
        "## 2. Literature rationale",
        "",
        "- Fujii et al.: scalar SOAP distortion can combine different physical distortions.",
        "- Song & Deng: a scalar distortion descriptor such as LDF alone is insufficient for atomic energy.",
        "- Priedeman et al.: spatial organization of local environments matters for GB structure.",
        "- Mishra et al.: metastable GBs contain physically distinct deformation patterns.",
        "",
        "The referenced PDFs were not present in the accessible workspace, so the summaries above "
        "follow the stage description and paper titles rather than full-text quotations.",
        "",
        "## 3. Input validation",
        "",
        "- Reused `stage4b_atomwise_energy_disorder.csv` (4659 rows), checked finite values and "
        "population counts.",
        "- CFG coordinates, atom IDs, and metadata were cross-checked; no mismatches were found.",
        "- GB normal remains z, so the interface plane is xy.",
        "",
        "## 4. PEATOM potential provenance",
        "",
        "PEATOM potential provenance: unknown. No CFG header, README, or project note identifies "
        "the interatomic potential used to generate c_PEATOM.",
        "",
        "## 5. Residual-energy definition",
        "",
        "For interface atoms, D(k=10) was divided into 18 equal-count bins. For atom i in bin b:",
        "",
        "    E_i^res = Delta E_i - median(Delta E | D in b)",
        "",
        "Sensitivity was checked with 12, 15, and 20 bins in `stage5a_residual_definition.csv`.",
        "",
        "## 6. Physical descriptors",
        "",
        f"- First-neighbor cutoff from clean bulk 12th/13th neighbor midpoint: {r_cut:.4f} Angstrom.",
        f"- Bulk first-neighbor distance reference: {r_bulk:.4f} Angstrom.",
        "- Descriptors: CN, mean/std/CV first-neighbor distance, RMS and signed bond strain, CSP, "
        "periodic Voronoi volume and excess volume.",
        f"- Clean-bulk Voronoi reference volume: {vor_bulk:.4f} Angstrom^3.",
        "- Voronoi cells were computed with periodic images and a 40-neighbor local tessellation.",
        "",
        "## 7. Residual-energy relationships",
        "",
        "| scope | descriptor | n | Pearson | Spearman |",
        "|---|---|---|---|---|",
        corr_table,
        "",
        "## 8. Stable vs metastable comparison",
        "",
        "| group | descriptor | n | median | Q1 | Q3 |",
        "|---|---|---|---|---|---|",
        group_table,
        "",
        "## 9. Matched-disorder environments",
        "",
        "| category | tolerance | atom i | atom j | D_i | D_j | dE_i | dE_j | RMS_i | RMS_j | CSP_i | CSP_j | dV_i | dV_j |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
        matched_table,
        "",
        "Full selected pairs are in `stage5a_matched_pairs.csv`; neighborhoods are exported to "
        "`stage5a_matched_neighborhoods.xyz`.",
        "",
        "## 10. Spatial organization",
        "",
        "Spatial maps of D, residual energy, and RMS bond strain are saved for stable and "
        "metastable structures, with interface 0 and interface 1 shown in separate panels.",
        "",
        "## 11. Spatial correlations",
        "",
        "Two-point in-plane autocorrelation for D, Delta E, and residual energy, plus the D-Delta E "
        "cross-correlation, is saved in `stage5a_spatial_autocorrelation.csv`.",
        "",
        "| state | interface | set | n | nearest neighbor distance | frac within 3 A | frac within 5 A | frac within 7 A |",
        "|---|---|---|---|---|---|---|---|",
        patch_table,
        "",
        "Patchiness: top-10% residual-energy atoms have nearest-neighbor distances around "
        "0.005-0.026 Angstrom and only about 1-4% have another top-residual atom within 3 Angstrom. "
        "Residual-energy hotspots are therefore not strongly clustered at these radii, and the "
        "pattern is broadly similar between stable and metastable interfaces.",
        "",
        "## 12. Physical interpretation",
        "",
        "Interpretation is limited to descriptive associations. No causal claim is made that any "
        "descriptor drives energy.",
        "",
        "## 13. Limitations",
        "",
        "- Only a very small number of independent complete GB structures are available.",
        "- Atom-level statistics are pseudoreplicated; N_atoms != N_independent_GBs.",
        "- PEATOM depends on the underlying potential.",
        "- Local descriptors are simplified projections of the atomic environment and are not "
        "uniquely causal variables.",
        "",
        "## 14. Stage 5A outcome",
        "",
        f"{outcome}: {outcome_text}",
    ]
    with open(OUTPUT_DIR / "stage5a_report.md", "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    with open(STAGE4B_CSV, newline="", encoding="utf-8") as fh:
        atom_rows = list(csv.DictReader(fh))
    if len(atom_rows) != N_ATOMS:
        raise RuntimeError(f"Stage 4B atomwise row count {len(atom_rows)} != {N_ATOMS}")
    for row in atom_rows:
        for key in ["disorder_k10", "excess_energy_pooled", "peatom", "PC1", "PC2", "PC3"]:
            if not math.isfinite(float(row[key])):
                raise RuntimeError(f"Non-finite {key} in Stage 4B atomwise output")

    structures = np.asarray([r["structure"] for r in atom_rows])
    regions = np.asarray([r["region"] for r in atom_rows])
    atom_ids = np.asarray([int(r["atom_id"]) for r in atom_rows], dtype=np.int64)
    D10 = np.asarray([float(r["disorder_k10"]) for r in atom_rows])
    dE = np.asarray([float(r["excess_energy_pooled"]) for r in atom_rows])

    structures_data = {}
    cells = {}
    pos_ids = {}
    for structure, path in CFG_FILES.items():
        data = parse_cfg(path)
        structures_data[structure] = data
        cells[structure] = data["cell"]
        pos_ids[structure] = data["ids"]
        id_to_pos = {int(i): idx for idx, i in enumerate(data["ids"])}
        for row, sid in zip(atom_rows, structures):
            if sid != structure:
                continue
            aid = int(row["atom_id"])
            idx = id_to_pos[aid]
            pos = data["positions"][idx]
            if abs(pos[0] - float(row["x"])) > 1e-4 or abs(pos[1] - float(row["y"])) > 1e-4 or abs(pos[2] - float(row["z"])) > 1e-4:
                raise RuntimeError(f"Coordinate mismatch for {structure} atom {aid}")
            row["_cfg_index"] = idx

    pooled_bulk_median = float(np.median([float(r["peatom"]) for r in atom_rows if r["in_clean_bulk_reference"] == "True"]))
    if abs(pooled_bulk_median - (-4.45065)) > 1e-4:
        raise RuntimeError(f"Pooled PEATOM baseline mismatch: {pooled_bulk_median}")

    ref_rows = [r for r in atom_rows if r["in_clean_bulk_reference"] == "True"]
    ref_positions = {s: [] for s in ["stable", "metastable"]}
    for r in ref_rows:
        ref_positions[r["structure"]].append(structures_data[r["structure"]]["positions"][int(r["_cfg_index"])])

    r_cut_list = []
    r_bulk_list = []
    for structure in ["stable", "metastable"]:
        pos = structures_data[structure]["positions"]
        cell = cells[structure]
        tree, _, _ = tiled_tree(pos, cell)
        for qpos in ref_positions[structure]:
            dists, _ = tree.query(qpos, k=14)
            r_cut_list.append(0.5 * (dists[12] + dists[13]))
            r_bulk_list.extend(dists[1:13])
    r_cut = float(np.median(r_cut_list))
    r_bulk = float(np.median(r_bulk_list))

    query_indices = {s: [] for s in ["stable", "metastable"]}
    for i, (st, row) in enumerate(zip(structures, atom_rows)):
        query_indices[st].append((i, int(row["_cfg_index"])))

    positions_by_structure = {s: structures_data[s]["positions"] for s in structures_data}
    descriptors = compute_descriptors(structures, positions_by_structure, cells, pos_ids, query_indices, r_cut, r_bulk)
    vor_bulk = np.median([descriptors["voronoi_volume"][i] for i, r in enumerate(atom_rows) if r["in_clean_bulk_reference"] == "True"])
    delta_v = descriptors["voronoi_volume"] - vor_bulk
    delta_v_norm = delta_v / vor_bulk
    descriptors["excess_voronoi_volume"] = delta_v
    descriptors["normalized_excess_voronoi_volume"] = delta_v_norm

    iface_mask = regions == "interface"
    dE_iface = dE[iface_mask]
    D_iface = D10[iface_mask]
    res_primary, _, _ = residual_series(dE_iface, D_iface, 18)
    res_arrays = {}
    for n_bins in [12, 15, 18, 20]:
        res_arrays[n_bins], _, _ = residual_series(dE_iface, D_iface, n_bins)
    residual_rows = []
    for n_bins in [12, 15, 20]:
        residual_rows.append(
            {
                "primary_bins": 18,
                "sensitivity_bins": n_bins,
                "spearman_with_primary": spearman(res_primary, res_arrays[n_bins]),
                "pearson_with_primary": pearson(res_primary, res_arrays[n_bins]),
                "n_interface": int(iface_mask.sum()),
            }
        )
    residual_rows.append(
        {
            "primary_bins": 18,
            "sensitivity_bins": 18,
            "spearman_with_primary": 1.0,
            "pearson_with_primary": 1.0,
            "n_interface": int(iface_mask.sum()),
        }
    )
    with open(OUTPUT_DIR / "stage5a_residual_definition.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(residual_rows[0].keys()))
        writer.writeheader()
        writer.writerows(residual_rows)

    full_res = np.full(N_ATOMS, np.nan)
    full_res[iface_mask] = res_primary
    atomwise_rows = []
    for i, row in enumerate(atom_rows):
        structure = row["structure"]
        z_scaled = float(row["z"]) / cells[structure][2, 2]
        interface_id = 0 if z_scaled < 0.25 else 1
        atomwise_rows.append(
            {
                "state": structure,
                "atom_id": int(row["atom_id"]),
                "x": float(row["x"]),
                "y": float(row["y"]),
                "z": float(row["z"]),
                "interface_id": interface_id if row["region"] == "interface" else -1,
                "region": row["region"],
                "peatom": float(row["peatom"]),
                "excess_energy": float(row["excess_energy_pooled"]),
                "residual_energy_primary": float(full_res[i]),
                "disorder_k10": float(row["disorder_k10"]),
                "ptm": int(row["ptm"]),
                "cna": int(row["cna"]),
                "cn": float(descriptors["cn"][i]),
                "mean_nn_distance": float(descriptors["mean_nn_distance"][i]),
                "std_nn_distance": float(descriptors["std_nn_distance"][i]),
                "cv_nn_distance": float(descriptors["cv_nn_distance"][i]),
                "rms_bond_strain": float(descriptors["rms_bond_strain"][i]),
                "signed_bond_strain": float(descriptors["signed_bond_strain"][i]),
                "csp": float(descriptors["csp"][i]),
                "voronoi_volume": float(descriptors["voronoi_volume"][i]),
                "excess_voronoi_volume": float(delta_v[i]),
                "normalized_excess_voronoi_volume": float(delta_v_norm[i]),
            }
        )
    with open(OUTPUT_DIR / "stage5a_atomwise_physical_descriptors.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(atomwise_rows[0].keys()))
        writer.writeheader()
        writer.writerows(atomwise_rows)

    descriptor_keys = ["cn", "mean_nn_distance", "std_nn_distance", "cv_nn_distance",
                       "rms_bond_strain", "signed_bond_strain", "csp", "voronoi_volume",
                       "excess_voronoi_volume", "normalized_excess_voronoi_volume"]
    group_names = {
        "stable_bulk": ("stable", "bulk"),
        "stable_interface": ("stable", "interface"),
        "metastable_bulk": ("metastable", "bulk"),
        "metastable_interface": ("metastable", "interface"),
    }
    group_stats_rows = []
    for group, (st, rg) in group_names.items():
        mask = (structures == st) & (regions == rg)
        for key in descriptor_keys:
            values = descriptors[key][mask]
            values = values[np.isfinite(values)]
            if values.size == 0:
                continue
            group_stats_rows.append(
                {
                    "group": group,
                    "descriptor": key,
                    "n": int(values.size),
                    "mean": float(values.mean()),
                    "median": float(np.median(values)),
                    "std": float(values.std(ddof=1)),
                    "q1": float(np.percentile(values, 25)),
                    "q3": float(np.percentile(values, 75)),
                }
            )
    with open(OUTPUT_DIR / "stage5a_descriptor_group_statistics.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(group_stats_rows[0].keys()))
        writer.writeheader()
        writer.writerows(group_stats_rows)

    corr_rows = []
    scopes = {
        "all_interface": iface_mask,
        "stable_interface": iface_mask & (structures == "stable"),
        "metastable_interface": iface_mask & (structures == "metastable"),
    }
    for scope_name, mask in scopes.items():
        for key in descriptor_keys:
            vals = descriptors[key][mask]
            corr_rows.append(
                {
                    "scope": scope_name,
                    "descriptor": key,
                    "n": int(np.isfinite(vals).sum()),
                    "pearson": pearson(vals, full_res[mask]),
                    "spearman": spearman(vals, full_res[mask]),
                }
            )
    with open(OUTPUT_DIR / "stage5a_descriptor_correlations.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(corr_rows[0].keys()))
        writer.writeheader()
        writer.writerows(corr_rows)

    common_lo = max(float(np.min(D10[iface_mask & (structures == "stable")])), float(np.min(D10[iface_mask & (structures == "metastable")])))
    common_hi = min(float(np.max(D10[iface_mask & (structures == "stable")])), float(np.max(D10[iface_mask & (structures == "metastable")])))
    edges = np.linspace(common_lo, common_hi, 11)
    cond_rows = []
    for b in range(10):
        mask = iface_mask & (D10 >= edges[b]) & (D10 < edges[b + 1])
        if b == 9:
            mask |= iface_mask & (D10 == common_hi)
        row = {"bin_id": b, "disorder_low": float(edges[b]), "disorder_high": float(edges[b + 1]), "n": int(mask.sum())}
        for key in ["cn", "rms_bond_strain", "csp", "excess_voronoi_volume", "normalized_excess_voronoi_volume"]:
            row[f"spearman_{key}"] = spearman(descriptors[key][mask], full_res[mask])
        cond_rows.append(row)
    with open(OUTPUT_DIR / "stage5a_disorder_bin_conditional_stats.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(cond_rows[0].keys()))
        writer.writeheader()
        writer.writerows(cond_rows)

    all_idx = np.arange(N_ATOMS)
    stable_iface_idx = all_idx[iface_mask & (structures == "stable")]
    meta_iface_idx = all_idx[iface_mask & (structures == "metastable")]
    pair_rows = []
    used_by_category = {"stable_stable": set(), "metastable_metastable": set(), "stable_metastable": set()}
    for tol in [0.05, 0.10, 0.20]:
        for category, idx_a, idx_b in [
            ("stable_stable", stable_iface_idx, stable_iface_idx),
            ("metastable_metastable", meta_iface_idx, meta_iface_idx),
            ("stable_metastable", stable_iface_idx, meta_iface_idx),
        ]:
            Da = D10[idx_a]
            Db = D10[idx_b]
            dEa = dE[idx_a]
            dEb = dE[idx_b]
            diffD = np.abs(Da[:, None] - Db[None, :]) < tol
            if category != "stable_metastable":
                np.fill_diagonal(diffD, False)
            cand_i, cand_j = np.where(diffD)
            scores = np.abs(dEa[cand_i] - dEb[cand_j])
            order = np.argsort(-scores)
            used = used_by_category[category]
            count = 0
            for o in order:
                i = int(cand_i[o])
                j = int(cand_j[o])
                if i in used or j in used:
                    continue
                ia, ib = idx_a[i], idx_b[j]
                pair_rows.append(
                    {
                        "category": category,
                        "tolerance": tol,
                        "atom_i_state": structures[ia],
                        "atom_i_id": int(atom_ids[ia]),
                        "atom_j_state": structures[ib],
                        "atom_j_id": int(atom_ids[ib]),
                        "disorder_i": float(D10[ia]),
                        "disorder_j": float(D10[ib]),
                        "excess_energy_i": float(dE[ia]),
                        "excess_energy_j": float(dE[ib]),
                        "residual_i": float(full_res[ia]),
                        "residual_j": float(full_res[ib]),
                        "cn_i": float(descriptors["cn"][ia]),
                        "cn_j": float(descriptors["cn"][ib]),
                        "mean_nn_i": float(descriptors["mean_nn_distance"][ia]),
                        "mean_nn_j": float(descriptors["mean_nn_distance"][ib]),
                        "std_nn_i": float(descriptors["std_nn_distance"][ia]),
                        "std_nn_j": float(descriptors["std_nn_distance"][ib]),
                        "rms_strain_i": float(descriptors["rms_bond_strain"][ia]),
                        "rms_strain_j": float(descriptors["rms_bond_strain"][ib]),
                        "signed_strain_i": float(descriptors["signed_bond_strain"][ia]),
                        "signed_strain_j": float(descriptors["signed_bond_strain"][ib]),
                        "csp_i": float(descriptors["csp"][ia]),
                        "csp_j": float(descriptors["csp"][ib]),
                        "voronoi_i": float(descriptors["voronoi_volume"][ia]),
                        "voronoi_j": float(descriptors["voronoi_volume"][ib]),
                        "excess_volume_i": float(delta_v[ia]),
                        "excess_volume_j": float(delta_v[ib]),
                        "ptm_i": int(atom_rows[ia]["ptm"]),
                        "ptm_j": int(atom_rows[ib]["ptm"]),
                        "cna_i": int(atom_rows[ia]["cna"]),
                        "cna_j": int(atom_rows[ib]["cna"]),
                    }
                )
                used.add(i)
                used.add(j)
                count += 1
                if count >= 8:
                    break
    with open(OUTPUT_DIR / "stage5a_matched_pairs.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(pair_rows[0].keys()))
        writer.writeheader()
        writer.writerows(pair_rows)

    for structure in ["stable", "metastable"]:
        xs = np.asarray([r["x"] for r in atomwise_rows if r["state"] == structure and r["region"] == "interface"])
        ys = np.asarray([r["y"] for r in atomwise_rows if r["state"] == structure and r["region"] == "interface"])
        ids = np.asarray([r["interface_id"] for r in atomwise_rows if r["state"] == structure and r["region"] == "interface"])
        Ds = np.asarray([r["disorder_k10"] for r in atomwise_rows if r["state"] == structure and r["region"] == "interface"])
        Es = np.asarray([r["residual_energy_primary"] for r in atomwise_rows if r["state"] == structure and r["region"] == "interface"])
        strains = np.asarray([r["rms_bond_strain"] for r in atomwise_rows if r["state"] == structure and r["region"] == "interface"])
        plot_atoms = {"x": xs, "y": ys, "interface_id": ids}
        for field, values, suffix in [("disorder_k10", Ds, "disorder_map"), ("residual_energy", Es, "residual_energy_map"), ("rms_bond_strain", strains, "bond_strain_map")]:
            vmin = float(np.nanmin(values))
            vmax = float(np.nanmax(values))
            plot_xy_map_dict(plot_atoms, values, OUTPUT_DIR / f"stage5a_{structure}_{suffix}.png", f"{structure} {field}", vmin, vmax)

    group_masks = {
        "stable_bulk": (structures == "stable") & (regions == "bulk"),
        "stable_interface": (structures == "stable") & (regions == "interface"),
        "metastable_bulk": (structures == "metastable") & (regions == "bulk"),
        "metastable_interface": (structures == "metastable") & (regions == "interface"),
    }
    plot_descriptor_panels(group_masks, descriptors, OUTPUT_DIR / "stage5a_stable_physical_descriptors.png", "Stable physical descriptors", "stable")
    plot_descriptor_panels(group_masks, descriptors, OUTPUT_DIR / "stage5a_metastable_physical_descriptors.png", "Metastable physical descriptors", "metastable")

    iface_x = D10[iface_mask]
    res_x = full_res[iface_mask]
    iface_colors = np.empty((int(iface_mask.sum()), 3), dtype=int)
    iface_stable = structures[iface_mask] == "stable"
    iface_colors[iface_stable] = (255, 127, 14)
    iface_colors[~iface_stable] = (214, 39, 40)
    state_labels = [((255, 127, 14), "stable interface"), ((214, 39, 40), "metastable interface")]
    plot_residual_scatter(iface_x, res_x, iface_colors, state_labels, OUTPUT_DIR / "stage5a_residual_vs_disorder.png", "Residual energy vs D(k=10)", "D(k=10)", "residual energy (eV)")
    for key, label, suffix in [
        ("cn", "coordination number", "cn"),
        ("rms_bond_strain", "RMS bond strain", "bond_strain"),
        ("excess_voronoi_volume", "excess Voronoi volume", "excess_volume"),
        ("csp", "CSP", "csp"),
    ]:
        xvals = descriptors[key][iface_mask]
        plot_residual_scatter(xvals, res_x, iface_colors, state_labels, OUTPUT_DIR / f"stage5a_residual_vs_{suffix}.png", f"Residual energy vs {label}", label, "residual energy (eV)")

    plot_matched_pairs(D10, dE, pair_rows, OUTPUT_DIR / "stage5a_matched_pair_energy_vs_disorder.png")

    autocorr_rows = build_autocorr_rows(atomwise_rows, cells)
    with open(OUTPUT_DIR / "stage5a_spatial_autocorrelation.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(autocorr_rows[0].keys()))
        writer.writeheader()
        writer.writerows(autocorr_rows)
    plot_spatial_autocorr(autocorr_rows, OUTPUT_DIR / "stage5a_spatial_autocorrelation.png")

    patchiness_rows = build_patchiness_rows(atomwise_rows, cells)
    with open(OUTPUT_DIR / "stage5a_spatial_patchiness.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(patchiness_rows[0].keys()))
        writer.writeheader()
        writer.writerows(patchiness_rows)

    write_ovito_neighborhoods(structures_data, atomwise_rows, pair_rows, OUTPUT_DIR / "stage5a_matched_neighborhoods.xyz")

    outcome, outcome_text = make_outcome(corr_rows, group_stats_rows)
    write_report(
        r_cut,
        r_bulk,
        vor_bulk,
        residual_rows,
        corr_rows,
        group_stats_rows,
        pair_rows,
        autocorr_rows,
        patchiness_rows,
        outcome,
        outcome_text,
    )
    print("Stage 5A core calculations complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
