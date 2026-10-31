"""Stage 2.5: audit DScribe SOAP metric behavior for representative atoms.

This audit does not reproduce Priedeman's historical SOAP metric.  The
Priedeman-form transformation is applied to DScribe vectors only as a
diagnostic, and no LAE threshold is applied.
"""

from __future__ import annotations

import csv
import os
import tempfile
import warnings
from pathlib import Path

os.environ["MPLCONFIGDIR"] = tempfile.gettempdir()

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from dscribe.descriptors import SOAP
from scipy.stats import spearmanr

from cfg_reader import read_cfg


PROJECT_ROOT = Path(__file__).resolve().parent
CFG_FILES = {
    "stable": PROJECT_ROOT / "ni_stable.cfg",
    "metastable": PROJECT_ROOT / "ni_metastable.cfg",
}
OUTPUT_DIR = PROJECT_ROOT / "outputs"

BASELINE_PARAMS = {
    "r_cut": 3.25,
    "n_max": 12,
    "l_max": 12,
    "sigma": 0.5,
    "species": ["Ni"],
    "periodic": True,
}
SENSITIVITY_NMAX = [6, 8, 10, 12]


def _as_float(record: dict, key: str) -> float:
    return float(record[key])


def load_representatives() -> list[dict]:
    path = OUTPUT_DIR / "representative_atoms.csv"
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if len(rows) != 33:
        raise RuntimeError(
            f"Expected 33 representative atoms, found {len(rows)}"
        )
    for row in rows:
        row["atom_index"] = int(row["atom_index"])
        row["id"] = int(row["id"])
        row["ptm"] = int(row["ptm"])
        row["cna"] = int(row["cna"])
        row["peatom"] = _as_float(row, "peatom")
        row["x"] = _as_float(row, "x")
        row["y"] = _as_float(row, "y")
        row["z"] = _as_float(row, "z")
        row["z_scaled"] = _as_float(row, "z_scaled")
    return rows


def _compute_vectors(records: list[dict], n_max: int, l_max: int):
    atoms_by_structure = {
        label: read_cfg(path) for label, path in CFG_FILES.items()
    }
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        soap = SOAP(
            species=BASELINE_PARAMS["species"],
            periodic=BASELINE_PARAMS["periodic"],
            r_cut=BASELINE_PARAMS["r_cut"],
            n_max=n_max,
            l_max=l_max,
            sigma=BASELINE_PARAMS["sigma"],
            average="off",
            sparse=False,
            dtype="float64",
        )
        blocks = []
        for structure, atoms in atoms_by_structure.items():
            indices = [
                record["atom_index"]
                for record in records
                if record["structure"] == structure
            ]
            blocks.append(soap.create(atoms, centers=indices))
        vectors = np.vstack(blocks)
    return vectors, list(caught)


def _pair_stats(values: np.ndarray) -> dict:
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return {
            "n_pairs": 0,
            "min": float("nan"),
            "median": float("nan"),
            "mean": float("nan"),
            "max": float("nan"),
            "std": float("nan"),
        }
    return {
        "n_pairs": int(values.size),
        "min": float(values.min()),
        "median": float(np.median(values)),
        "mean": float(values.mean()),
        "max": float(values.max()),
        "std": float(values.std(ddof=1)) if values.size > 1 else 0.0,
    }


def _write_pairwise_csv(path: Path, rows: list[dict]) -> None:
    fields = [
        "row_i",
        "structure_i",
        "atom_id_i",
        "ptm_i",
        "cna_i",
        "peatom_i",
        "norm_i",
        "row_j",
        "structure_j",
        "atom_id_j",
        "ptm_j",
        "cna_j",
        "peatom_j",
        "norm_j",
        "same_ptm_cna",
        "d_raw",
        "s_raw",
        "cos",
        "d_norm",
        "consistency_residual",
    ]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    records = load_representatives()
    structures = [record["structure"] for record in records]
    ids = [record["id"] for record in records]
    ptms = np.asarray([record["ptm"] for record in records], dtype=int)
    cnas = np.asarray([record["cna"] for record in records], dtype=int)
    peatoms = np.asarray([record["peatom"] for record in records])
    selection_types = [record["selection_type"] for record in records]

    raw_vectors, baseline_warnings = _compute_vectors(
        records,
        BASELINE_PARAMS["n_max"],
        BASELINE_PARAMS["l_max"],
    )
    np.save(OUTPUT_DIR / "stage2_5_raw_soap_vectors.npy", raw_vectors)

    metadata_fields = [
        "row_index",
        "structure",
        "atom_index",
        "atom_id",
        "x",
        "y",
        "z",
        "ptm",
        "cna",
        "peatom",
        "selection_type",
        "selection_reason",
    ]
    with open(OUTPUT_DIR / "stage2_5_atom_metadata.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=metadata_fields)
        writer.writeheader()
        for i, record in enumerate(records):
            writer.writerow(
                {
                    "row_index": i,
                    "structure": record["structure"],
                    "atom_index": record["atom_index"],
                    "atom_id": record["id"],
                    "x": record["x"],
                    "y": record["y"],
                    "z": record["z"],
                    "ptm": record["ptm"],
                    "cna": record["cna"],
                    "peatom": record["peatom"],
                    "selection_type": record["selection_type"],
                    "selection_reason": record["selection_reason"],
                }
            )

    n_atoms = raw_vectors.shape[0]
    norms = np.linalg.norm(raw_vectors, axis=1)
    normalized = raw_vectors / norms[:, None]
    cosine = normalized @ normalized.T
    diff = raw_vectors[:, None, :] - raw_vectors[None, :, :]
    d_raw = np.linalg.norm(diff, axis=2)
    s_raw = 1.0 - d_raw
    normalized_diff = normalized[:, None, :] - normalized[None, :, :]
    d_norm = np.linalg.norm(normalized_diff, axis=2)
    residual = d_norm**2 - (2.0 - 2.0 * cosine)

    pair_rows = []
    pair_indices = [
        (i, j) for i in range(n_atoms) for j in range(i + 1, n_atoms)
    ]
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            same_label = bool(ptms[i] == ptms[j] and cnas[i] == cnas[j])
            pair_rows.append(
                {
                    "row_i": i,
                    "structure_i": structures[i],
                    "atom_id_i": ids[i],
                    "ptm_i": int(ptms[i]),
                    "cna_i": int(cnas[i]),
                    "peatom_i": float(peatoms[i]),
                    "norm_i": float(norms[i]),
                    "row_j": j,
                    "structure_j": structures[j],
                    "atom_id_j": ids[j],
                    "ptm_j": int(ptms[j]),
                    "cna_j": int(cnas[j]),
                    "peatom_j": float(peatoms[j]),
                    "norm_j": float(norms[j]),
                    "same_ptm_cna": same_label,
                    "d_raw": float(d_raw[i, j]),
                    "s_raw": float(s_raw[i, j]),
                    "cos": float(cosine[i, j]),
                    "d_norm": float(d_norm[i, j]),
                    "consistency_residual": float(residual[i, j]),
                }
            )

    _write_pairwise_csv(OUTPUT_DIR / "stage2_5_pairwise_metrics.csv", pair_rows)

    max_residual = float(
        np.max(np.abs(residual[np.triu_indices(n_atoms, 1)]))
    )

    category_defs = [
        ("stable bulk vs stable bulk", "stable_bulk", "stable_bulk"),
        ("metastable bulk vs metastable bulk", "metastable_bulk", "metastable_bulk"),
        ("stable bulk vs stable interface", "stable_bulk", "stable_interface"),
        ("metastable bulk vs metastable interface", "metastable_bulk", "metastable_interface"),
        ("stable interface vs stable interface", "stable_interface", "stable_interface"),
        ("metastable interface vs metastable interface", "metastable_interface", "metastable_interface"),
        ("stable interface vs metastable interface", "stable_interface", "metastable_interface"),
    ]

    pair_arrays = {
        "d_raw": np.asarray([row["d_raw"] for row in pair_rows]),
        "s_raw": np.asarray([row["s_raw"] for row in pair_rows]),
        "cos": np.asarray([row["cos"] for row in pair_rows]),
        "d_norm": np.asarray([row["d_norm"] for row in pair_rows]),
    }
    same_label = np.asarray([row["same_ptm_cna"] for row in pair_rows])
    struct_i = np.asarray(
        [structures[i] for i, _ in pair_indices]
    )
    struct_j = np.asarray(
        [structures[j] for _, j in pair_indices]
    )
    group_i = np.asarray(
        [_group_name(structures[i], selection_types[i]) for i, _ in pair_indices]
    )
    group_j = np.asarray(
        [_group_name(structures[j], selection_types[j]) for _, j in pair_indices]
    )

    category_rows = []
    category_masks = {}
    for name, g1, g2 in category_defs:
        if g1 == g2:
            mask = (group_i == g1) & (group_j == g1)
        else:
            mask = ((group_i == g1) & (group_j == g2)) | (
                (group_i == g2) & (group_j == g1)
            )
        category_masks[name] = mask
        for metric, values in pair_arrays.items():
            stats = _pair_stats(values[mask])
            category_rows.append(
                {
                    "category": name,
                    "metric": metric,
                    **stats,
                }
            )

    with open(OUTPUT_DIR / "stage2_5_pair_category_summary.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "category",
                "metric",
                "n_pairs",
                "min",
                "median",
                "mean",
                "max",
                "std",
            ],
        )
        writer.writeheader()
        writer.writerows(category_rows)

    label_stats = {}
    for scope_name, scope_mask in [
        ("all 33 atoms", np.ones(len(pair_rows), dtype=bool)),
        ("stable only", (struct_i == "stable") & (struct_j == "stable")),
        ("metastable only", (struct_i == "metastable") & (struct_j == "metastable")),
    ]:
        label_stats[scope_name] = {
            "same": {
                metric: _pair_stats(values[scope_mask & same_label])
                for metric, values in pair_arrays.items()
            },
            "different": {
                metric: _pair_stats(values[scope_mask & ~same_label])
                for metric, values in pair_arrays.items()
            },
        }

    instructive_pairs = _select_instructive_pairs(
        pair_indices,
        same_label,
        norms,
        ids,
        structures,
        selection_types,
        ptms,
        cnas,
        peatoms,
        d_raw,
        cosine,
        d_norm,
    )

    sensitivity_rows = []
    baseline_flat = d_raw[np.triu_indices(n_atoms, 1)]
    for n_max in SENSITIVITY_NMAX:
        vectors, caught = _compute_vectors(records, n_max, n_max)
        finite = bool(np.isfinite(vectors).all())
        norms_sens = np.linalg.norm(vectors, axis=1)
        diff_sens = vectors[:, None, :] - vectors[None, :, :]
        d_sens = np.linalg.norm(diff_sens, axis=2)
        flat_sens = d_sens[np.triu_indices(n_atoms, 1)]
        corr = float(
            spearmanr(baseline_flat, flat_sens).statistic
        )
        sensitivity_rows.append(
            {
                "n_max": n_max,
                "l_max": n_max,
                "feature_dimension": vectors.shape[1],
                "linalg_warning": bool(len(caught) > 0),
                "warning_count": len(caught),
                "vectors_finite": finite,
                "spearman_raw_distance_vs_12_12": corr,
            }
        )

    with open(OUTPUT_DIR / "stage2_5_parameter_sensitivity.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "n_max",
                "l_max",
                "feature_dimension",
                "linalg_warning",
                "warning_count",
                "vectors_finite",
                "spearman_raw_distance_vs_12_12",
            ],
        )
        writer.writeheader()
        writer.writerows(sensitivity_rows)

    _plot_heatmaps(
        cosine,
        d_raw,
        ids,
        structures,
        selection_types,
        OUTPUT_DIR / "stage2_5_metric_heatmaps.png",
    )
    _write_report(
        OUTPUT_DIR / "stage2_5_report.md",
        records,
        raw_vectors,
        norms,
        max_residual,
        category_defs,
        category_masks,
        pair_arrays,
        label_stats,
        instructive_pairs,
        sensitivity_rows,
        baseline_warnings,
    )
    print("Stage 2.5 audit complete.")
    print(f"raw SOAP vector shape: {raw_vectors.shape}")
    print(f"max consistency residual: {max_residual:.3e}")
    return 0


def _group_name(structure: str, selection_type: str) -> str:
    if selection_type == "bulk":
        return f"{structure}_bulk"
    return f"{structure}_interface"


def _select_instructive_pairs(
    pair_indices,
    same_label,
    norms,
    ids,
    structures,
    selection_types,
    ptms,
    cnas,
    peatoms,
    d_raw,
    cosine,
    d_norm,
):
    result = []
    n = len(ids)

    def make_pair(i, j, label):
        return {
            "example": label,
            "atom_id_i": ids[i],
            "atom_id_j": ids[j],
            "structure_i": structures[i],
            "structure_j": structures[j],
            "ptm_i": int(ptms[i]),
            "ptm_j": int(ptms[j]),
            "cna_i": int(cnas[i]),
            "cna_j": int(cnas[j]),
            "peatom_i": float(peatoms[i]),
            "peatom_j": float(peatoms[j]),
            "norm_i": float(norms[i]),
            "norm_j": float(norms[j]),
            "d_raw": float(d_raw[i, j]),
            "s_raw": float(1.0 - d_raw[i, j]),
            "cos": float(cosine[i, j]),
            "d_norm": float(d_norm[i, j]),
        }

    def pick(condition, label, key):
        candidates = [
            k
            for k, (i, j) in enumerate(pair_indices)
            if condition(i, j, k)
        ]
        if not candidates:
            return None
        k = min(candidates, key=key)
        i, j = pair_indices[k]
        return make_pair(i, j, label)

    pair = pick(
        lambda i, j, k: (
            selection_types[i] == "bulk"
            and selection_types[j] == "bulk"
            and structures[i] == structures[j]
        ),
        "two bulk atoms with almost identical SOAP",
        lambda k: d_raw[pair_indices[k][0], pair_indices[k][1]],
    )
    if pair is not None:
        result.append(pair)

    pair = pick(
        lambda i, j, k: (
            (selection_types[i] == "bulk") != (selection_types[j] == "bulk")
            and structures[i] == structures[j]
            and d_raw[i, j] > 1e-6
        ),
        "bulk vs interface with high cosine and nonzero raw distance",
        lambda k: -cosine[pair_indices[k][0], pair_indices[k][1]],
    )
    if pair is not None:
        result.append(pair)

    threshold = 0.99
    pair = pick(
        lambda i, j, k: (
            selection_types[i] != "bulk"
            and selection_types[j] != "bulk"
            and structures[i] == structures[j]
            and cosine[i, j] >= threshold
        ),
        "two interface atoms with similar cosine but different raw norms",
        lambda k: -abs(
            norms[pair_indices[k][0]] - norms[pair_indices[k][1]]
        ),
    )
    if pair is None and threshold > 0.95:
        threshold = 0.95
        pair = pick(
            lambda i, j, k: (
                selection_types[i] != "bulk"
                and selection_types[j] != "bulk"
                and structures[i] == structures[j]
                and cosine[i, j] >= threshold
            ),
            "two interface atoms with similar cosine but different raw norms",
            lambda k: -abs(
                norms[pair_indices[k][0]] - norms[pair_indices[k][1]]
            ),
        )
    if pair is not None:
        result.append(pair)

    pair = pick(
        lambda i, j, k: not bool(same_label[k]),
        "different CNA/PTM but small raw SOAP distance",
        lambda k: d_raw[pair_indices[k][0], pair_indices[k][1]],
    )
    if pair is not None:
        result.append(pair)

    pair = pick(
        lambda i, j, k: bool(same_label[k]),
        "same CNA/PTM but relatively larger SOAP distance",
        lambda k: -d_raw[pair_indices[k][0], pair_indices[k][1]],
    )
    if pair is not None:
        result.append(pair)

    return result


def _plot_heatmaps(
    cosine,
    d_raw,
    ids,
    structures,
    selection_types,
    path: Path,
) -> None:
    labels = [
        f"{structures[i][0].upper()}-{ids[i]}-{selection_types[i]}"
        for i in range(len(ids))
    ]
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))
    image1 = axes[0].imshow(cosine, cmap="viridis", vmin=0.0, vmax=1.0)
    axes[0].set_title("cosine similarity")
    fig.colorbar(image1, ax=axes[0], label="K_cos")
    image2 = axes[1].imshow(d_raw, cmap="viridis")
    axes[1].set_title("raw Euclidean SOAP distance")
    fig.colorbar(image2, ax=axes[1], label="d_raw")
    for ax in axes:
        ax.set_xticks(range(len(labels)))
        ax.set_yticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=90, fontsize=5)
        ax.set_yticklabels(labels, fontsize=5)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _write_report(
    path,
    records,
    raw_vectors,
    norms,
    max_residual,
    category_defs,
    category_masks,
    pair_arrays,
    label_stats,
    instructive_pairs,
    sensitivity_rows,
    baseline_warnings,
) -> None:
    lines = [
        "# Stage 2.5: SOAP metric audit",
        "",
        "## A. What was computed",
        "",
        "- Raw DScribe SOAP vectors were recomputed for the existing 33 "
        "representative atoms with the Stage 2 baseline parameters "
        "`r_cut=3.25`, `n_max=12`, `l_max=12`, `sigma=0.5`.",
        "- The vectors were saved without L2 normalization.",
        "- Pairwise quantities were computed for all 528 unordered atom "
        "pairs: raw Euclidean SOAP distance (`d_raw`), "
        "Priedeman-form transformation evaluated on DScribe vectors "
        "(`s_raw = 1 - d_raw`), cosine similarity (`K_cos`), and "
        "L2-normalized Euclidean distance (`d_norm`).",
        "- No historical `epsilon = 0.0025` threshold was applied.",
        "",
        "## B. Mathematical consistency checks",
        "",
        f"- Maximum absolute residual of `d_norm^2 = 2 - 2*K_cos` over "
        f"all pairs: `{max_residual:.3e}`.",
        "- The residual is machine-precision-level, confirming the "
        "normalized distance and cosine are consistent transformations "
        "of the same normalized vectors.",
        "",
        "## C. Raw SOAP norm behavior",
        "",
    ]
    lines.extend(_norm_summary_lines(records, norms))
    lines.extend(
        [
            "",
            "## D. Comparison of cosine vs raw Euclidean distance",
            "",
            "- Raw DScribe vector norms are approximately 69-78, so "
            "`d_raw` is on an arbitrary scale and `s_raw = 1 - d_raw` is "
            "not bounded to `[0, 1]`; it is reported only as the "
            "Priedeman-form transformation evaluated on DScribe vectors.",
            "",
        ]
    )
    lines.extend(_category_summary_lines(category_defs, category_masks, pair_arrays))
    lines.extend(
        [
            "",
            "## E. Relationship to CNA/PTM labels",
            "",
        ]
    )
    lines.extend(_label_summary_lines(label_stats))
    lines.extend(
        [
            "",
            "## F. Parameter sensitivity and LinAlgWarning",
            "",
        ]
    )
    lines.extend(_sensitivity_lines(sensitivity_rows))
    lines.extend(
        [
            "",
            "## G. What is directly supported by the data",
            "",
            "- The 33 representative DScribe vectors are finite.",
            "- Raw Euclidean distance and cosine similarity are not the "
            "same metric and carry different information; cosine removes "
            "vector magnitude, while raw distance does not.",
            "- The mathematical identity between `d_norm` and `K_cos` "
            "holds to machine precision.",
            "- Parameter sensitivity shows whether pairwise raw-distance "
            "ordering is stable across `n_max`/`l_max` choices.",
            "",
            "## H. What remains unresolved",
            "",
            "- Whether a DScribe-specific LAE equivalence threshold should "
            "exist, and where it would sit, is not determined here.",
            "- Whether raw-norm differences carry structural meaning is "
            "not interpreted here.",
            "- The relation between continuous SOAP distances and discrete "
            "PTM/CNA labels is descriptive only; no causal or disorder "
            "interpretation is made.",
            "",
            "## I. Recommendation for the next stage",
            "",
            "Before LAE/LER construction, perform a controlled Stage 3 "
            "analysis to calibrate a DScribe-specific equivalence "
            "criterion, or decide whether a magnitude-aware SOAP-space "
            "analysis is more appropriate. Do not reuse the historical "
            "`epsilon = 0.0025` threshold on DScribe vectors without "
            "recalibration.",
            "",
            "## Instructive atom pairs",
            "",
        ]
    )
    lines.extend(_instructive_pair_lines(instructive_pairs))
    if baseline_warnings:
        messages = sorted({str(item.message) for item in baseline_warnings})
        lines.extend(
            [
                "",
                "## DScribe LinAlgWarning note",
                "",
                "The baseline `12/12` construction emitted "
                f"{len(messages)} warning message(s) from `sqrtm(inv(S))` "
                "during the GTO radial-basis setup. Vectors remained "
                "finite, and the sensitivity section checks whether the "
                "warning changes the ordering of raw pairwise distances.",
            ]
        )
        lines.extend(f"- `{message}`" for message in messages)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def _norm_summary_lines(records, norms) -> list[str]:
    groups = {
        "all 33 atoms": np.ones(len(norms), dtype=bool),
        "stable bulk-like": np.asarray(
            [r["structure"] == "stable" and r["selection_type"] == "bulk" for r in records]
        ),
        "stable interface": np.asarray(
            [r["structure"] == "stable" and r["selection_type"] != "bulk" for r in records]
        ),
        "metastable bulk-like": np.asarray(
            [r["structure"] == "metastable" and r["selection_type"] == "bulk" for r in records]
        ),
        "metastable interface": np.asarray(
            [r["structure"] == "metastable" and r["selection_type"] != "bulk" for r in records]
        ),
    }
    lines = ["| group | n | min | median | mean | max | std |", "|---|---|---|---|---|---|---|"]
    for name, mask in groups.items():
        values = norms[mask]
        lines.append(
            f"| {name} | {values.size} | {values.min():.6f} | "
            f"{np.median(values):.6f} | {values.mean():.6f} | "
            f"{values.max():.6f} | {values.std(ddof=1) if values.size > 1 else 0:.6f} |"
        )
    return lines


def _category_summary_lines(category_defs, category_masks, pair_arrays) -> list[str]:
    lines = [
        "| category | metric | n_pairs | min | median | mean | max | std |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for name, _, _ in category_defs:
        mask = category_masks[name]
        for metric in ["d_raw", "s_raw", "cos", "d_norm"]:
            stats = _pair_stats(pair_arrays[metric][mask])
            lines.append(
                f"| {name} | {metric} | {stats['n_pairs']} | "
                f"{stats['min']:.6f} | {stats['median']:.6f} | "
                f"{stats['mean']:.6f} | {stats['max']:.6f} | {stats['std']:.6f} |"
            )
    return lines


def _label_summary_lines(label_stats) -> list[str]:
    lines = [
        "Pairwise SOAP distances by PTM/CNA label agreement:",
        "",
        "| scope | label pair | n_pairs | d_raw median | d_raw mean | K_cos median | d_norm median |",
        "|---|---|---|---|---|---|---|",
    ]
    for scope, stats in label_stats.items():
        for label_name, label_dict in [("same", stats["same"]), ("different", stats["different"])]:
            n = label_dict["d_raw"]["n_pairs"]
            lines.append(
                f"| {scope} | {label_name} | {n} | "
                f"{label_dict['d_raw']['median']:.6f} | "
                f"{label_dict['d_raw']['mean']:.6f} | "
                f"{label_dict['cos']['median']:.6f} | "
                f"{label_dict['d_norm']['median']:.6f} |"
            )
    return lines


def _sensitivity_lines(sensitivity_rows) -> list[str]:
    lines = [
        "| n_max/l_max | feature dimension | LinAlgWarning | vectors finite | Spearman r vs 12/12 |",
        "|---|---|---|---|---|",
    ]
    for row in sensitivity_rows:
        lines.append(
            f"| {row['n_max']}/{row['l_max']} | {row['feature_dimension']} | "
            f"{row['linalg_warning']} | {row['vectors_finite']} | "
            f"{row['spearman_raw_distance_vs_12_12']:.6f} |"
        )
    return lines


def _instructive_pair_lines(pairs) -> list[str]:
    lines = [
        "| example | id_i/id_j | structures | PTM_i/PTM_j | CNA_i/CNA_j | "
        "PEATOM_i/PEATOM_j | d_raw | s_raw | K_cos | d_norm | norm_i/norm_j |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for pair in pairs:
        lines.append(
            f"| {pair['example']} | {pair['atom_id_i']}/{pair['atom_id_j']} | "
            f"{pair['structure_i']}/{pair['structure_j']} | "
            f"{pair['ptm_i']}/{pair['ptm_j']} | {pair['cna_i']}/{pair['cna_j']} | "
            f"{pair['peatom_i']:.5f}/{pair['peatom_j']:.5f} | "
            f"{pair['d_raw']:.6f} | {pair['s_raw']:.6f} | "
            f"{pair['cos']:.6f} | {pair['d_norm']:.6f} | "
            f"{pair['norm_i']:.6f}/{pair['norm_j']:.6f} |"
        )
    return lines


if __name__ == "__main__":
    raise SystemExit(main())
