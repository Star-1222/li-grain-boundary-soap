"""Stage 2.75: full-population SOAP variability and metric audit.

This stage expands Stage 2.5 from 33 representative atoms to all
interface atoms plus reproducible bulk samples.  It does not define
LAEs, thresholds, disorder, clustering, or dimensionality reduction.
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

from cfg_reader import get_auxiliary, read_cfg
from stage2_profiles import (
    compute_axis_profile,
    find_interface_bins,
    identify_gb_normal,
)


PROJECT_ROOT = Path(__file__).resolve().parent
CFG_FILES = {
    "stable": PROJECT_ROOT / "ni_stable.cfg",
    "metastable": PROJECT_ROOT / "ni_metastable.cfg",
}
OUTPUT_DIR = PROJECT_ROOT / "outputs"

SOAP_PARAMS = {
    "r_cut": 3.25,
    "n_max": 12,
    "l_max": 12,
    "sigma": 0.5,
    "species": ["Ni"],
    "periodic": True,
}
BULK_SAMPLE_SIZE = 500
RANDOM_SEED = 20260817
PAIR_SAMPLE_SIZE = 100_000
MIN_GROUP_COUNT = 10
DISTANCE_TOLERANCE = 1e-10


def _interface_mask_for_structure(atoms) -> tuple[np.ndarray, str, list]:
    profiles = {
        axis: compute_axis_profile(atoms, i) for i, axis in enumerate("xyz")
    }
    normal = identify_gb_normal(profiles)
    profile = profiles[normal]
    bin_mask, runs = find_interface_bins(profile)
    atom_mask = bin_mask[profile["bin_indices"]]
    return atom_mask, normal, runs


def _summary(values: np.ndarray) -> dict:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {
            "N": 0,
            "min": float("nan"),
            "p5": float("nan"),
            "p25": float("nan"),
            "median": float("nan"),
            "mean": float("nan"),
            "p75": float("nan"),
            "p95": float("nan"),
            "max": float("nan"),
            "std": float("nan"),
        }
    return {
        "N": int(values.size),
        "min": float(values.min()),
        "p5": float(np.percentile(values, 5)),
        "p25": float(np.percentile(values, 25)),
        "median": float(np.median(values)),
        "mean": float(values.mean()),
        "p75": float(np.percentile(values, 75)),
        "p95": float(np.percentile(values, 95)),
        "max": float(values.max()),
        "std": float(values.std(ddof=1)) if values.size > 1 else 0.0,
    }


def _sample_within(rng, indices: np.ndarray, count: int) -> tuple[np.ndarray, np.ndarray]:
    m = len(indices)
    total = m * (m - 1) // 2
    count = min(count, total)
    flat = rng.choice(total, size=count, replace=False)
    i = (
        m
        - 2
        - np.floor(
            np.sqrt(-8 * flat + 4 * m * (m - 1) - 7) / 2 - 0.5
        ).astype(int)
    )
    before_i = i * m - i * (i + 1) // 2
    j = flat - before_i + i + 1
    return indices[i], indices[j]


def _sample_cross(rng, first: np.ndarray, second: np.ndarray, count: int) -> tuple[np.ndarray, np.ndarray]:
    count = min(count, len(first) * len(second))
    pairs = np.empty((0, 2), dtype=int)
    while len(pairs) < count:
        need = count - len(pairs)
        i = rng.integers(0, len(first), size=need * 2)
        j = rng.integers(0, len(second), size=need * 2)
        new = np.unique(np.stack([i, j], axis=1), axis=0)
        pairs = np.unique(np.vstack([pairs, new]), axis=0)
    pairs = pairs[:count]
    return first[pairs[:, 0]], second[pairs[:, 1]]


def _plot_norm_distributions(groups: dict, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 6))
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    for (label, values), color in zip(groups.items(), colors):
        ax.hist(values, bins=50, alpha=0.45, density=True, label=label, color=color)
    ax.set_xlabel("raw SOAP vector L2 norm")
    ax.set_ylabel("density")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _plot_distance_distributions(groups: dict, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 6))
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    for (label, values), color in zip(groups.items(), colors):
        ax.hist(values, bins=60, alpha=0.45, density=True, label=label, color=color)
    ax.set_xlabel("d_to_bulk_mean (raw Euclidean SOAP distance)")
    ax.set_ylabel("density")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _plot_decomposition(
    d_mag: np.ndarray,
    d_ang: np.ndarray,
    structures: list[str],
    regions: list[str],
    path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    markers = {"bulk": ".", "interface": "x"}
    for structure in ["stable", "metastable"]:
        for region in ["bulk", "interface"]:
            mask = np.asarray(
                [s == structure and r == region for s, r in zip(structures, regions)]
            )
            ax.scatter(
                d_mag[mask],
                d_ang[mask],
                s=6,
                alpha=0.4,
                marker=markers[region],
                label=f"{structure} {region}",
            )
    ax.set_xlabel("d_mag (magnitude contribution, raw units)")
    ax.set_ylabel("d_ang (angular contribution, raw units)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _plot_distance_vs_peatom(
    d_to_bulk: np.ndarray,
    peatoms: np.ndarray,
    structures: list[str],
    regions: list[str],
    path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    markers = {"bulk": ".", "interface": "x"}
    for structure in ["stable", "metastable"]:
        for region in ["bulk", "interface"]:
            mask = np.asarray(
                [s == structure and r == region for s, r in zip(structures, regions)]
            )
            ax.scatter(
                peatoms[mask],
                d_to_bulk[mask],
                s=6,
                alpha=0.4,
                marker=markers[region],
                label=f"{structure} {region}",
            )
    ax.set_xlabel("c_PEATOM")
    ax.set_ylabel("d_to_bulk_mean")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _plot_distance_vs_z(
    d_to_bulk: np.ndarray,
    z_positions: np.ndarray,
    structures: list[str],
    regions: list[str],
    path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    markers = {"bulk": ".", "interface": "x"}
    for structure in ["stable", "metastable"]:
        for region in ["bulk", "interface"]:
            mask = np.asarray(
                [s == structure and r == region for s, r in zip(structures, regions)]
            )
            ax.scatter(
                z_positions[mask],
                d_to_bulk[mask],
                s=6,
                alpha=0.4,
                marker=markers[region],
                label=f"{structure} {region}",
            )
    ax.set_xlabel("z position (Angstrom)")
    ax.set_ylabel("d_to_bulk_mean")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RANDOM_SEED)

    atoms_by_structure = {}
    interface_counts = {}
    records: list[dict] = []
    cat_indices = {
        "stable_bulk": [],
        "stable_interface": [],
        "metastable_bulk": [],
        "metastable_interface": [],
    }
    indices_by_structure = {}

    for structure, path in CFG_FILES.items():
        atoms = read_cfg(path)
        atoms_by_structure[structure] = atoms
        interface_mask, normal, runs = _interface_mask_for_structure(atoms)
        if normal != "z":
            raise RuntimeError(
                f"Expected z as GB normal for {structure}, got {normal}"
            )
        interface_idx = np.flatnonzero(interface_mask)
        bulk_pool = np.flatnonzero(~interface_mask)
        bulk_idx = rng.choice(bulk_pool, size=BULK_SAMPLE_SIZE, replace=False)
        interface_counts[structure] = int(interface_idx.size)
        order = np.concatenate([interface_idx, bulk_idx])
        indices_by_structure[structure] = order
        positions = atoms.positions
        scaled = atoms.get_scaled_positions()
        ids = get_auxiliary(atoms, "id").astype(int)
        ptms = get_auxiliary(atoms, "c_PTM_1").astype(int)
        cnas = get_auxiliary(atoms, "c_CNA").astype(int)
        peatoms = get_auxiliary(atoms, "c_PEATOM")

        for atom_index in order:
            region = "interface" if interface_mask[atom_index] else "bulk"
            cat_indices[f"{structure}_{region}"].append(len(records))
            records.append(
                {
                    "structure": structure,
                    "region": region,
                    "atom_id": int(ids[atom_index]),
                    "x": float(positions[atom_index, 0]),
                    "y": float(positions[atom_index, 1]),
                    "z": float(positions[atom_index, 2]),
                    "z_scaled": float(scaled[atom_index, 2]),
                    "ptm": int(ptms[atom_index]),
                    "cna": int(cnas[atom_index]),
                    "peatom": float(peatoms[atom_index]),
                }
            )

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        soap = SOAP(
            species=SOAP_PARAMS["species"],
            periodic=SOAP_PARAMS["periodic"],
            r_cut=SOAP_PARAMS["r_cut"],
            n_max=SOAP_PARAMS["n_max"],
            l_max=SOAP_PARAMS["l_max"],
            sigma=SOAP_PARAMS["sigma"],
            average="off",
            sparse=False,
            dtype="float64",
        )
        blocks = []
        for structure in CFG_FILES:
            indices = indices_by_structure[structure]
            blocks.append(soap.create(atoms_by_structure[structure], centers=indices))
        vectors = np.vstack(blocks)

    np.save(OUTPUT_DIR / "stage2_75_soap_vectors.npy", vectors)

    structures = [record["structure"] for record in records]
    regions = [record["region"] for record in records]
    ptms = np.asarray([record["ptm"] for record in records], dtype=int)
    cnas = np.asarray([record["cna"] for record in records], dtype=int)
    peatoms = np.asarray([record["peatom"] for record in records])
    z_positions = np.asarray([record["z"] for record in records])

    norms = np.linalg.norm(vectors, axis=1)
    normalized = vectors / norms[:, None]
    bulk_masks = {
        structure: np.asarray(
            [s == structure and r == "bulk" for s, r in zip(structures, regions)]
        )
        for structure in CFG_FILES
    }
    bulk_means = np.vstack(
        [vectors[bulk_masks[structure]].mean(axis=0) for structure in CFG_FILES]
    )
    np.save(OUTPUT_DIR / "stage2_75_bulk_reference_vectors.npy", bulk_means)

    bulk_norms = np.linalg.norm(bulk_means, axis=1)
    n_atoms = vectors.shape[0]
    d_to_bulk = np.zeros(n_atoms)
    cos_bulk = np.zeros(n_atoms)
    d_mag_sq = np.zeros(n_atoms)
    d_ang_sq = np.zeros(n_atoms)
    d_mag = np.zeros(n_atoms)
    d_ang = np.zeros(n_atoms)
    f_mag = np.full(n_atoms, np.nan)
    f_ang = np.full(n_atoms, np.nan)
    residual = np.zeros(n_atoms)

    for i, structure in enumerate(structures):
        structure_index = 0 if structure == "stable" else 1
        b = bulk_means[structure_index]
        norm_b = bulk_norms[structure_index]
        a = vectors[i]
        d_to_bulk[i] = np.linalg.norm(a - b)
        cos_bulk[i] = float(np.dot(a, b) / (norms[i] * norm_b))
        d_mag_sq[i] = (norms[i] - norm_b) ** 2
        d_ang_sq[i] = max(
            0.0, 2.0 * norms[i] * norm_b * (1.0 - cos_bulk[i])
        )
        d_mag[i] = np.sqrt(d_mag_sq[i])
        d_ang[i] = np.sqrt(d_ang_sq[i])
        residual[i] = (
            d_to_bulk[i] ** 2 - (d_mag_sq[i] + d_ang_sq[i])
        )
        if d_to_bulk[i] > DISTANCE_TOLERANCE:
            f_mag[i] = d_mag_sq[i] / d_to_bulk[i] ** 2
            f_ang[i] = d_ang_sq[i] / d_to_bulk[i] ** 2

    max_decomposition_residual = float(np.max(np.abs(residual)))

    for i, record in enumerate(records):
        record["soap_norm"] = float(norms[i])
        record["d_to_bulk_mean"] = float(d_to_bulk[i])
        record["cos_to_bulk_mean"] = float(cos_bulk[i])
        record["d_mag_sq"] = float(d_mag_sq[i])
        record["d_ang_sq"] = float(d_ang_sq[i])
        record["d_mag"] = float(d_mag[i])
        record["d_ang"] = float(d_ang[i])
        record["f_mag"] = float(f_mag[i])
        record["f_ang"] = float(f_ang[i])

    metadata_fields = [
        "structure",
        "region",
        "atom_id",
        "x",
        "y",
        "z",
        "z_scaled",
        "ptm",
        "cna",
        "peatom",
        "soap_norm",
        "d_to_bulk_mean",
        "cos_to_bulk_mean",
        "d_mag_sq",
        "d_ang_sq",
        "d_mag",
        "d_ang",
        "f_mag",
        "f_ang",
    ]
    with open(OUTPUT_DIR / "stage2_75_atom_metadata.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=metadata_fields)
        writer.writeheader()
        writer.writerows(records)

    category_defs = [
        ("stable bulk", "stable_bulk"),
        ("stable interface", "stable_interface"),
        ("metastable bulk", "metastable_bulk"),
        ("metastable interface", "metastable_interface"),
    ]
    metric_arrays = {
        "soap_norm": norms,
        "d_to_bulk_mean": d_to_bulk,
        "cos_to_bulk_mean": cos_bulk,
        "d_mag_sq": d_mag_sq,
        "d_ang_sq": d_ang_sq,
        "d_mag": d_mag,
        "d_ang": d_ang,
        "f_mag": f_mag,
        "f_ang": f_ang,
        "peatom": peatoms,
    }
    population_rows = []
    norm_groups = {}
    distance_groups = {}
    for label, key in category_defs:
        structure_name, region_name = label.split(" ", 1)
        mask = np.asarray(
            [
                s == structure_name and r == region_name
                for s, r in zip(structures, regions)
            ]
        )
        norm_groups[label] = norms[mask]
        distance_groups[label] = d_to_bulk[mask]
        for metric, values in metric_arrays.items():
            stats = _summary(values[mask])
            population_rows.append({"category": label, "metric": metric, **stats})

    with open(OUTPUT_DIR / "stage2_75_population_summary.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "category",
                "metric",
                "N",
                "min",
                "p5",
                "p25",
                "median",
                "mean",
                "p75",
                "p95",
                "max",
                "std",
            ],
        )
        writer.writeheader()
        writer.writerows(population_rows)

    bulk_variability = {}
    for structure in CFG_FILES:
        mask = bulk_masks[structure]
        x = vectors[mask]
        diff = x[:, None, :] - x[None, :, :]
        pair_distances = np.linalg.norm(diff, axis=2)
        upper = pair_distances[np.triu_indices(len(x), 1)]
        bulk_variability[structure] = {
            "median": float(np.median(upper)),
            "mean": float(upper.mean()),
            "std": float(upper.std(ddof=1)),
            "p5": float(np.percentile(upper, 5)),
            "p95": float(np.percentile(upper, 95)),
        }

    ptm_cna_rows = []
    for structure in CFG_FILES:
        struct_mask = np.asarray([s == structure for s in structures])
        group_specs = [
            ("ptm", ptms),
            ("cna", cnas),
            ("ptm_cna", np.asarray([f"PTM{p}-CNA{c}" for p, c in zip(ptms, cnas)])),
        ]
        for group_type, group_values in group_specs:
            for group_value in np.unique(group_values[struct_mask]):
                mask = struct_mask & (group_values == group_value)
                if int(mask.sum()) < MIN_GROUP_COUNT:
                    continue
                ptm_cna_rows.append(
                    {
                        "structure": structure,
                        "group_type": group_type,
                        "group_value": str(group_value),
                        "N": int(mask.sum()),
                        "soap_norm_median": float(np.median(norms[mask])),
                        "soap_norm_iqr": float(
                            np.percentile(norms[mask], 75)
                            - np.percentile(norms[mask], 25)
                        ),
                        "d_to_bulk_median": float(np.median(d_to_bulk[mask])),
                        "d_to_bulk_iqr": float(
                            np.percentile(d_to_bulk[mask], 75)
                            - np.percentile(d_to_bulk[mask], 25)
                        ),
                        "cos_to_bulk_median": float(np.median(cos_bulk[mask])),
                        "cos_to_bulk_iqr": float(
                            np.percentile(cos_bulk[mask], 75)
                            - np.percentile(cos_bulk[mask], 25)
                        ),
                        "peatom_median": float(np.median(peatoms[mask])),
                        "peatom_iqr": float(
                            np.percentile(peatoms[mask], 75)
                            - np.percentile(peatoms[mask], 25)
                        ),
                    }
                )

    with open(OUTPUT_DIR / "stage2_75_ptm_cna_summary.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "structure",
                "group_type",
                "group_value",
                "N",
                "soap_norm_median",
                "soap_norm_iqr",
                "d_to_bulk_median",
                "d_to_bulk_iqr",
                "cos_to_bulk_median",
                "cos_to_bulk_iqr",
                "peatom_median",
                "peatom_iqr",
            ],
        )
        writer.writeheader()
        writer.writerows(ptm_cna_rows)

    pair_categories = [
        ("stable bulk-bulk", "stable_bulk", "stable_bulk"),
        ("stable bulk-interface", "stable_bulk", "stable_interface"),
        ("stable interface-interface", "stable_interface", "stable_interface"),
        ("metastable bulk-bulk", "metastable_bulk", "metastable_bulk"),
        ("metastable bulk-interface", "metastable_bulk", "metastable_interface"),
        ("metastable interface-interface", "metastable_interface", "metastable_interface"),
        ("stable-interface vs metastable-interface", "stable_interface", "metastable_interface"),
    ]
    pairwise_summary_rows = []
    with open(OUTPUT_DIR / "stage2_75_pairwise_sample.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "category",
                "row_i",
                "row_j",
                "d_raw",
                "cos",
                "d_norm",
            ]
        )
        for category, key1, key2 in pair_categories:
            idx1 = np.asarray(cat_indices[key1], dtype=int)
            idx2 = np.asarray(cat_indices[key2], dtype=int)
            if key1 == key2:
                pair_i, pair_j = _sample_within(rng, idx1, PAIR_SAMPLE_SIZE)
            else:
                pair_i, pair_j = _sample_cross(rng, idx1, idx2, PAIR_SAMPLE_SIZE)
            vi = vectors[pair_i]
            vj = vectors[pair_j]
            d_raw = np.linalg.norm(vi - vj, axis=1)
            cos = np.einsum("ij,ij->i", vi, vj) / (
                norms[pair_i] * norms[pair_j]
            )
            d_norm = np.linalg.norm(
                normalized[pair_i] - normalized[pair_j], axis=1
            )
            for k in range(len(pair_i)):
                writer.writerow(
                    [
                        category,
                        int(pair_i[k]),
                        int(pair_j[k]),
                        f"{d_raw[k]:.9f}",
                        f"{cos[k]:.9f}",
                        f"{d_norm[k]:.9f}",
                    ]
                )
            for metric, values in [
                ("d_raw", d_raw),
                ("cos", cos),
                ("d_norm", d_norm),
            ]:
                stats = _summary(values)
                pairwise_summary_rows.append(
                    {"category": category, "metric": metric, **stats}
                )

    with open(OUTPUT_DIR / "stage2_75_pairwise_summary.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "category",
                "metric",
                "N",
                "min",
                "p5",
                "p25",
                "median",
                "mean",
                "p75",
                "p95",
                "max",
                "std",
            ],
        )
        writer.writeheader()
        writer.writerows(pairwise_summary_rows)

    _plot_norm_distributions(
        norm_groups,
        OUTPUT_DIR / "stage2_75_norm_distributions.png",
    )
    _plot_distance_distributions(
        distance_groups,
        OUTPUT_DIR / "stage2_75_bulk_distance_distributions.png",
    )
    _plot_decomposition(
        d_mag,
        d_ang,
        structures,
        regions,
        OUTPUT_DIR / "stage2_75_distance_decomposition.png",
    )
    _plot_distance_vs_peatom(
        d_to_bulk,
        peatoms,
        structures,
        regions,
        OUTPUT_DIR / "stage2_75_distance_vs_PEATOM.png",
    )
    _plot_distance_vs_z(
        d_to_bulk,
        z_positions,
        structures,
        regions,
        OUTPUT_DIR / "stage2_75_distance_vs_z.png",
    )

    _write_report(
        OUTPUT_DIR / "stage2_75_report.md",
        interface_counts,
        records,
        vectors,
        bulk_means,
        norms,
        d_to_bulk,
        cos_bulk,
        d_mag,
        d_ang,
        f_mag,
        f_ang,
        peatoms,
        max_decomposition_residual,
        population_rows,
        bulk_variability,
        ptm_cna_rows,
        pairwise_summary_rows,
        caught,
    )
    print("Stage 2.75 complete.")
    print(f"vector matrix shape: {vectors.shape}")
    print(f"interface counts: {interface_counts}")
    print(f"max decomposition residual: {max_decomposition_residual:.3e}")
    return 0


def _write_report(
    path,
    interface_counts,
    records,
    vectors,
    bulk_means,
    norms,
    d_to_bulk,
    cos_bulk,
    d_mag,
    d_ang,
    f_mag,
    f_ang,
    peatoms,
    max_residual,
    population_rows,
    bulk_variability,
    ptm_cna_rows,
    pairwise_summary_rows,
    caught,
) -> None:
    structures = [record["structure"] for record in records]
    regions = [record["region"] for record in records]
    stable_bulk = np.asarray(
        [s == "stable" and r == "bulk" for s, r in zip(structures, regions)]
    )
    stable_interface = np.asarray(
        [s == "stable" and r == "interface" for s, r in zip(structures, regions)]
    )
    meta_bulk = np.asarray(
        [s == "metastable" and r == "bulk" for s, r in zip(structures, regions)]
    )
    meta_interface = np.asarray(
        [s == "metastable" and r == "interface" for s, r in zip(structures, regions)]
    )

    lines = [
        "# Stage 2.75: full-population SOAP variability audit",
        "",
        "## A. Population definition",
        "",
        f"- Interface masks were reused from Stage 2 profiling (200 bins, "
        f"non-dominant fraction > 0.1, min 50 atoms per bin), not "
        f"redefined.",
        f"- Stable interface atoms: {interface_counts['stable']}.",
        f"- Metastable interface atoms: {interface_counts['metastable']}.",
        f"- Bulk reference samples: {sum(int(mask.sum()) for mask in [stable_bulk, meta_bulk])} "
        f"atoms total, 500 per structure, drawn with `np.random.default_rng("
        f"{RANDOM_SEED})`.",
        "",
        "## B. SOAP calculation",
        "",
        f"- Raw DScribe SOAP vectors: shape `{vectors.shape}`, feature "
        f"dimension {vectors.shape[1]}.",
        "- Parameters: `r_cut=3.25`, `n_max=12`, `l_max=12`, "
        "`sigma=0.5`, species `[\"Ni\"]`, periodic.",
        "- Vectors were saved without L2 normalization.",
        "",
        "## C. Bulk reference variability",
        "",
        "Pairwise raw Euclidean SOAP distances among the 500 bulk atoms:",
        "",
        "| structure | median | mean | std | p5 | p95 |",
        "|---|---|---|---|---|---|",
    ]
    for structure in ["stable", "metastable"]:
        row = bulk_variability[structure]
        lines.append(
            f"| {structure} | {row['median']:.6f} | {row['mean']:.6f} | "
            f"{row['std']:.6f} | {row['p5']:.6f} | {row['p95']:.6f} |"
        )
    lines.extend(
        [
            "",
            "Stable-bulk and metastable-bulk distributions are compared "
            "descriptively in the population summary; no identity "
            "assumption is made.",
            "",
            "## D. SOAP norm distributions",
            "",
            "| category | N | min | p5 | p25 | median | mean | p75 | p95 | max | std |",
            "|---|---|---|---|---|---|---|---|---|---|---|",
        ]
    )
    for row in population_rows:
        if row["metric"] == "soap_norm":
            lines.append(
                f"| {row['category']} | {row['N']} | {row['min']:.6f} | "
                f"{row['p5']:.6f} | {row['p25']:.6f} | {row['median']:.6f} | "
                f"{row['mean']:.6f} | {row['p75']:.6f} | {row['p95']:.6f} | "
                f"{row['max']:.6f} | {row['std']:.6f} |"
            )
    lines.extend(
        [
            "",
            "Lower SOAP norm is not interpreted as greater disorder.",
            "",
            "## E. Distance from bulk reference",
            "",
            "| category | N | min | p5 | p25 | median | mean | p75 | p95 | max | std |",
            "|---|---|---|---|---|---|---|---|---|---|---|",
        ]
    )
    for row in population_rows:
        if row["metric"] in ("d_to_bulk_mean", "cos_to_bulk_mean"):
            lines.append(
                f"| {row['category']} {row['metric']} | {row['N']} | "
                f"{row['min']:.6f} | {row['p5']:.6f} | {row['p25']:.6f} | "
                f"{row['median']:.6f} | {row['mean']:.6f} | {row['p75']:.6f} | "
                f"{row['p95']:.6f} | {row['max']:.6f} | {row['std']:.6f} |"
            )
    lines.extend(
        [
            "",
            "## F. Magnitude/angular decomposition",
            "",
            f"- Maximum absolute residual of "
            f"`d_raw^2 = d_mag^2 + d_ang^2`: `{max_residual:.3e}`.",
            "- `f_mag` and `f_ang` are reported for atoms with "
            f"`d_raw > {DISTANCE_TOLERANCE}`.",
            "- These fractions are mathematical decompositions only and "
            "are not labeled as disorder.",
            "",
            "| category | metric | N | median | mean |",
            "|---|---|---|---|---|",
        ]
    )
    for row in population_rows:
        if row["metric"] in (
            "d_mag_sq",
            "d_ang_sq",
            "d_mag",
            "d_ang",
            "f_mag",
            "f_ang",
        ):
            lines.append(
                f"| {row['category']} {row['metric']} | {row['N']} | "
                f"{row['median']:.6f} | {row['mean']:.6f} |"
            )
    lines.extend(
        [
            "",
            "## G. Relationship with PTM/CNA",
            "",
            "Groups with at least "
            f"{MIN_GROUP_COUNT} atoms are summarized in "
            "`stage2_75_ptm_cna_summary.csv`. The data are described in "
            "terms of narrow, broad, or overlapping SOAP distributions; "
            "SOAP is not said to validate PTM/CNA.",
            "",
            "| structure | group | N | norm median | norm IQR | d median | d IQR | cos median | cos IQR | PEATOM median | PEATOM IQR |",
            "|---|---|---|---|---|---|---|---|---|---|---|",
        ]
    )
    for row in ptm_cna_rows:
        group_display = row["group_value"]
        if row["group_type"] == "ptm":
            group_display = f"PTM {row['group_value']}"
        elif row["group_type"] == "cna":
            group_display = f"CNA {row['group_value']}"
        lines.append(
            f"| {row['structure']} | {group_display} | {row['N']} | "
            f"{row['soap_norm_median']:.4f} | {row['soap_norm_iqr']:.4f} | "
            f"{row['d_to_bulk_median']:.4f} | {row['d_to_bulk_iqr']:.4f} | "
            f"{row['cos_to_bulk_median']:.4f} | {row['cos_to_bulk_iqr']:.4f} | "
            f"{row['peatom_median']:.4f} | {row['peatom_iqr']:.4f} |"
        )
    lines.extend(
        [
            "",
            "## H. Relationship with PEATOM",
            "",
            "- Population summaries include PEATOM statistics per "
            "bulk/interface category.",
            "- `d_to_bulk_mean` vs `c_PEATOM` is plotted in "
            "`stage2_75_distance_vs_PEATOM.png`; the relationship is "
            "descriptive.",
            "",
            "## I. Stable vs metastable descriptive comparison",
            "",
            "The population summary and pairwise sample allow descriptive "
            "comparison of SOAP norm, distance-to-bulk-reference, cosine, "
            "decomposition fractions, and PEATOM. No claim of more or "
            "less disorder is made.",
            "",
            "## J. Numerical checks and limitations",
            "",
            f"- Raw vectors are finite: {bool(np.isfinite(vectors).all())}.",
            f"- Decomposition residual max: `{max_residual:.3e}`.",
            f"- Random seed used for bulk samples and pair sampling: "
            f"`{RANDOM_SEED}`.",
            f"- Pair sample size per category: up to {PAIR_SAMPLE_SIZE}.",
        ]
    )
    if caught:
        messages = sorted({str(item.message) for item in caught})
        lines.append(
            f"- DScribe LinAlgWarning occurred during construction "
            f"({len(messages)} message(s)); vectors remained finite. "
            "The warning was already characterized in Stage 2.5."
        )
        lines.extend(f"  - `{message}`" for message in messages)
    lines.extend(
        [
            "",
            "## K. What Stage 2.75 supports",
            "",
            "- Raw SOAP vectors for the full interface populations and "
            "reproducible bulk samples.",
            "- Descriptive distributions for SOAP norms, distances to "
            "bulk reference vectors, cosine similarities, and "
            "magnitude/angular decomposition fractions.",
            "- Empirical bulk SOAP variability baselines for both "
            "structures.",
            "- Pairwise population sampling with documented random seed.",
            "",
            "## L. What Stage 2.75 does NOT support",
            "",
            "- It does not define LAEs, thresholds, LER, PCA/UMAP/t-SNE, "
            "clusters, or disorder.",
            "- It does not claim that stable or metastable is more "
            "ordered/disordered.",
            "- It does not claim that discrete PTM/CNA classes are "
            "validated by SOAP.",
            "",
            "## M. Recommendation for Stage 3",
            "",
            "Use the bulk-variability baseline and distance decomposition "
            "to design a controlled, DScribe-specific LAE equivalence "
            "criterion before constructing LAEs/LER. Do not reuse the "
            "historical `epsilon = 0.0025` on DScribe vectors without "
            "recalibration.",
        ]
    )
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
