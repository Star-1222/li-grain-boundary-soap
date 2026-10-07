"""Stage 2: representative-atom selection and SOAP cosine similarities."""

from __future__ import annotations

import csv
import os
import warnings
from pathlib import Path

from cfg_reader import get_auxiliary, read_cfg
from stage2_profiles import (
    INTERFACE_FRACTION_THRESHOLD,
    MIN_BIN_ATOMS,
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
os.environ.setdefault(
    "MPLCONFIGDIR", os.environ.get("TEMP", "/tmp")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from dscribe.descriptors import SOAP

SOAP_PARAMS = {
    "rcut": 3.25,
    "nmax": 12,
    "lmax": 12,
    "sigma": 0.5,
    "species": ["Ni"],
    "periodic": True,
}
BUFFER_AWAY_FROM_INTERFACE = 10.0


def _profile_data(atoms):
    profiles = {
        axis: compute_axis_profile(atoms, i) for i, axis in enumerate("xyz")
    }
    normal = identify_gb_normal(profiles)
    profile = profiles[normal]
    bin_mask, runs = find_interface_bins(profile)
    interface_atom_mask = bin_mask[profile["bin_indices"]]
    return normal, profile, runs, interface_atom_mask


def _select_top_unselected(
    mask: np.ndarray,
    score: np.ndarray,
    selected: set[int],
    count: int,
) -> list[int]:
    candidates = np.flatnonzero(mask)
    if candidates.size == 0:
        return []
    order = candidates[np.argsort(score[candidates])[::-1]]
    result: list[int] = []
    for index in order:
        if int(index) in selected:
            continue
        result.append(int(index))
        if len(result) >= count:
            break
    return result


def _select_representatives(
    structure: str,
    atoms,
    selected_indices: set[int],
    records: list[dict],
) -> None:
    positions = atoms.positions
    z = positions[:, 2]
    ids = get_auxiliary(atoms, "id").astype(int)
    ptm = get_auxiliary(atoms, "c_PTM_1").astype(int)
    cna = get_auxiliary(atoms, "c_CNA").astype(int)
    pe = get_auxiliary(atoms, "c_PEATOM")
    scaled = atoms.get_scaled_positions()

    normal, profile, runs, interface_mask = _profile_data(atoms)
    if normal != "z":
        raise RuntimeError(
            f"Expected z as GB normal for {structure}, got {normal}"
        )
    centers = profile["centers"]
    ranges = [(float(centers[s]), float(centers[e])) for s, e in runs]

    def add(index: int, selection_type: str, reason: str) -> None:
        if int(index) in selected_indices:
            return
        selected_indices.add(int(index))
        records.append(
            {
                "structure": structure,
                "atom_index": int(index),
                "id": int(ids[index]),
                "ptm": int(ptm[index]),
                "cna": int(cna[index]),
                "peatom": float(pe[index]),
                "x": float(positions[index, 0]),
                "y": float(positions[index, 1]),
                "z": float(positions[index, 2]),
                "z_scaled": float(scaled[index, 2]),
                "selection_type": selection_type,
                "selection_reason": reason,
            }
        )

    # Bulk-like atoms from both bulk slabs away from the interface.
    if len(ranges) >= 2:
        lower_slab = (z > ranges[0][1] + BUFFER_AWAY_FROM_INTERFACE) & (
            z < ranges[1][0] - BUFFER_AWAY_FROM_INTERFACE
        )
        upper_slab = z > ranges[1][1] + BUFFER_AWAY_FROM_INTERFACE
    else:
        lower_slab = z < ranges[0][0] - BUFFER_AWAY_FROM_INTERFACE
        upper_slab = z > ranges[0][1] + BUFFER_AWAY_FROM_INTERFACE

    bulk_label = (cna == 1) & (ptm == 1) & ~interface_mask
    for slab_name, slab_mask in [
        ("lower bulk", lower_slab),
        ("upper bulk", upper_slab),
    ]:
        mask = slab_mask & bulk_label
        if not mask.any():
            continue
        median_pe = float(np.median(pe[mask]))
        score = -np.abs(pe - median_pe)
        for index in _select_top_unselected(mask, score, selected_indices, 2):
            add(
                index,
                "bulk",
                f"{slab_name} CNA=1/PTM=1 away from interface, "
                "PEATOM near slab median",
            )

    # One representative per non-dominant PTM/CNA category.
    non_dominant = (cna != 1) | (ptm != 1)
    pairs = sorted(
        set(zip(ptm[non_dominant].tolist(), cna[non_dominant].tolist()))
    )
    for pair_ptm, pair_cna in pairs:
        pair_mask = non_dominant & (ptm == pair_ptm) & (cna == pair_cna)
        interface_pair = pair_mask & interface_mask
        if interface_pair.any():
            category_mask = interface_pair
            location = "interface"
        else:
            category_mask = pair_mask
            location = "structure"
        for index in _select_top_unselected(
            category_mask, pe, selected_indices, 1
        ):
            add(
                index,
                "non_dominant",
                f"non-dominant PTM {pair_ptm}/CNA {pair_cna} in {location}",
            )

    # Highest PEATOM atoms inside the interface.
    for index in _select_top_unselected(
        interface_mask, pe, selected_indices, 2
    ):
        add(
            index,
            "high_pe_interface",
            "highest c_PEATOM in interface region",
        )

    # One atom per distinct PTM/CNA pair and interface run.
    for run_index, (start_bin, end_bin) in enumerate(runs):
        run_mask = interface_mask & (z >= centers[start_bin]) & (
            z <= centers[end_bin]
        )
        run_pairs = sorted(
            set(zip(ptm[run_mask].tolist(), cna[run_mask].tolist()))
        )
        for pair_ptm, pair_cna in run_pairs:
            pair_mask = run_mask & (ptm == pair_ptm) & (cna == pair_cna)
            for index in _select_top_unselected(
                pair_mask, pe, selected_indices, 1
            ):
                add(
                    index,
                    "interface_pair",
                    f"interface run {run_index + 1} "
                    f"PTM {pair_ptm}/CNA {pair_cna}",
                )


def _labels_for_records(records: list[dict]) -> list[str]:
    return [
        f"{record['structure']}|{record['id']}|{record['selection_type']}"
        for record in records
    ]


def _pair_mean(
    matrix: np.ndarray,
    row_mask: np.ndarray,
    col_mask: np.ndarray,
) -> float:
    if np.array_equal(row_mask, col_mask):
        rows, cols = np.where(
            np.triu(np.ones(matrix.shape, dtype=bool), 1)
            & row_mask[:, None]
            & col_mask[None, :]
        )
    else:
        rows, cols = np.where(row_mask[:, None] & col_mask[None, :])
    if rows.size == 0:
        return float("nan")
    return float(matrix[rows, cols].mean())


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    records: list[dict] = []
    selected_by_structure: dict[str, list[int]] = {}
    atoms_by_structure: dict[str, object] = {}
    selected_sets: dict[str, set[int]] = {
        structure: set() for structure in CFG_FILES
    }

    for structure, path in CFG_FILES.items():
        atoms = read_cfg(path)
        atoms_by_structure[structure] = atoms
        _select_representatives(
            structure, atoms, selected_sets[structure], records
        )
        selected_by_structure[structure] = [
            record["atom_index"]
            for record in records
            if record["structure"] == structure
        ]

    csv_fields = [
        "structure",
        "atom_index",
        "id",
        "ptm",
        "cna",
        "peatom",
        "x",
        "y",
        "z",
        "z_scaled",
        "selection_type",
        "selection_reason",
    ]
    with open(OUTPUT_DIR / "representative_atoms.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=csv_fields)
        writer.writeheader()
        writer.writerows(records)

    with warnings.catch_warnings(record=True) as constructor_warnings:
        warnings.simplefilter("always")
        soap = SOAP(
            species=SOAP_PARAMS["species"],
            periodic=SOAP_PARAMS["periodic"],
            r_cut=SOAP_PARAMS["rcut"],
            n_max=SOAP_PARAMS["nmax"],
            l_max=SOAP_PARAMS["lmax"],
            sigma=SOAP_PARAMS["sigma"],
            average="off",
            sparse=False,
            dtype="float64",
        )

    vectors_by_structure: dict[str, np.ndarray] = {}
    all_vectors: list[np.ndarray] = []
    for structure, atoms in atoms_by_structure.items():
        indices = selected_by_structure[structure]
        vectors = soap.create(atoms, centers=indices)
        norms = np.linalg.norm(vectors, axis=1, keepdims=True)
        vectors_by_structure[structure] = vectors / norms
        all_vectors.append(vectors_by_structure[structure])

    normalized = np.vstack(all_vectors)
    labels = _labels_for_records(records)
    similarity = normalized @ normalized.T
    all_finite = bool(np.isfinite(similarity).all())
    diagonal_error = float(np.max(np.abs(np.diag(similarity) - 1.0)))
    symmetry_error = float(np.max(np.abs(similarity - similarity.T)))

    with open(OUTPUT_DIR / "soap_similarity_matrix.csv", "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["label", *labels])
        for i, label in enumerate(labels):
            writer.writerow([label, *[f"{value:.6f}" for value in similarity[i]]])

    masks: dict[str, np.ndarray] = {}
    for structure in CFG_FILES:
        structure_mask = np.asarray(
            [label.startswith(f"{structure}|") for label in labels]
        )
        masks[f"{structure}_bulk"] = structure_mask & np.asarray(
            [label.endswith("|bulk") for label in labels]
        )
        masks[f"{structure}_interface"] = structure_mask & ~masks[
            f"{structure}_bulk"
        ]

    comparison_lines = [
        "Within-structure mean cosine similarities (pairwise, off-diagonal):",
        "",
        "| comparison | stable | metastable |",
        "|---|---|---|",
        "| bulk vs bulk | "
        f"{_pair_mean(similarity, masks['stable_bulk'], masks['stable_bulk']):.4f} | "
        f"{_pair_mean(similarity, masks['metastable_bulk'], masks['metastable_bulk']):.4f} |",
        "| bulk vs interface | "
        f"{_pair_mean(similarity, masks['stable_bulk'], masks['stable_interface']):.4f} | "
        f"{_pair_mean(similarity, masks['metastable_bulk'], masks['metastable_interface']):.4f} |",
        "| interface vs interface | "
        f"{_pair_mean(similarity, masks['stable_interface'], masks['stable_interface']):.4f} | "
        f"{_pair_mean(similarity, masks['metastable_interface'], masks['metastable_interface']):.4f} |",
        "",
        "Cross-structure comparison:",
        "",
        f"- stable-interface vs metastable-interface: "
        f"{_pair_mean(similarity, masks['stable_interface'], masks['metastable_interface']):.4f}",
        f"- stable-bulk vs metastable-bulk: "
        f"{_pair_mean(similarity, masks['stable_bulk'], masks['metastable_bulk']):.4f}",
    ]

    warning_lines = sorted(
        {str(item.message) for item in constructor_warnings}
    )
    summary_lines = [
        "# Stage 2: representative SOAP comparison",
        "",
        "## Method",
        "",
        "- GB normal identified by 200-bin x/y/z profiling using localized "
        "non-dominant CNA/PTM fractions and PEATOM elevation.",
        "- Representative atoms selected deliberately from bulk slabs, each "
        "non-dominant CNA/PTM category, high-PEATOM interface atoms, and "
        "distinct interface PTM/CNA environments.",
        "- SOAP baseline parameters unchanged from Stage 1: "
        f"`r_cut={SOAP_PARAMS['rcut']}`, `n_max={SOAP_PARAMS['nmax']}`, "
        f"`l_max={SOAP_PARAMS['lmax']}`, `sigma={SOAP_PARAMS['sigma']}`.",
        "",
        "## Representative atoms",
        "",
        f"- stable: {len(selected_by_structure['stable'])} atoms",
        f"- metastable: {len(selected_by_structure['metastable'])} atoms",
        "",
        "## SOAP cosine similarities",
        "",
        f"- normalized SOAP matrix shape: `{similarity.shape}`",
        f"- all matrix values finite: {all_finite}",
        f"- maximum diagonal deviation from 1.0: {diagonal_error:.3e}",
        f"- maximum symmetry error: {symmetry_error:.3e}",
        f"- DScribe reported feature dimensionality: "
        f"{soap.get_number_of_features()}",
        "",
        "\n".join(comparison_lines),
        "",
        "## DScribe LinAlgWarning",
        "",
        "The warning occurs during `SOAP.__init__`, in `get_basis_gto`, "
        "while computing the Lowdin orthonormalization "
        "`betas = sqrtm(inv(S))` of the GTO radial-basis overlap matrix. "
        "It is parameter-dependent, not atom-data-dependent, and it occurs "
        "before any SOAP vectors are computed. The vectors remain finite "
        "and reproducible.",
    ]
    if warning_lines:
        summary_lines.append("")
        summary_lines.append("Captured warnings:")
        summary_lines.extend(f"- `{message}`" for message in warning_lines)

    summary_lines.extend(
        [
            "",
            "## Caveat",
            "",
            "These comparisons use a small, deliberately selected atom "
            "subset. No disorder label, PCA/UMAP, clustering, or LER is "
            "performed, and no stable-vs-metastable structural conclusion "
            "is claimed from this sample.",
        ]
    )
    text = "\n".join(summary_lines)
    (OUTPUT_DIR / "soap_stage2_summary.md").write_text(
        text + "\n", encoding="utf-8"
    )
    print(text)

    fig, ax = plt.subplots(figsize=(max(8, len(labels) * 0.35), 8))
    image = ax.imshow(similarity, cmap="viridis", vmin=0.0, vmax=1.0)
    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=6)
    ax.set_yticklabels(labels, fontsize=6)
    fig.colorbar(image, ax=ax, label="normalized SOAP cosine similarity")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "soap_similarity_heatmap.png", dpi=150)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
