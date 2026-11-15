"""Stage 2: locate the GB normal empirically by binning atoms."""

from __future__ import annotations

import os
from pathlib import Path

from cfg_reader import get_auxiliary, read_cfg


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

N_BINS = 200
INTERFACE_FRACTION_THRESHOLD = 0.1
MIN_BIN_ATOMS = 50


def compute_axis_profile(atoms, axis: int) -> dict:
    """Bin atoms along one Cartesian axis and summarize CNA/PTM/PEATOM."""

    positions = atoms.positions[:, axis]
    cna = get_auxiliary(atoms, "c_CNA").astype(int)
    ptm = get_auxiliary(atoms, "c_PTM_1").astype(int)
    pe = get_auxiliary(atoms, "c_PEATOM")

    lo = float(positions.min())
    hi = float(positions.max())
    edges = np.linspace(lo, hi, N_BINS + 1)
    bin_indices = np.clip(
        np.digitize(positions, edges[1:-1]), 0, N_BINS - 1
    )
    centers = 0.5 * (edges[:-1] + edges[1:])

    unique_cna = np.unique(cna)
    unique_ptm = np.unique(ptm)
    cna_counts = np.zeros((N_BINS, len(unique_cna)), dtype=int)
    ptm_counts = np.zeros((N_BINS, len(unique_ptm)), dtype=int)
    bin_counts = np.zeros(N_BINS, dtype=int)
    non_dominant_counts = np.zeros(N_BINS, dtype=int)
    pe_values: list[list[float]] = [[] for _ in range(N_BINS)]
    non_dominant = (cna != 1) | (ptm != 1)

    for atom_index, bin_index in enumerate(bin_indices):
        bin_counts[bin_index] += 1
        cna_counts[bin_index, np.searchsorted(unique_cna, cna[atom_index])] += 1
        ptm_counts[bin_index, np.searchsorted(unique_ptm, ptm[atom_index])] += 1
        if non_dominant[atom_index]:
            non_dominant_counts[bin_index] += 1
        pe_values[bin_index].append(float(pe[atom_index]))

    non_dominant_fraction = np.divide(
        non_dominant_counts,
        bin_counts,
        out=np.full(N_BINS, np.nan),
        where=bin_counts > 0,
    )

    pe_mean = np.full(N_BINS, np.nan)
    pe_median = np.full(N_BINS, np.nan)
    pe_std = np.full(N_BINS, np.nan)
    for i in range(N_BINS):
        values = np.asarray(pe_values[i], dtype=float)
        if values.size:
            pe_mean[i] = values.mean()
            pe_median[i] = np.median(values)
            pe_std[i] = values.std()

    return {
        "axis": axis,
        "centers": centers,
        "bin_indices": bin_indices,
        "bin_counts": bin_counts,
        "cna_counts": cna_counts,
        "ptm_counts": ptm_counts,
        "unique_cna": unique_cna,
        "unique_ptm": unique_ptm,
        "non_dominant_fraction": non_dominant_fraction,
        "pe_mean": pe_mean,
        "pe_median": pe_median,
        "pe_std": pe_std,
    }


def find_interface_bins(
    profile: dict,
    *,
    fraction_threshold: float = INTERFACE_FRACTION_THRESHOLD,
    min_bin_atoms: int = MIN_BIN_ATOMS,
) -> tuple[np.ndarray, list[tuple[int, int]]]:
    """Return a boolean bin mask and contiguous interface-bin runs."""

    mask = (
        (profile["non_dominant_fraction"] > fraction_threshold)
        & (profile["bin_counts"] >= min_bin_atoms)
    )
    runs: list[tuple[int, int]] = []
    start: int | None = None
    for i, flag in enumerate(mask):
        if flag and start is None:
            start = i
        elif not flag and start is not None:
            runs.append((start, i - 1))
            start = None
    if start is not None:
        runs.append((start, len(mask) - 1))
    return mask, runs


def identify_gb_normal(profiles: dict[str, dict]) -> str:
    """Choose the axis with the most concentrated interface-like bins."""

    scores: dict[str, tuple[int, float]] = {}
    for axis, profile in profiles.items():
        mask, _ = find_interface_bins(profile)
        max_fraction = float(np.nanmax(profile["non_dominant_fraction"]))
        scores[axis] = (int(mask.sum()), max_fraction)
    return max(scores, key=lambda axis: scores[axis])


def _interface_runs_as_ranges(
    profile: dict, runs: list[tuple[int, int]]
) -> list[tuple[float, float]]:
    centers = profile["centers"]
    return [(float(centers[start]), float(centers[end])) for start, end in runs]


def _atom_interface_mask(profile: dict, bin_mask: np.ndarray) -> np.ndarray:
    return bin_mask[profile["bin_indices"]]


def plot_spatial_profile(
    label: str,
    profile: dict,
    runs: list[tuple[int, int]],
    output_path: Path,
) -> None:
    """Plot CNA/PTM and PEATOM profiles along the identified GB normal."""

    centers = profile["centers"]
    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

    for start, end in runs:
        for ax in axes:
            ax.axvspan(centers[start], centers[end], color="0.85", zorder=0)

    axes[0].plot(
        centers,
        profile["non_dominant_fraction"],
        color="#1f77b4",
        label="non-dominant CNA/PTM fraction",
    )
    axes[0].set_ylabel("non-dominant fraction")
    axes[0].legend(loc="upper right", fontsize=8)

    for i, value in enumerate(profile["unique_ptm"]):
        axes[1].plot(
            centers,
            profile["ptm_counts"][:, i] / np.maximum(profile["bin_counts"], 1),
            label=f"PTM {value}",
        )
    for i, value in enumerate(profile["unique_cna"]):
        axes[1].plot(
            centers,
            profile["cna_counts"][:, i] / np.maximum(profile["bin_counts"], 1),
            linestyle="--",
            label=f"CNA {value}",
        )
    axes[1].set_ylabel("fraction")
    axes[1].legend(loc="upper right", fontsize=8, ncol=2)

    axes[2].plot(centers, profile["pe_mean"], color="#d62728", label="PEATOM mean")
    axes[2].fill_between(
        centers,
        profile["pe_mean"] - profile["pe_std"],
        profile["pe_mean"] + profile["pe_std"],
        color="#d62728",
        alpha=0.2,
        label="PEATOM +/- 1 std",
    )
    axes[2].set_xlabel("position along z (Angstrom)")
    axes[2].set_ylabel("c_PEATOM")
    axes[2].legend(loc="lower right", fontsize=8)

    fig.suptitle(f"{label}: spatial profile along z")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def summarize_structure(label: str) -> list[str]:
    lines: list[str] = []
    atoms = read_cfg(CFG_FILES[label])
    profiles = {
        axis: compute_axis_profile(atoms, i) for i, axis in enumerate("xyz")
    }
    normal = identify_gb_normal(profiles)
    normal_profile = profiles[normal]
    normal_mask, normal_runs = find_interface_bins(normal_profile)
    atom_interface = _atom_interface_mask(normal_profile, normal_mask)

    cna = get_auxiliary(atoms, "c_CNA").astype(int)
    ptm = get_auxiliary(atoms, "c_PTM_1").astype(int)
    pe = get_auxiliary(atoms, "c_PEATOM")
    non_dominant = (cna != 1) | (ptm != 1)
    bulk_pe_reference = float(np.median(pe[~non_dominant]))

    lines.append(f"=== {label} ===")
    for axis in "xyz":
        mask, runs = find_interface_bins(profiles[axis])
        lines.append(
            f"axis {axis}: interface-like bins={int(mask.sum())}, "
            f"max non-dominant fraction="
            f"{np.nanmax(profiles[axis]['non_dominant_fraction']):.4f}"
        )
    lines.append(
        f"identified GB normal from binning: {normal} "
        "(localized interfacial CNA/PTM deviations)"
    )
    ranges = _interface_runs_as_ranges(normal_profile, normal_runs)
    lines.append(
        "interface-like z ranges (Angstrom): "
        + "; ".join(f"{lo:.3f}-{hi:.3f}" for lo, hi in ranges)
    )
    lines.append(
        f"interface atom count: {int(atom_interface.sum())} "
        f"({100 * atom_interface.mean():.2f}% of structure)"
    )
    lines.append(
        f"bulk PEATOM median (dominant CNA/PTM outside interface): "
        f"{bulk_pe_reference:.5f}"
    )
    if atom_interface.any():
        lines.append(
            "interface PEATOM: "
            f"median={np.median(pe[atom_interface]):.5f}, "
            f"mean={pe[atom_interface].mean():.5f}, "
            f"max={pe[atom_interface].max():.5f}"
        )
    lines.append("")

    plot_spatial_profile(
        label,
        normal_profile,
        normal_runs,
        OUTPUT_DIR / f"spatial_profile_{label}.png",
    )
    return lines


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    lines: list[str] = [
        "Stage 2 spatial profiling",
        f"bin count per axis: {N_BINS}",
        "interface-like bin definition: non-dominant CNA/PTM fraction "
        f"> {INTERFACE_FRACTION_THRESHOLD} and bin atom count >= "
        f"{MIN_BIN_ATOMS}",
        "non-dominant means CNA != 1 or PTM != 1 (the dominant class "
        "inferred from the data)",
        "",
    ]
    for label in CFG_FILES:
        lines.extend(summarize_structure(label))
    text = "\n".join(lines)
    print(text)
    (OUTPUT_DIR / "spatial_profile_summary.txt").write_text(
        text + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
