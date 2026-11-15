"""Small first DScribe SOAP test on a few atoms from each Ni CFG.

This deliberately computes SOAP only for a small subset of atoms.  The
parameters are the literature/reproduction baseline from Priedeman et
al. (2018).  DScribe's implementation is not assumed to reproduce the
historical SOAP implementation numerically.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
from ase import Atoms
from dscribe.descriptors import SOAP

from cfg_reader import read_cfg


PROJECT_ROOT = Path(__file__).resolve().parent
CFG_FILES = {
    "stable": PROJECT_ROOT / "ni_stable.cfg",
    "metastable": PROJECT_ROOT / "ni_metastable.cfg",
}
OUTPUT_DIR = PROJECT_ROOT / "outputs"

# Literature-derived SOAP parameters used by Priedeman et al. as a
# reproduction/learning baseline.  These are choices C, not assumptions
# derived from the CFG files.
SOAP_PARAMS = {
    "rcut": 3.25,
    "nmax": 12,
    "lmax": 12,
    "sigma": 0.5,
    "periodic": True,
}
SPECIES = ["Ni"]

# Small test subset only; do not run SOAP on all ~47k atoms yet.
TEST_INDICES = [0, 1, 2, 3, 4]


def _cosine_similarity(vectors: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(vectors, axis=1)
    normalized = vectors / norms[:, None]
    return normalized @ normalized.T


def run_soap_test(path: Path, label: str, soap: SOAP) -> list[str]:
    lines: list[str] = []
    atoms = read_cfg(path)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        first = soap.create(atoms, centers=TEST_INDICES)
        second = soap.create(atoms, centers=TEST_INDICES)

    lines.append(f"=== SOAP test: {label} ===")
    lines.append(
        f"parameter set: rcut={SOAP_PARAMS['rcut']} A, "
        f"nmax={SOAP_PARAMS['nmax']}, lmax={SOAP_PARAMS['lmax']}, "
        f"sigma={SOAP_PARAMS['sigma']} A, species={SPECIES}, "
        f"periodic={SOAP_PARAMS['periodic']}"
    )
    lines.append(f"test atom indices: {TEST_INDICES}")
    lines.append(f"SOAP matrix shape: {first.shape}")
    lines.append(
        "shape matches DScribe reported dimensionality "
        f"({soap.get_number_of_features()}): "
        f"{bool(first.shape[1] == soap.get_number_of_features())}"
    )
    lines.append(f"all finite: {bool(np.isfinite(first).all())}")
    lines.append(f"reproducible: {bool(np.array_equal(first, second))}")
    if not np.array_equal(first, second):
        lines.append(
            "maximum absolute reproducibility difference: "
            f"{float(np.max(np.abs(first - second))):.3e}"
        )

    similarity = _cosine_similarity(first)
    lines.append("cosine similarity matrix (rows/cols follow test indices):")
    for row in similarity:
        lines.append(
            "  " + "  ".join(f"{value:.4f}" for value in row)
        )
    lines.append(
        "off-diagonal similarity range: "
        f"{float(np.min(similarity[~np.eye(len(similarity), dtype=bool)])):.4f} - "
        f"{float(np.max(similarity[~np.eye(len(similarity), dtype=bool)])):.4f}"
    )
    unique_warnings = sorted(
        {str(item.message) for item in caught if issubclass(item.category, Warning)}
    )
    if unique_warnings:
        lines.append("DScribe/NumPy warnings:")
        for message in unique_warnings:
            lines.append(f"  {message}")
    lines.append("")

    np.savez_compressed(
        OUTPUT_DIR / f"soap_vectors_{label}.npz",
        vectors=first,
        indices=np.asarray(TEST_INDICES),
        rcut=SOAP_PARAMS["rcut"],
        nmax=SOAP_PARAMS["nmax"],
        lmax=SOAP_PARAMS["lmax"],
        sigma=SOAP_PARAMS["sigma"],
    )
    return lines


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with warnings.catch_warnings(record=True) as constructor_warnings:
        warnings.simplefilter("always")
        soap = SOAP(
            species=SPECIES,
            periodic=SOAP_PARAMS["periodic"],
            r_cut=SOAP_PARAMS["rcut"],
            n_max=SOAP_PARAMS["nmax"],
            l_max=SOAP_PARAMS["lmax"],
            sigma=SOAP_PARAMS["sigma"],
            average="off",
            sparse=False,
            dtype="float64",
        )
    expected_features = soap.get_number_of_features()

    summary_lines = [
        "DScribe SOAP stage-1 test",
        f"SOAP feature dimensionality reported by DScribe: {expected_features}",
        "",
    ]
    constructor_messages = sorted(
        {str(item.message) for item in constructor_warnings}
    )
    if constructor_messages:
        summary_lines.append("DScribe construction warnings:")
        summary_lines.extend(f"  {message}" for message in constructor_messages)
        summary_lines.append("")
    for label, path in CFG_FILES.items():
        if not path.exists():
            summary_lines.append(f"Missing input: {path}")
            continue
        summary_lines.extend(run_soap_test(path, label, soap))

    text = "\n".join(summary_lines)
    print(text)
    (OUTPUT_DIR / "soap_test_summary.txt").write_text(
        text + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
