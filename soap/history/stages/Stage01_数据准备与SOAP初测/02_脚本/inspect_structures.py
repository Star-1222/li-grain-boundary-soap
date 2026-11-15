"""Validate the CFG parser and summarize both Ni grain-boundary structures."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from cfg_reader import auxiliary_names, read_cfg


PROJECT_ROOT = Path(__file__).resolve().parent
CFG_FILES = {
    "stable": PROJECT_ROOT / "ni_stable.cfg",
    "metastable": PROJECT_ROOT / "ni_metastable.cfg",
}
OUTPUT_DIR = PROJECT_ROOT / "outputs"


def _format_unique(values: np.ndarray) -> str:
    unique = np.unique(values)
    if unique.size <= 12:
        return ", ".join(str(int(v)) for v in unique)
    return f"{unique.size} unique values; first 12: " + ", ".join(
        str(int(v)) for v in unique[:12]
    )


def _field_summary(name: str, values: np.ndarray) -> str:
    if not np.all(np.isfinite(values)):
        return f"{name}: contains non-finite values"
    if np.all(values == np.floor(values)):
        return (
            f"{name}: integer field; unique {_format_unique(values)}; "
            f"min={int(values.min())}, max={int(values.max())}"
        )
    return (
        f"{name}: numeric field; min={values.min():.6g}, "
        f"max={values.max():.6g}, unique={np.unique(values).size}"
    )


def inspect_one(path: Path, label: str) -> list[str]:
    lines: list[str] = []
    atoms = read_cfg(path)

    scaled = atoms.get_scaled_positions()
    cartesian = atoms.positions
    cell = atoms.cell[:]
    lengths = atoms.cell.lengths()

    lines.append(f"=== {label}: {path.name} ===")
    lines.append(
        "header atom count: "
        f"{atoms.info.get('cfg_natoms_header', 'not recorded')}"
    )
    lines.append(f"parsed atom count: {len(atoms)}")
    lines.append(f"elements: {sorted(set(atoms.get_chemical_symbols()))}")
    lines.append(f"source: {atoms.info.get('cfg_source', 'not recorded')}")
    lines.append(
        "entry_count: "
        f"{atoms.info.get('entry_count', 'not recorded')} "
        "(3 coordinates + auxiliary fields)"
    )
    lines.append(
        "coordinate convention used for ASE Atoms: "
        f"{atoms.info.get('cfg_coordinate_convention', 'not recorded')}"
    )
    lines.append(
        "cell (H0 rows, Angstrom):\n"
        + "\n".join("  [" + ", ".join(f"{v:.6f}" for v in row) + "]" for row in cell)
    )
    lines.append(
        "cell lengths (Angstrom): "
        + ", ".join(f"{v:.6f}" for v in lengths)
    )
    lines.append(
        "scaled coordinate ranges: "
        + ", ".join(
            f"{axis}=[{scaled[:, axis].min():.8g}, {scaled[:, axis].max():.8g}]"
            for axis in range(3)
        )
    )
    lines.append(
        "cartesian coordinate ranges (Angstrom): "
        + ", ".join(
            f"{axis}=[{cartesian[:, axis].min():.6f}, {cartesian[:, axis].max():.6f}]"
            for axis in range(3)
        )
    )
    lines.append(
        "scaled-coordinate sanity check: "
        + (
            "PASS (all scaled coordinates in [0, 1] within 1e-8)"
            if np.all(scaled >= -1e-8) and np.all(scaled <= 1.0 + 1e-8)
            else "NOTE (some scaled coordinates outside [0, 1])"
        )
    )
    lines.append(
        "ASE pbc: "
        + ", ".join("True" if flag else "False" for flag in atoms.pbc)
    )
    lines.append(
        "PBC consistency: "
        + ("PASS" if np.all(atoms.pbc) else "FAIL")
        + " (ASE pbc = ["
        + ", ".join("True" if flag else "False" for flag in atoms.pbc)
        + "])"
    )

    names = auxiliary_names(atoms)
    lines.append(f"auxiliary field names ({len(names)}): {', '.join(names)}")
    lines.append("auxiliary field summaries:")
    values = atoms.arrays["auxiliary_values"]
    for i, name in enumerate(names):
        lines.append("  " + _field_summary(name, values[:, i]))

    return lines


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_lines: list[str] = []

    for label, path in CFG_FILES.items():
        if not path.exists():
            print(f"Missing input: {path}", file=sys.stderr)
            return 1
        summary_lines.extend(inspect_one(path, label))
        summary_lines.append("")

    text = "\n".join(summary_lines)
    print(text)
    (OUTPUT_DIR / "inspection_summary.txt").write_text(
        text + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
