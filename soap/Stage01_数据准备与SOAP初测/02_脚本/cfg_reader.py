"""Reader for LAMMPS-style CFG files that produce ASE Atoms objects.

The parser works from the header fields rather than from fixed line
numbers.  It preserves the numeric per-atom auxiliary columns as a 2D
array and keeps the auxiliary names in ``atoms.info``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from ase import Atoms


_PARTICLES_RE = re.compile(r"^\s*Number\s+of\s+particles\s*=\s*(\d+)\s*$")
_H0_RE = re.compile(
    r"^\s*H0\(\s*(\d+)\s*,\s*(\d+)\s*\)\s*=\s*([^\sA]+)"
)
_ENTRY_COUNT_RE = re.compile(r"^\s*entry_count\s*=\s*(\d+)\s*$")
_AUXILIARY_RE = re.compile(
    r"^\s*auxiliary\[\s*(\d+)\s*\]\s*=\s*(\S+)\s*$"
)


@dataclass(frozen=True)
class CfgHeader:
    """Fields read from the CFG header."""

    natoms: int
    scale: float
    cell: np.ndarray
    entry_count: int
    auxiliary_names: tuple[str, ...]
    has_velocities: bool


class CfgFormatError(ValueError):
    """Raised when a file does not look like a supported CFG file."""


def _parse_float_field(text: str, pattern: re.Pattern[str], label: str) -> float:
    match = pattern.match(text)
    if match is None:
        raise CfgFormatError(f"Could not parse {label} from: {text!r}")
    return float(match.group(1))


def parse_header(lines: list[str]) -> tuple[CfgHeader, int]:
    """Parse a CFG header and return the header plus data-start line index."""

    if not lines:
        raise CfgFormatError("Empty CFG file")

    natoms_match = _PARTICLES_RE.match(lines[0])
    if natoms_match is None:
        raise CfgFormatError(
            f"First line is not 'Number of particles = N': {lines[0]!r}"
        )
    natoms = int(natoms_match.group(1))

    scale = 1.0
    cell = np.zeros((3, 3), dtype=float)
    entry_count: int | None = None
    auxiliary: dict[int, str] = {}
    has_velocities = False
    data_start: int | None = None

    for idx, line in enumerate(lines):
        stripped = line.strip()

        scale_match = re.match(r"^\s*A\s*=\s*([^\s]+)", stripped)
        if scale_match is not None:
            scale = float(scale_match.group(1))

        h0_match = _H0_RE.match(stripped)
        if h0_match is not None:
            i = int(h0_match.group(1)) - 1
            j = int(h0_match.group(2)) - 1
            if not (0 <= i < 3 and 0 <= j < 3):
                raise CfgFormatError(f"Unexpected H0 index in: {stripped!r}")
            cell[i, j] = float(h0_match.group(3))

        if stripped.startswith(".") and "VELOCITY" in stripped.upper():
            has_velocities = ".NO_VELOCITY." not in stripped.upper()

        entry_match = _ENTRY_COUNT_RE.match(stripped)
        if entry_match is not None:
            entry_count = int(entry_match.group(1))

        aux_match = _AUXILIARY_RE.match(stripped)
        if aux_match is not None:
            auxiliary[int(aux_match.group(1))] = aux_match.group(2)

        if data_start is None and entry_count is not None and auxiliary:
            expected_last_index = entry_count - 3 - 1
            if max(auxiliary) == expected_last_index:
                data_start = idx + 1

    if entry_count is None:
        raise CfgFormatError("Missing 'entry_count' in CFG header")
    if data_start is None:
        raise CfgFormatError(
            "Could not locate the end of the auxiliary-name block; "
            "entry_count and auxiliary[] fields are inconsistent"
        )

    auxiliary_names = tuple(
        auxiliary[i] for i in sorted(auxiliary)
    )
    if len(auxiliary_names) != entry_count - 3:
        raise CfgFormatError(
            f"entry_count={entry_count} but parsed "
            f"{len(auxiliary_names)} auxiliary fields"
        )

    header = CfgHeader(
        natoms=natoms,
        scale=scale,
        cell=cell,
        entry_count=entry_count,
        auxiliary_names=auxiliary_names,
        has_velocities=has_velocities,
    )
    return header, data_start


def _is_fractional_coordinates(positions: np.ndarray, tolerance: float = 1e-8) -> bool:
    if positions.size == 0:
        return True
    return bool(
        np.all(positions >= -tolerance) and np.all(positions <= 1.0 + tolerance)
    )


def read_cfg(
    filename: str | Path,
    *,
    coordinate_convention: str = "auto",
) -> Atoms:
    """Read a CFG file into an ASE Atoms object.

    Parameters
    ----------
    filename:
        Path to the CFG file.
    coordinate_convention:
        ``"scaled"``, ``"cartesian"``, or ``"auto"``.  ``"auto"`` uses
        scaled coordinates when every coordinate lies in ``[0, 1]``.

    Returns
    -------
    Atoms
        An ASE structure.  Per-atom auxiliary columns are stored in
        ``atoms.arrays["auxiliary_values"]`` and their names in
        ``atoms.info["auxiliary_names"]``.
    """

    path = Path(filename)
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    header, data_start = parse_header(lines)

    natoms = header.natoms
    record_lines = 3 if not header.has_velocities else 6
    expected_data_lines = natoms * record_lines
    available_data_lines = len(lines) - data_start
    if available_data_lines != expected_data_lines:
        raise CfgFormatError(
            f"Expected {expected_data_lines} data lines for {natoms} atoms, "
            f"found {available_data_lines}"
        )

    elements: list[str] = []
    masses: list[float] = []
    values: list[list[float]] = []
    for atom_index in range(natoms):
        base = data_start + atom_index * record_lines
        try:
            mass = float(lines[base])
            element = lines[base + 1].strip()
            values_line = lines[base + 2]
        except IndexError as exc:
            raise CfgFormatError(
                f"Incomplete record for atom {atom_index + 1}"
            ) from exc
        if not element:
            raise CfgFormatError(
                f"Empty element symbol for atom {atom_index + 1}"
            )

        record = [float(token) for token in values_line.split()]
        if len(record) != header.entry_count:
            raise CfgFormatError(
                f"Atom {atom_index + 1} has {len(record)} numeric entries; "
                f"expected {header.entry_count}"
            )

        elements.append(element)
        masses.append(mass)
        values.append(record)

    values_array = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(values_array)):
        raise CfgFormatError("Non-finite numeric entries found in CFG data")

    if coordinate_convention not in {"auto", "scaled", "cartesian"}:
        raise ValueError(
            "coordinate_convention must be 'auto', 'scaled', or 'cartesian'"
        )

    fractional = values_array[:, :3]
    if coordinate_convention == "cartesian" or (
        coordinate_convention == "auto" and not _is_fractional_coordinates(fractional)
    ):
        positions = fractional * header.scale
    else:
        positions = fractional @ header.cell

    atoms = Atoms(
        symbols=elements,
        positions=positions,
        cell=header.cell,
        pbc=True,
    )
    atoms.set_masses(np.asarray(masses, dtype=float))

    auxiliary_values = values_array[:, 3:]
    atoms.arrays["auxiliary_values"] = auxiliary_values
    atoms.info["auxiliary_names"] = list(header.auxiliary_names)
    atoms.info["entry_count"] = header.entry_count
    atoms.info["cfg_natoms_header"] = header.natoms
    atoms.info["cfg_source"] = str(path.resolve())
    used_fractional = coordinate_convention != "cartesian" and (
        coordinate_convention == "scaled"
        or _is_fractional_coordinates(fractional)
    )
    atoms.info["cfg_coordinate_convention"] = (
        "scaled" if used_fractional else "cartesian"
    )
    return atoms


def get_auxiliary(atoms: Atoms, name: str) -> np.ndarray:
    """Return one auxiliary column by its CFG field name."""

    names = atoms.info.get("auxiliary_names")
    if not names:
        raise KeyError("Atoms object has no stored auxiliary field names")
    values = atoms.arrays["auxiliary_values"]
    try:
        index = names.index(name)
    except ValueError as exc:
        raise KeyError(f"Unknown auxiliary field {name!r}") from exc
    return values[:, index]


def auxiliary_names(atoms: Atoms) -> list[str]:
    """Return the auxiliary field names stored on an Atoms object."""

    return list(atoms.info.get("auxiliary_names", []))


def list_auxiliary_fields(atoms: Atoms) -> list[tuple[str, np.ndarray]]:
    """Return ``(name, column)`` pairs for all auxiliary fields."""

    values = atoms.arrays["auxiliary_values"]
    return [
        (name, values[:, i])
        for i, name in enumerate(auxiliary_names(atoms))
    ]


__all__ = [
    "CfgFormatError",
    "CfgHeader",
    "auxiliary_names",
    "get_auxiliary",
    "list_auxiliary_fields",
    "parse_header",
    "read_cfg",
]
