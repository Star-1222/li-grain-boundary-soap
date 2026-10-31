"""Configuration parameters and grain boundary database.

All physical and numerical parameters are defined here.
No hardcoded values should appear elsewhere in the code.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any


# ── Material parameters ──────────────────────────────────────────────

MATERIAL: str = "Li"
CRYSTAL: str = "bcc"
LATTICE_CONSTANT: float = 3.51  # Angstrom (initial guess, override for specific model)

# ── Supercell and numerical parameters ────────────────────────────────

SUPERCELL: tuple[int, int, int] = (20, 12, 12)
"""Repetitions of the oriented unit cell in (x, y, z).
x: normal to GB plane (thickness direction).
y: in-plane, perpendicular to tilt axis.
z: tilt axis direction.
"""

OVERLAP_CUTOFF: float = 2.3  # Angstrom — atoms closer than this are merged

VACUUM: float = 0.0  # Angstrom of vacuum on each side (non-periodic surfaces)

PERIODIC: tuple[bool, bool, bool] = (False, True, True)
"""x is non-periodic (free surfaces), y and z are periodic (GB plane)."""

# ── Output settings ──────────────────────────────────────────────────

OUTPUT_DIR: str = "output"
OUTPUT_FORMATS: tuple[str, ...] = ("lammps", "xyz", "cif", "grainid_xyz")

# ── Reserved interface: GrainBoundary dataclass ──────────────────────

@dataclass
class GrainBoundary:
    """Represents a physically meaningful grain boundary.

    This dataclass defines the standard interface for a grain boundary.
    Future extensions (e.g., twist boundaries, other crystal systems)
    should implement this same interface.
    """
    name: str
    sigma: int
    rotation_axis: tuple[int, int, int]
    misorientation_angle: float
    gb_plane: tuple[int, int, int]
    description: str = ""


# ── Grain boundary database ──────────────────────────────────────────

# Each entry represents a physically meaningful tilt grain boundary
# in BCC, derived from the CSL (coincident site lattice) model.
#
# Fields:
#   name   — short label used for selection and file naming
#   sigma  — CSL sigma value (inverse density of coincidence sites)
#   axis   — rotation (tilt) axis in Miller indices
#   angle  — misorientation angle in degrees
#   plane  — GB plane (boundary plane) in Miller indices

GRAIN_BOUNDARIES: list[dict[str, Any]] = [
    {
        "name": "Sigma5_310",
        "sigma": 5,
        "axis": (0, 0, 1),
        "angle": 36.87,
        "plane": (3, 1, 0),
        "description": "Sigma5 [001] tilt GB, (310) plane",
    },
    {
        "name": "Sigma5_210",
        "sigma": 5,
        "axis": (0, 0, 1),
        "angle": 53.13,
        "plane": (2, 1, 0),
        "description": "Sigma5 [001] tilt GB, (210) plane",
    },
    {
        "name": "Sigma13_510",
        "sigma": 13,
        "axis": (0, 0, 1),
        "angle": 22.62,
        "plane": (5, 1, 0),
        "description": "Sigma13 [001] tilt GB, (510) plane",
    },
    {
        "name": "Sigma17_410",
        "sigma": 17,
        "axis": (0, 0, 1),
        "angle": 28.07,
        "plane": (4, 1, 0),
        "description": "Sigma17 [001] tilt GB, (410) plane",
    },
    {
        "name": "Sigma17_530",
        "sigma": 17,
        "axis": (0, 0, 1),
        "angle": 61.93,
        "plane": (5, 3, 0),
        "description": "Sigma17 [001] tilt GB, (530) plane",
    },
    {
        "name": "Sigma29_520",
        "sigma": 29,
        "axis": (0, 0, 1),
        "angle": 43.60,
        "plane": (5, 2, 0),
        "description": "Sigma29 [001] tilt GB, (520) plane",
    },
]


def get_default_config() -> dict[str, Any]:
    """Return a mutable copy of the default configuration dictionary.

    Users can override individual values after calling this function.
    """
    return {
        "material": MATERIAL,
        "crystal": CRYSTAL,
        "lattice_constant": LATTICE_CONSTANT,
        "supercell": SUPERCELL,
        "overlap_cutoff": OVERLAP_CUTOFF,
        "vacuum": VACUUM,
        "periodic": PERIODIC,
        "output_dir": OUTPUT_DIR,
        "output_formats": OUTPUT_FORMATS,
    }


def get_gb_names() -> list[str]:
    """Return the list of available grain boundary names."""
    return [gb["name"] for gb in GRAIN_BOUNDARIES]


def get_gb_by_name(name: str) -> dict[str, Any] | None:
    """Look up a grain boundary entry by its name."""
    for gb in GRAIN_BOUNDARIES:
        if gb["name"] == name:
            return gb
    return None