"""Utility and validation helpers."""

from __future__ import annotations

import sys

from .config import get_gb_by_name, get_gb_names, GRAIN_BOUNDARIES


def validate_config(config: dict) -> list[str]:
    """Check that *config* contains all required keys.

    Returns a list of error messages (empty if the config is valid).
    """
    errors: list[str] = []
    required = [
        "material", "crystal", "lattice_constant",
        "supercell", "overlap_cutoff",
    ]
    for key in required:
        if key not in config:
            errors.append(f"Missing required config key: '{key}'")
    return errors


def list_available_gbs() -> None:
    """Print a formatted table of all available grain boundaries."""
    print(f"\n{'Name':<20} {'Sigma':<8} {'Axis':<12} {'Angle':<10} {'Plane':<12}")
    print("-" * 62)
    for gb in GRAIN_BOUNDARIES:
        axis_str = f"[{gb['axis'][0]}{gb['axis'][1]}{gb['axis'][2]}]"
        plane_str = f"({gb['plane'][0]}{gb['plane'][1]}{gb['plane'][2]})"
        print(f"{gb['name']:<20} {gb['sigma']:<8} {axis_str:<12} "
              f"{gb['angle']:<10.2f} {plane_str:<12}")
    print()


def identify_boundary_atoms(atoms) -> list[int]:
    """(Placeholder) Identify atoms belonging to the grain boundary region.

    This function is designed as an extension point for future work.
    When fully implemented, it should return the indices of atoms that
    are structurally distinct from the bulk (i.e. the GB region).

    Possible detection strategies (to be implemented later):
        - Coordination number analysis
        - Common Neighbor Analysis (CNA) / Polyhedral Template
          Matching (PTM)
        - SOAP similarity against a reference bulk environment
        - Local structural disorder / centrosymmetry parameter

    Parameters
    ----------
    atoms : Atoms
        Atomic structure (expected to have a ``grain_id`` array).

    Returns
    -------
    list[int]
        Indices of boundary atoms.

    Raises
    ------
    NotImplementedError
        Always — this is a placeholder for future implementation.
    """
    raise NotImplementedError(
        "identify_boundary_atoms is not yet implemented.\n"
        "Future versions will detect GB atoms using methods such as:\n"
        "  - coordination number analysis\n"
        "  - Common Neighbor Analysis (CNA)\n"
        "  - Polyhedral Template Matching (PTM)\n"
        "  - SOAP similarity (requires DScribe or quippy)\n"
        "  - centrosymmetry / local disorder parameter\n"
        "Implement the detection logic in this function when ready."
    )


def resolve_gb(name: str) -> dict | None:
    """Resolve a GB name and exit with a helpful message if not found."""
    gb = get_gb_by_name(name)
    if gb is not None:
        return gb
    print(f"Error: grain boundary '{name}' not found in database.")
    print(f"Available grain boundaries: {', '.join(get_gb_names())}")
    sys.exit(1)