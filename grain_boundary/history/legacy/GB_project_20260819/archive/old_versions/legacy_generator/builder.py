"""Core bicrystal builder using ASE + NumPy.

Builds a bicrystal from a grain boundary descriptor by:

1. Creating an oriented unit cell (GB normal along x, tilt axis along z).
2. Replicating to form a supercell.
3. Splitting into left and right grains.
4. Rotating each grain by ±θ/2 around the tilt axis.
5. Reassembling and removing overlapping atoms.
6. Centering the GB in the cell.

Usage is intentionally decoupled from any specific GB — pass the
parameters and get back an ``ase.Atoms`` object.
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
from ase import Atoms
from ase.build import bulk, cut

from .config import get_default_config
from .geometry import (
    center_gb,
    orientation_matrix,
    remove_overlapping_atoms,
    rotation_matrix_around_axis,
)

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════
#  GrainBoundaryDescriptor  (lightweight parameter container)
# ═══════════════════════════════════════════════════════════════════════

class GrainBoundaryDescriptor:
    """Carries the geometric parameters of a single grain boundary.

    This is the internal representation; construct from a database dict
    via ``cls.from_dict()`` or create one directly.
    """

    def __init__(
        self,
        name: str,
        sigma: int,
        axis: tuple[int, int, int],
        angle: float,
        plane: tuple[int, int, int],
        description: str = "",
    ) -> None:
        self.name = name
        self.sigma = sigma
        self.axis = axis
        self.angle = angle  # degrees
        self.plane = plane
        self.description = description

    # ── factory ─────────────────────────────────────────────────────

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> GrainBoundaryDescriptor:
        """Build from a dictionary (e.g. from the database in config.py)."""
        return cls(
            name=d["name"],
            sigma=d["sigma"],
            axis=d["axis"],
            angle=d["angle"],
            plane=d["plane"],
            description=d.get("description", ""),
        )

    def __repr__(self) -> str:
        return (
            f"GrainBoundaryDescriptor("
            f"name={self.name!r}, sigma={self.sigma}, "
            f"axis={self.axis}, angle={self.angle}, "
            f"plane={self.plane})"
        )


# ═══════════════════════════════════════════════════════════════════════
#  BiCrystalBuilder
# ═══════════════════════════════════════════════════════════════════════

class BiCrystalBuilder:
    """Builds and exports bicrystal models.

    Parameters
    ----------
    config : dict
        Configuration dictionary (use ``config.get_default_config()``
        as a starting point and override as needed).
    """

    def __init__(self, config: dict[str, Any] | None = None) -> None:
        self.config = get_default_config() if config is None else config

        self.lattice_constant: float = self.config["lattice_constant"]
        self.crystal: str = self.config["crystal"]
        self.material: str = self.config["material"]
        self.supercell: tuple[int, int, int] = self.config["supercell"]
        self.overlap_cutoff: float = self.config["overlap_cutoff"]
        self.vacuum: float = self.config.get("vacuum", 0.0)
        self.periodic: tuple[bool, bool, bool] = self.config["periodic"]

    # ── public API ──────────────────────────────────────────────────

    def build(self, gb: GrainBoundaryDescriptor) -> Atoms:
        """Build a bicrystal for the given grain boundary descriptor.

        Parameters
        ----------
        gb : GrainBoundaryDescriptor
            Grain boundary parameters (name, sigma, axis, angle, plane).

        Returns
        -------
        atoms : Atoms
            Bicrystal structure.
        """
        logger.info("Building bicrystal: %s", gb.name)
        logger.info("  plane=%s, axis=%s, angle=%.2f deg, sigma=%d",
                     gb.plane, gb.axis, gb.angle, gb.sigma)

        # Step 1 – Create an oriented unit cell
        oriented_unit = self._create_oriented_cell(gb)
        logger.info("  Oriented unit cell: %d atoms, cell=%s",
                     len(oriented_unit),
                     oriented_unit.cell.lengths())

        # Step 2 – Replicate to form a supercell
        supercell_atoms = oriented_unit.repeat(self.supercell)
        logger.info("  Supercell: %d atoms, cell=%s",
                     len(supercell_atoms),
                     supercell_atoms.cell.lengths())

        # Step 3 – Split and rotate grains
        bicrystal = self._split_and_rotate(supercell_atoms, gb)
        logger.info("  After split/rotate: %d atoms", len(bicrystal))

        # Step 4 – Remove overlapping atoms at the GB interface
        bicrystal = remove_overlapping_atoms(bicrystal, self.overlap_cutoff)
        logger.info("  After overlap removal: %d atoms", len(bicrystal))

        # Step 5 – Centre the GB in the cell
        bicrystal = center_gb(bicrystal)

        # Log grain_id statistics
        if "grain_id" in bicrystal.arrays:
            g1 = int((bicrystal.arrays["grain_id"] == 1).sum())
            g2 = int((bicrystal.arrays["grain_id"] == 2).sum())
            logger.info("  Grain 1: %d atoms, Grain 2: %d atoms", g1, g2)

        # Step 6 – Set periodicity
        p = self.periodic
        bicrystal.pbc = (p[0], p[1], p[2])

        # Step 7 – Add vacuum if requested
        if self.vacuum > 0.0:
            bicrystal.center(vacuum=self.vacuum, axis=0)

        logger.info("  Final: %d atoms, cell=%s",
                     len(bicrystal),
                     bicrystal.cell.lengths())
        return bicrystal

    # ── internal helpers ────────────────────────────────────────────

    def _create_oriented_cell(self, gb: GrainBoundaryDescriptor) -> Atoms:
        """Create a unit cell aligned with the GB geometry.

        The cell is oriented so that:
            x  —— GB plane normal (plane normal direction)
            y  —— in-plane, perpendicular to the tilt axis
            z  —— tilt axis (rotation axis)
        """
        h, k, l = gb.plane

        # Start with a conventional cubic BCC cell
        bulk_atoms = bulk(self.material, self.crystal, a=self.lattice_constant)

        # Use `cut` to reorient the cell.
        # The new basis vectors in terms of the old cubic basis are:
        #   a_new = (h, k, l)   — GB plane normal
        #   b_new = (-k, h, 0)  — in-plane, perpendicular to tilt axis
        #   c_new = (0, 0, 1)   — tilt axis
        oriented = cut(
            bulk_atoms,
            a=(h, k, l),
            b=(-k, h, 0),
            c=(0, 0, 1),
        )
        # `cut` may produce atoms outside the cell; wrap them back.
        oriented.wrap()
        return oriented

    def _split_and_rotate(
        self,
        atoms: Atoms,
        gb: GrainBoundaryDescriptor,
    ) -> Atoms:
        """Split the supercell at the x-midplane and rotate each half.

        The left grain is rotated by -θ/2 around the z-axis;
        the right grain is rotated by +θ/2 around the z-axis.
        """
        frac = atoms.get_scaled_positions()
        x_frac = frac[:, 0]
        mid = 0.5 * (x_frac.max() + x_frac.min())

        left_mask = x_frac <= mid
        right_mask = ~left_mask

        # ── Assign grain_id labels ────────────────────────────
        # Grain 1 = left half, Grain 2 = right half.
        # Stored as an ASE atom array: survives overlap removal,
        # centering, wrap, and subset operations automatically.
        grain_ids = np.zeros(len(atoms), dtype=int)
        grain_ids[left_mask] = 1
        grain_ids[right_mask] = 2
        atoms.new_array("grain_id", grain_ids)

        # Build rotation matrices
        theta = gb.angle  # degrees
        axis = np.array(gb.axis, dtype=float)

        R_left = rotation_matrix_around_axis(axis, -theta / 2.0, deg=True)
        R_right = rotation_matrix_around_axis(axis, theta / 2.0, deg=True)

        # Reference point for rotation: centre of the simulation box
        centre = 0.5 * (atoms.cell[0] + atoms.cell[1] + atoms.cell[2])

        pos = atoms.get_positions()

        # Rotate left grain
        pos_left = pos[left_mask] - centre
        pos[left_mask] = centre + (R_left @ pos_left.T).T

        # Rotate right grain
        pos_right = pos[right_mask] - centre
        pos[right_mask] = centre + (R_right @ pos_right.T).T

        atoms.set_positions(pos)
        atoms.wrap()
        return atoms