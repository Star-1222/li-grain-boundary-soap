"""Li_Bicrystal_Generator: A reusable framework for generating BCC lithium bicrystal models."""

__version__ = "1.0.0"

from .builder import BiCrystalBuilder
from .config import (
    GRAIN_BOUNDARIES,
    GrainBoundary,
    get_default_config,
)
from .exporter import export_lammps, export_xyz, export_cif, export_grainid_xyz, export
from .geometry import remove_overlapping_atoms
from .utils import identify_boundary_atoms