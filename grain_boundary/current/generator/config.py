"""Sigma5 (310)[001] symmetric tilt grain boundary -- construction parameters.

All crystallographic choices are made explicit here so that the geometry can
be audited without running any external tool.

* BCC Li, conventional cubic cell with lattice constant ``a`` (default 3.506186475312 A,
  from the matched bulk reference for the ZuoChenLi SNAP potential).
* Tilt axis [001], misorientation angle theta = 2*atan(1/3) = 36.86989765 deg
  (Sigma5), boundary plane {310} (a mirror plane of the bicrystal).
* Orientation matrices (rows = crystal axes expressed in the lab frame):

      M1 = Rz(+theta/2)     grain 1, [100] tilted +18.4349 deg from GB normal
      M2 = Rz(-theta/2)     grain 2, [100] tilted -18.4349 deg (mirror-related)

  so the relative rotation is  R = M2 . M1^T = Rz(-theta)  about [001].
* Orthogonal periodic box (no tilt factors):

      Lx = K * sqrt(10) * a     GB normal, two half-grains of Lx/2 each
      Ly = J * sqrt(10) * a     in-plane, J CSL periods along [130]
      Lz = L * a                in-plane, L periods along [001]

  Defaults K=6, J=6, L=20 -> 66.53 x 66.53 x 70.12 A3, 14400 atoms after
  removal of the 480 coincident GB-plane atoms (see builder.py).
"""
from __future__ import annotations

import math

SIGMA = 5
TILT_AXIS: tuple[int, int, int] = (0, 0, 1)  # [001]
GB_PLANE_LABEL = "310"  # {310} family plane
GB_PLANE_MILLER_1 = "[3,-1,0]"  # GB normal in crystal-1 frame
GB_PLANE_MILLER_2 = "[3,1,0]"  # GB normal in crystal-2 frame
THETA_DEG = 36.86989764584402  # 2*atan(1/3)
THETA = 2.0 * math.atan(1.0 / 3.0)

DEFAULT_A = 3.506186475312  # Angstrom, matched bulk reference
DEFAULT_K, DEFAULT_J, DEFAULT_L = 6, 6, 20  # box periods

MASS_LI = 6.941  # g/mol
ELEMENT = "Li"

# OpenKIM model used as the "tentative" potential for the current LAMMPS stages.
KIM_MODEL = "SNAP_ZuoChenLi_2019_Li__MO_732106099012_000"

# Legacy (broken) pipeline facts, used only for the explicit regression check.
LEGACY_ATOMS_BEFORE = 3600
LEGACY_ATOMS_AFTER = 3224
LEGACY_REMOVAL_PCT = 100.0 * (3600 - 3224) / 3600  # 10.444...


def sqrt10() -> float:
    return math.sqrt(10.0)
