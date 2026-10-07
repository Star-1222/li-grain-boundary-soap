"""1. Conventional BCC cell construction."""
import math

import numpy as np

from generator.config import DEFAULT_A
from generator.geometry import conventional_bcc_basis, box_dims


def test_basis_has_two_sites():
    basis = conventional_bcc_basis()
    assert basis.shape == (2, 3)
    assert np.allclose(basis[0], [0.0, 0.0, 0.0])
    assert np.allclose(basis[1], [0.5, 0.5, 0.5])


def test_density_of_conventional_cell():
    a = DEFAULT_A
    basis = conventional_bcc_basis()
    # one conventional cubic cell holds exactly 2 BCC atoms
    density = len(basis) / (a ** 3)
    assert abs(density - 2.0 / a ** 3) < 1e-15


def test_nn_distance_in_conventional_cell():
    a = DEFAULT_A
    basis = conventional_bcc_basis()
    # nearest neighbours of a BCC lattice sit at a*sqrt(3)/2
    d = math.sqrt(3.0) / 2.0 * a
    assert abs(d - 3.039749167283379) < 1e-9


def test_box_dims_formula():
    a = DEFAULT_A
    Lx, Ly, Lz = box_dims(a, 6, 3, 10)
    assert abs(Lx - 6 * math.sqrt(10) * a) < 1e-9
    assert abs(Ly - 3 * math.sqrt(10) * a) < 1e-9
    assert abs(Lz - 10 * a) < 1e-9
