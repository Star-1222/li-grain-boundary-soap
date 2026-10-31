"""4. Grain-boundary periodicity: y and z in-plane periods match."""
import math

import pytest

from generator.builder import build_bicrystal
from generator.config import DEFAULT_A
from generator.geometry import (lattice_indices_grain1, lattice_indices_grain2,
                                plane_spacing)


@pytest.fixture(scope="module")
def full():
    return build_bicrystal()


def test_in_plane_periods_values(full):
    a = full.a
    assert abs(math.sqrt(10) * a - 11.099594587191012) < 1e-9
    assert abs(a - 3.51) < 1e-12


def test_all_atoms_on_both_lattices(full):
    a = full.a
    d = plane_spacing(a)
    miss1 = miss2 = 0
    for at in full.atoms:
        if at.grain_id == 1:
            if lattice_indices_grain1(at.x, at.y, at.z, a, d) is None:
                miss1 += 1
        else:
            if lattice_indices_grain2(at.x, at.y, at.z, a, d) is None:
                miss2 += 1
    assert miss1 == 0
    assert miss2 == 0


def test_y_periodicity_wrap(full):
    # translating any grain-1 atom by the y-period must land on another site
    a = full.a
    d = plane_spacing(a)
    y_period = math.sqrt(10) * a
    g1 = [at for at in full.atoms if at.grain_id == 1]
    hits = 0
    for at in g1:
        y2 = (at.y + y_period) % full.box[1]
        if lattice_indices_grain1(at.x, y2, at.z, a, d) is not None:
            hits += 1
    assert hits == len(g1)
