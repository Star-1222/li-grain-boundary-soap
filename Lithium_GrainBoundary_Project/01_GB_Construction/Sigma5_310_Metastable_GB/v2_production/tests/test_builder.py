"""3. Bicrystal builder: counts, mirror symmetry, GB equivalence."""
import math

import pytest

from generator.builder import COINC_TOL, GB_TOL, build_bicrystal
from generator.config import DEFAULT_A, DEFAULT_J, DEFAULT_K, DEFAULT_L, THETA_DEG
from generator.geometry import pbc_dist


@pytest.fixture(scope="module")
def full():
    return build_bicrystal()


def test_full_model_counts(full):
    assert full.meta["sites_per_grain"] == 1860
    assert full.meta["total_before_removal"] == 3720
    assert full.meta["coincident_pairs_gb1"] == 60
    assert full.meta["coincident_pairs_gb2"] == 60
    assert full.meta["removed_total"] == 120
    assert full.meta["final_total"] == 3600
    assert full.meta["final_grain1"] == 1800
    assert full.meta["final_grain2"] == 1800
    assert len(full.atoms) == 3600


def test_mirror_symmetry(full):
    Lx, Ly, Lz = full.box
    # every non-GB atom must have an exact mirror partner in the other grain.
    # The only exceptions are the 120 atoms that survived the removal of their
    # coincident twin on a GB plane (removed-twin survivors).
    g1 = {(round(a.x, 6), round(a.y, 6), round(a.z, 6)) for a in full.atoms if a.grain_id == 1}
    g2 = {(round(a.x, 6), round(a.y, 6), round(a.z, 6)) for a in full.atoms if a.grain_id == 2}
    unmatched = 0
    for a in full.atoms:
        key = (round((Lx - a.x) % Lx, 6), round(a.y, 6), round(a.z, 6))
        partner_set = g2 if a.grain_id == 1 else g1
        if key not in partner_set:
            unmatched += 1
    assert unmatched == 120, f"expected 120 unmatched mirror twins, got {unmatched}"


def test_removed_atoms_only_on_gb_planes(full):
    Lx, Ly, Lz = full.box
    assert len(full.removal_log) == 120
    for rec in full.removal_log:
        assert abs(rec["x"]) < GB_TOL or abs(rec["x"] - Lx / 2) < GB_TOL
        assert rec["nearest_distance_A"] < COINC_TOL


def test_no_overlaps_remain(full):
    Lx, Ly, Lz = full.box
    pos = [(a.x, a.y, a.z) for a in full.atoms]
    # brute force on the plane subset + cell list would be overkill here;
    # check the known global minimum equals 2d
    from generator.geometry import CellList
    cl = CellList(pos, Lx, Ly, Lz, 3.0)
    gmin = cl.min_distance()
    assert abs(gmin - 2 * DEFAULT_A / math.sqrt(10)) < 1e-6


def test_two_gbs_equivalent(full):
    Lx, Ly, Lz = full.box
    from generator.verify import _cross_gb_spectra
    s1, s2 = _cross_gb_spectra(full)
    assert len(s1) == len(s2) == 900
    assert s1 == pytest.approx(s2, abs=1e-9)
