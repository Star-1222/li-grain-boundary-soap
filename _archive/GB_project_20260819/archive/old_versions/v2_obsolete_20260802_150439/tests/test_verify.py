"""10. Full automatic verification passes on the built model."""
import pytest

from generator.builder import build_bicrystal
from generator.verify import run_all_checks


@pytest.fixture(scope="module")
def results():
    bc = build_bicrystal()
    return run_all_checks(bc)


def test_all_checks_pass(results):
    assert results["passed"] is True


def test_orientation_values(results):
    o = results["orientation"]
    assert abs(o["misorientation_angle_deg"] - 36.86989764584402) < 1e-6
    assert abs(o["cubic_symmetry_min_angle_deg"] - 36.86989764584402) < 1e-6
    assert abs(o["rotation_axis"][2]) > 0.9999


def test_counts_and_density(results):
    c = results["counts"]
    assert c["final_total"] == 3600
    assert c["final_grain1"] == c["final_grain2"] == 1800
    assert abs(c["density_atoms_per_A3"] - 2.0 / 3.51 ** 3) < 1e-12


def test_distances(results):
    d = results["distances"]
    assert abs(d["global_min_A"] - 2.219918917438202) < 1e-6
    assert abs(d["nn_bcc_A"] - 3.039749167283379) < 1e-9


def test_bulk_bcc(results):
    b = results["bulk_bcc"]
    assert b["sampled"] == 2520
    assert b["failed"] == 0


def test_two_gbs(results):
    g = results["two_gbs"]
    assert g["spectra_identical"] is True
    assert g["cross_gb_pairs_below_4p5A_gb1"] == 900
    assert g["cross_gb_pairs_below_4p5A_gb2"] == 900


def test_no_voids(results):
    v = results["voids"]
    assert v["empty_bins"] == []
    assert v["max_empty_gap_A"] < 1.0


def test_removal_counts(results):
    r = results["removal"]
    assert r["removed_total"] == 120
    assert r["removed_grain1"] == r["removed_grain2"] == 60
    assert r["gb1_removed"] == r["gb2_removed"] == 60


def test_periodicity(results):
    p = results["periodicity"]
    assert p["grain1_misses"] == 0
    assert p["grain2_misses"] == 0
