"""5. Overlap removal only affects cross-grain interface atoms."""
import pytest

from generator.builder import build_bicrystal


@pytest.fixture(scope="module")
def full():
    return build_bicrystal()


def test_removal_log_has_required_fields(full):
    required = {"atom_id", "grain_id", "x", "y", "z",
                "nearest_distance_A", "gb", "partner_atom_id"}
    for rec in full.removal_log:
        assert required.issubset(set(rec.keys()))
        assert rec["gb"] in ("GB1", "GB2")
        assert rec["grain_id"] in (1, 2)


def test_removed_atoms_are_cross_grain_pairs(full):
    # per GB the removed atoms come 30/30 from each grain (mirror balance)
    from collections import Counter
    per_gb = Counter((rec["gb"], rec["grain_id"]) for rec in full.removal_log)
    assert per_gb == {("GB1", 1): 30, ("GB1", 2): 30,
                      ("GB2", 1): 30, ("GB2", 2): 30}


def test_no_interior_atom_removed(full):
    Lx, Ly, Lz = full.box
    for rec in full.removal_log:
        # interior means more than 1 A away from both GB planes
        d1 = min(abs(rec["x"]), abs(rec["x"] - Lx))
        d2 = abs(rec["x"] - Lx / 2)
        assert d1 < 1.0 or d2 < 1.0
