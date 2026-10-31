"""11. The broken legacy 3600 -> 3224 / 10.44% scenario is absent."""
from generator.builder import build_bicrystal
from generator.config import LEGACY_ATOMS_BEFORE, LEGACY_ATOMS_AFTER, LEGACY_REMOVAL_PCT
from generator.verify import check_legacy_regression


def test_legacy_constants():
    assert LEGACY_ATOMS_BEFORE == 3600
    assert LEGACY_ATOMS_AFTER == 3224
    assert abs(LEGACY_REMOVAL_PCT - 10.444444444444445) < 1e-9


def test_new_model_does_not_reproduce_legacy():
    bc = build_bicrystal()
    res = check_legacy_regression(bc)
    assert res["before"] == 3720
    assert res["after"] == 3600
    assert res["removed"] == res["expected_removed"] == 120
    assert res["legacy_bad"] is False


def test_small_model_uses_size_dependent_removal_count():
    bc = build_bicrystal(K=3, J=3, L=5)
    res = check_legacy_regression(bc)
    assert res["before"] == 960
    assert res["after"] == 900
    assert res["removed"] == res["expected_removed"] == 60
    assert abs(res["removal_pct"] - 6.25) < 1e-12
    assert res["legacy_bad"] is False
