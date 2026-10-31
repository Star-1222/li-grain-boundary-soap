"""6. grain_id preservation across exports."""
import pytest

from generator.builder import build_bicrystal
from generator.exporters import (write_extxyz, write_grain_id_map,
                                 write_lammps_data, write_xyz_with_grain_id)


@pytest.fixture(scope="module")
def small(tmp_path_factory):
    bc = build_bicrystal(K=3, J=3, L=5)
    d = tmp_path_factory.mktemp("grain_id")
    write_lammps_data(bc, d / "s.data")
    write_extxyz(bc, d / "s.xyz")
    write_xyz_with_grain_id(bc, d / "s_grain_id.xyz")
    write_grain_id_map(bc, d / "s.grain_map")
    return bc, d


def test_map_counts_match_atoms(small):
    bc, d = small
    lines = [ln for ln in (d / "s.grain_map").read_text().splitlines() if not ln.startswith("#") and ln.strip()]
    assert len(lines) == len(bc.atoms)
    assert set(int(ln.split()[1]) for ln in lines) <= {1, 2}


def test_grain_id_xyz_consistent(small):
    bc, d = small
    lines = (d / "s_grain_id.xyz").read_text().splitlines()
    assert int(lines[0]) == len(bc.atoms)
    ids = [int(ln.split()[4]) for ln in lines[2:]]
    assert set(ids) <= {1, 2}


def test_extxyz_grain_id_property(small):
    bc, d = small
    lines = (d / "s.xyz").read_text().splitlines()
    assert int(lines[0]) == len(bc.atoms)
    assert "grain_id" in lines[1]
    ids = [int(ln.split()[4]) for ln in lines[2:]]
    assert set(ids) <= {1, 2}
