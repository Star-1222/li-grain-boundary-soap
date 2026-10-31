"""7. Output files contain the same number of atoms."""
import pytest

from generator.builder import build_bicrystal
from generator.exporters import (write_cif, write_extxyz, write_grain_id_map,
                                 write_lammps_data, write_xyz_with_grain_id)


@pytest.fixture(scope="module")
def exported(tmp_path_factory):
    bc = build_bicrystal(K=3, J=3, L=5)
    d = tmp_path_factory.mktemp("outputs")
    write_lammps_data(bc, d / "m.data")
    write_extxyz(bc, d / "m.xyz")
    write_cif(bc, d / "m.cif")
    write_xyz_with_grain_id(bc, d / "m_grain_id.xyz")
    write_grain_id_map(bc, d / "m.grain_map")
    return bc, d


def _count_lammps_atoms(path):
    n = None
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if "atoms" in line and "atom types" not in line and not line.startswith("#"):
            n = int(line.split()[0])
    return n


def _count_cif_atoms(path):
    return sum(1 for ln in path.read_text().splitlines()
               if ln.startswith("Li"))


def test_all_exporters_same_count(exported):
    bc, d = exported
    n = len(bc.atoms)
    assert _count_lammps_atoms(d / "m.data") == n
    assert int((d / "m.xyz").read_text().splitlines()[0]) == n
    assert _count_cif_atoms(d / "m.cif") == n
    assert int((d / "m_grain_id.xyz").read_text().splitlines()[0]) == n
    map_lines = [ln for ln in (d / "m.grain_map").read_text().splitlines()
                 if ln.strip() and not ln.startswith("#")]
    assert len(map_lines) == n
