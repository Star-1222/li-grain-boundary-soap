"""8. Final data file velocities must all be zero."""
import re
from pathlib import Path

import pytest

from generator.builder import build_bicrystal
from generator.exporters import write_lammps_data

LAMMPS_DIR = Path(__file__).resolve().parents[1] / "lammps"


@pytest.fixture(scope="module")
def final_data(tmp_path_factory):
    bc = build_bicrystal()  # full model, zero velocities by construction
    p = tmp_path_factory.mktemp("final") / "final.data"
    write_lammps_data(bc, p)
    return p


def test_exported_data_has_zero_velocities(final_data):
    text = final_data.read_text()
    assert "Velocities" in text
    lines = text.splitlines()
    idx = lines.index("Velocities")
    for line in lines[idx + 2:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) == 4 and parts[0].isdigit():
            assert all(abs(float(v)) < 1e-12 for v in parts[1:]), line


def test_final_script_zeroes_velocities_before_minimize():
    script = (LAMMPS_DIR / "in.4_final").read_text()
    m_zero = script.index("velocity all set 0.0 0.0 0.0")
    m_min = script.index("minimize")
    assert m_zero < m_min
