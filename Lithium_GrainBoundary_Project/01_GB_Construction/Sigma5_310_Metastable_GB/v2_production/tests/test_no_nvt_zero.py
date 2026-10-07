"""9. NVT with target 0 K is forbidden in the workflow scripts."""
import re
from pathlib import Path

import pytest

LAMMPS_DIR = Path(__file__).resolve().parents[1] / "lammps"
ALL_IN = sorted(LAMMPS_DIR.glob("in.*"))


def test_quench_targets_1K_not_0K():
    quench = (LAMMPS_DIR / "in.3_quench").read_text()
    assert "fix 1 all nvt temp 300.0 1.0 0.1" in quench
    assert "fix 1 all nvt temp 300.0 0.0" not in quench


@pytest.mark.parametrize("script", [p.name for p in ALL_IN])
def test_no_nvt_temp_target_zero(script):
    text = (LAMMPS_DIR / script).read_text()
    # any "nvt temp A B" where B == 0.0 is forbidden
    for m in re.finditer(r"nvt\s+temp\s+([0-9.]+)\s+([0-9.]+)", text):
        assert float(m.group(2)) > 0.0, f"{script}: NVT target 0 K forbidden"


def test_anneal_targets_300K():
    anneal = (LAMMPS_DIR / "in.2_anneal").read_text()
    assert "fix 1 all nvt temp 300.0 300.0 0.1" in anneal
