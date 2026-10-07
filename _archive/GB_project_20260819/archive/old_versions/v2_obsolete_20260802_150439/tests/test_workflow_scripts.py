"""12. Workflow script sanity: in.resume differs from stage scripts."""
from pathlib import Path

LAMMPS_DIR = Path(__file__).resolve().parents[1] / "lammps"


def test_resume_is_not_a_copy_of_main_scripts():
    resume = (LAMMPS_DIR / "in.resume").read_text()
    for name in ("in.1_minimize", "in.2_anneal", "in.3_quench", "in.4_final"):
        other = (LAMMPS_DIR / name).read_text()
        # the resume script must not contain the fixed stage commands
        assert other != resume
    assert "read_restart" in resume
    assert "nsteps" in resume


def test_dump2data_uses_add_yes_not_empty_box():
    text = (LAMMPS_DIR / "in.dump2data").read_text()
    assert "read_dump" in text
    assert "add yes" in text
    assert "box yes" in text
    # the legacy empty-box pattern is gone
    assert "create_box" not in text


def test_each_stage_uses_its_own_log():
    run_all = (LAMMPS_DIR / "run_all.sh").read_text()
    for stage in ("stage0_check", "stage1_minimize", "stage2_anneal",
                  "stage3_quench", "stage4_final"):
        assert f"logs/{stage}.lammps" in run_all


def test_stage_outputs_do_not_overwrite():
    # stage 2/3/4 write distinct data/restart/dump names
    anneal = (LAMMPS_DIR / "in.2_anneal").read_text()
    quench = (LAMMPS_DIR / "in.3_quench").read_text()
    final = (LAMMPS_DIR / "in.4_final").read_text()
    assert "stage2_annealed.data" in anneal and "stage2_annealed.restart" in anneal
    assert "stage3_quenched.data" in quench and "stage3_quenched.restart" in quench
    assert "stage4_final.data" in final and "stage4_final.restart" in final
    assert "stage2_anneal.dump" in anneal
    assert "stage3_quench.dump" in quench
    assert "stage4_final.dump" in final
