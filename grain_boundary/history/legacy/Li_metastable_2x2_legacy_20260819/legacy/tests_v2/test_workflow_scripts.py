"""Workflow script regression tests for the verified staged workflow."""
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LAMMPS_DIR = ROOT / "lammps"


def read(name):
    return (LAMMPS_DIR / name).read_text(encoding="utf-8")


def active_lines(name):
    return "\n".join(
        line for line in read(name).splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    )


def test_obsolete_generic_resume_is_retired():
    assert not (LAMMPS_DIR / "in.resume").exists()
    assert (LAMMPS_DIR / "in.2_anneal_resume").is_file()
    assert (LAMMPS_DIR / "in.3_quench_resume").is_file()
    assert (LAMMPS_DIR / "resume.sh").is_file()


def test_stage2_resume_preserves_target_timestep():
    text = active_lines("in.2_anneal_resume")
    assert "read_restart ${restart}" in text
    assert "fix 1 all nvt temp 300.0 300.0 0.1" in text
    assert "run 20000 upto" in text
    assert "dump_modify 1 append yes sort id" in text


def test_stage3_resume_preserves_complete_temperature_ramp():
    text = active_lines("in.3_quench_resume")
    assert "read_restart ${restart}" in text
    assert "fix 1 all nvt temp 300.0 1.0 0.1" in text
    assert "run 50000 upto start 0 stop 50000" in text
    assert "dump_modify 1 append yes sort id" in text


def test_dump2data_uses_valid_template_replacement_workflow():
    text = active_lines("in.dump2data")
    assert "read_data ${template}" in text
    assert "read_dump" in text
    assert "box yes" in text
    assert "replace yes" in text
    assert "trim yes" in text
    assert "add keep" in text
    assert "format native" in text
    assert "create_box" not in text
    assert "format line" not in text


def test_run_all_uses_mpi_and_supplies_stage_inputs():
    text = read("run_all.sh")
    assert 'MPI_NP="${MPI_NP:-12}"' in text
    assert '"$MPIEXEC" -np "$MPI_NP" "$LMP"' in text
    assert '-var input_data "$OUTDIR/Sigma5_310.data"' in text
    assert '-var input_restart "$OUTDIR/stage3_quenched.restart"' in text

    for stage in (
        "stage0_check",
        "stage1_minimize",
        "stage2_anneal",
        "stage3_quench",
        "stage4_final",
    ):
        assert stage in text

    assert '-log "$LOGDIR/$stage.lammps"' in text


def test_resume_launcher_routes_each_verified_checkpoint():
    text = read("resume.sh")
    assert "stage1_checkpoint.data" in text
    assert "stage2_anneal.checkpoint.*" in text
    assert "stage3_quench.checkpoint.*" in text
    assert "stage4_checkpoint.restart" in text
    assert "in.2_anneal_resume" in text
    assert "in.3_quench_resume" in text

    for marker in (
        "stage1_minimize.done",
        "stage2_anneal.done",
        "stage3_quench.done",
        "stage4_final.done",
    ):
        assert marker in text


def test_stage_outputs_have_distinct_names():
    anneal = read("in.2_anneal")
    quench = read("in.3_quench")
    final = read("in.4_final")

    assert "stage2_annealed.data" in anneal
    assert "stage2_annealed.restart" in anneal
    assert "stage2_anneal.dump" in anneal

    assert "stage3_quenched.data" in quench
    assert "stage3_quenched.restart" in quench
    assert "stage3_quench.dump" in quench

    assert "stage4_final.data" in final
    assert "stage4_final.restart" in final
    assert "STAGE4_FINAL_DONE" in final


def test_final_minimization_uses_verified_cell_and_checkpoint_policy():
    text = active_lines("in.4_final")
    assert "read_restart ${input_restart}" in text
    assert "fix 1 all box/relax x 0.0" in text
    assert "box/relax y" not in text
    assert "box/relax z" not in text
    assert "stage4_checkpoint.data" in text
    assert "stage4_checkpoint.restart" in text
    assert "${fmax_saved} < 1.0e-6" in text
    assert "${pxx_saved} < 1.0" in text
    assert "${pxx_saved} > -1.0" in text
