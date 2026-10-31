"""Protocol and provenance tests. Fixtures are synthetic and are not MD results."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
CURRENT = ROOT / "grain_boundary/current"
sys.path.insert(0, str(CURRENT))
import checkpoint
import run as workflow
from generator.builder import build_bicrystal
from generator.verify import check_counts, check_periodicity


def make_candidate(base):
    run = base / "candidate"
    stage = run / "relax_aniso"
    stage.mkdir(parents=True)
    checkpoint.write_json(run / "run.json", {"protocol_id": checkpoint.PROTOCOL, "status": "candidate_pending_review"})
    for name in ("candidate.data", "candidate.restart", "trajectory.dump"):
        (stage / name).write_text("Synthetic test fixture; not a scientific result.\n", encoding="utf-8")
    header = "Step Atoms Temp PotEng Press Pxx Pyy Pzz Pxy Pxz Pyz Lx Ly Lz\n"
    rows = [f"{step} 14400 10 -27200 0 0 0 0 0 0 0 66.5 66.5 70.1\n" for step in range(0, 5001, 200)]
    (stage / "log.lammps").write_text(header + "".join(rows) + "CURRENT_RELAX_ANISO_DONE\n", encoding="utf-8")
    for name in ("model.json", "PROTOCOL.md", "lammps/relax_aniso.in"):
        target = run / "source" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((CURRENT / name).read_bytes())
    return run


def approve_fixture(run):
    review_path = run / "review/review.json"
    review = json.loads(review_path.read_text(encoding="utf-8"))
    review.update(approved=True, reviewed_by="synthetic-test-reviewer")
    evidence = run / "review/synthetic_structure_evidence.txt"
    evidence.write_text("Synthetic fixture only; no real structure evaluated.\n", encoding="utf-8")
    for key, value in review["checks"].items():
        value.update(accepted=True, criterion_and_result="Synthetic fixture acceptance for provenance test only.",
                     evidence=["review/synthetic_structure_evidence.txt" if key == "structure" else "review/thermo_summary.json"])
    checkpoint.write_json(review_path, review)
    return review_path


class StandardWorkflowTests(unittest.TestCase):
    def test_current_default_is_14400_with_balanced_periodic_grains(self):
        model = build_bicrystal()
        self.assertEqual(check_counts(model)["final_total"], 14400)
        self.assertEqual(model.meta["final_grain1"], 7200)
        self.assertEqual(model.meta["removed_total"], 480)
        self.assertTrue(check_periodicity(model)["ok"])
        self.assertAlmostEqual(model.a, workflow.MODEL["a_A"])

    def test_current_chain_has_one_aniso_stage_after_quench(self):
        plan = workflow.command_plan(Path("unused"), None, "lmp", 6, "mpirun", 5000)
        self.assertEqual([stage["stage"] for stage in plan], ["relax_initial", "anneal", "quench", "relax_aniso"])
        last = plan[-1]["command"]
        self.assertIn(str(Path("unused/quench/final.restart")), last)
        text = (CURRENT / "lammps/relax_aniso.in").read_text(encoding="utf-8")
        self.assertIn("temp 10.0 10.0 0.1 aniso 0.0 0.0 ${pdamp}", text)
        self.assertIn("pxy pxz pyz", text)
        self.assertNotIn("cna/atom", text)

    def test_from_quench_dry_run_starts_no_process_and_creates_no_outputs(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            restart = base / "quenched.restart"
            restart.write_text("Synthetic test fixture", encoding="utf-8")
            destination = base / "new run"
            result = subprocess.run([
                sys.executable, str(CURRENT / "run.py"), "--out", str(destination),
                "--from-quench", str(restart), "--dry-run",
            ], capture_output=True, text=True, check=True)
            plan = json.loads(result.stdout)
            self.assertFalse(plan["build_initial"])
            self.assertEqual([stage["stage"] for stage in plan["stages"]], ["relax_aniso"])
            self.assertFalse(destination.exists())

    def test_inspection_counts_actual_tail_samples_and_does_not_approve(self):
        with tempfile.TemporaryDirectory() as temp:
            run = make_candidate(Path(temp))
            checkpoint.inspect(run)
            summary = json.loads((run / "review/thermo_summary.json").read_text())
            self.assertEqual(summary["total_samples"], 26)
            self.assertEqual(summary["window_samples"], 13)
            self.assertEqual(summary["window_steps"], [2600.0, 5000.0])
            self.assertAlmostEqual(summary["window_duration_ps"], 2.4)
            review = json.loads((run / "review/review.json").read_text())
            self.assertFalse(review["approved"])
            self.assertTrue(all(not item["accepted"] for item in review["checks"].values()))

    def test_unreviewed_candidate_cannot_be_saved(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            run = make_candidate(base)
            checkpoint.inspect(run)
            with self.assertRaisesRegex(ValueError, "Explicit approval"):
                checkpoint.promote(run, run / "review/review.json", base / "standard")
            self.assertFalse((base / "standard").exists())

    def test_changed_candidate_invalidates_review(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            run = make_candidate(base)
            checkpoint.inspect(run)
            review = approve_fixture(run)
            (run / "relax_aniso/candidate.restart").write_text("Changed fixture", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "changed after inspection"):
                checkpoint.promote(run, review, base / "standard")

    def test_structure_review_requires_evidence(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            run = make_candidate(base)
            checkpoint.inspect(run)
            review_path = approve_fixture(run)
            review = json.loads(review_path.read_text())
            review["checks"]["structure"]["evidence"] = []
            checkpoint.write_json(review_path, review)
            with self.assertRaisesRegex(ValueError, "Incomplete structure review"):
                checkpoint.promote(run, review_path, base / "standard")

    def test_reviewed_save_preserves_state_and_evidence_and_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            run = make_candidate(base)
            checkpoint.inspect(run)
            review = approve_fixture(run)
            destination = base / "standard"
            checkpoint.promote(run, review, destination)
            self.assertEqual((destination / "standard_start.restart").read_bytes(), (run / "relax_aniso/candidate.restart").read_bytes())
            manifest = json.loads((destination / "STANDARD_START.json").read_text())
            for name, checksum in manifest["files"].items():
                self.assertEqual(checkpoint.digest(destination / name), checksum)
            with self.assertRaisesRegex(ValueError, "must be new"):
                checkpoint.promote(run, review, destination)

    def test_evidence_cannot_escape_candidate_directory(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaisesRegex(ValueError, "leaves candidate"):
                checkpoint.inside(Path(temp), "../outside.txt")
            with self.assertRaisesRegex(ValueError, "must be relative"):
                checkpoint.inside(Path(temp), str((Path(temp) / "inside.txt").resolve()))


if __name__ == "__main__":
    unittest.main()
