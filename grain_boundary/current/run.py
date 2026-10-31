"""Run the adopted Li workflow, or start its aniso stage from a quench restart.

--dry-run prints the command plan without creating files or starting MD.
Outputs are candidates until checkpoint.py records a reviewed standard start.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

CURRENT = Path(__file__).resolve().parent
REPOSITORY = CURRENT.parents[1]
MODEL = json.loads((CURRENT / "model.json").read_text(encoding="utf-8"))
STAGES = (("relax_initial", "RELAX_INITIAL"), ("anneal", "ANNEAL"),
          ("quench", "QUENCH"), ("relax_aniso", "RELAX_ANISO"))


def validate_destination(destination: Path) -> Path:
    destination = destination.resolve()
    if destination.exists():
        raise ValueError(f"Destination already exists: {destination}")
    if destination.is_relative_to(REPOSITORY) and not destination.is_relative_to(REPOSITORY / "runs"):
        raise ValueError("New outputs inside the repository must be under runs/")
    if destination == REPOSITORY / "runs":
        raise ValueError("Choose a named run below runs/")
    if any(char in str(destination) for char in ('"', "\n", "\r")):
        raise ValueError("Output path cannot contain quotes or newlines")
    return destination


def build(destination: Path) -> None:
    from generator.generate import main as generate
    result = generate([
        "--a", str(MODEL["a_A"]), "--k", str(MODEL["K"]),
        "--j", str(MODEL["J"]), "--l", str(MODEL["L"]),
        "--outdir", str(destination),
    ])
    if result:
        raise RuntimeError("Initial geometry verification failed")


def command_plan(destination: Path, restart: Path | None, lmp: str,
                 ranks: int, mpiexec: str, steps: int) -> list[dict]:
    if ranks < 1 or steps < 400 or steps % 200:
        raise ValueError("ranks must be positive; steps must be >=400 and divisible by 200")
    stages = STAGES[-1:] if restart else STAGES
    source = restart or destination / "build/Sigma5_310_K6_J6_L20.data"
    prefix = [lmp] if ranks == 1 else [mpiexec, "-np", str(ranks), lmp]
    plan = []
    for name, marker in stages:
        out = destination / name
        command = prefix + [
            "-in", str(destination / "source/lammps" / f"{name}.in"),
            "-var", "input", str(source), "-var", "outdir", str(out),
            "-var", "kim_model", MODEL["kim_model"],
            "-var", "seed", str(MODEL["seed"]),
            "-var", "pdamp", str(MODEL["aniso_pdamp_ps"]),
            "-var", "nsteps", str(steps), "-log", str(out / "log.lammps"),
        ]
        plan.append({"stage": name, "command": command, "marker": f"CURRENT_{marker}_DONE"})
        source = out / "final.restart"
    return plan


def execute(destination: Path, restart: Path | None, plan: list[dict]) -> None:
    destination.mkdir(parents=True, exist_ok=False)
    sources = []
    for source in sorted(CURRENT.rglob("*")):
        if not source.is_file() or "__pycache__" in source.parts:
            continue
        relative = source.relative_to(CURRENT)
        target = destination / "source" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        sources.append({"path": relative.as_posix(), "sha256": hashlib.sha256(target.read_bytes()).hexdigest()})
    provenance = {"protocol_id": MODEL["protocol_id"], "status": "running",
                  "model": MODEL, "sources": sources, "commands": plan,
                  "quench_restart": None}
    if restart:
        provenance["quench_restart"] = {"path": str(restart), "sha256": hashlib.sha256(restart.read_bytes()).hexdigest()}
    (destination / "run.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    if not restart:
        build(destination / "build")
    environment = dict(os.environ, OMP_NUM_THREADS="1")
    for item in plan:
        stage = destination / item["stage"]
        stage.mkdir()
        print(f"Running: {item['stage']}", flush=True)
        with (stage / "screen.txt").open("w", encoding="utf-8") as screen:
            subprocess.run(item["command"], cwd=destination, env=environment,
                           stdout=screen, stderr=subprocess.STDOUT, check=True)
        log = (stage / "log.lammps").read_text(encoding="utf-8", errors="replace")
        if item["marker"] not in log:
            raise RuntimeError(f"Missing completion marker: {item['stage']}")
        basename = "candidate" if item["stage"] == "relax_aniso" else "final"
        for suffix in ("data", "restart"):
            output = stage / f"{basename}.{suffix}"
            if not output.is_file() or not output.stat().st_size:
                raise RuntimeError(f"Missing output: {output}")
    provenance["status"] = "candidate_pending_review"
    (destination / "run.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print("Candidate complete. Run checkpoint.py inspect; then review before saving a standard start.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--from-quench", type=Path)
    parser.add_argument("--lmp", default="lmp")
    parser.add_argument("--np", type=int, default=6)
    parser.add_argument("--mpiexec", default="mpirun")
    parser.add_argument("--steps", type=int, default=MODEL["aniso_initial_steps"])
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    try:
        destination = validate_destination(args.out)
        restart = args.from_quench.resolve() if args.from_quench else None
        if restart and (not restart.is_file() or any(c in str(restart) for c in ('"', "\n", "\r"))):
            raise ValueError("Quench restart must be an existing file with a valid path")
        plan = command_plan(destination, restart, args.lmp, args.np, args.mpiexec, args.steps)
        if args.dry_run:
            print(json.dumps({"protocol_id": MODEL["protocol_id"], "build_initial": restart is None,
                              "model": MODEL, "stages": plan}, indent=2))
        else:
            execute(destination, restart, plan)
    except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"error: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
