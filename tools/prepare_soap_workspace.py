"""Prepare an isolated flat workspace for the archived Ni Stage01-06 scripts.

Copies original inputs and scripts without executing them or seeding old results.
Requires only the Python standard library. Existing destinations are refused.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
PROJECT = "Ni_SOAP_Method_Exploration_by_Stage"
STEPS = (
    ("inspect_structures.py", "soap_test.py", "stage1_qc.py"),
    ("stage2_profiles.py", "stage2_soap.py", "stage2_5_audit.py", "stage2_75_full.py"),
    ("stage3a_lae_threshold.py", "stage3b_lae_ler.py"),
    ("stage4a_continuous_disorder.py", "stage4b_disorder_energy.py"),
    ("stage5a_physical_spatial_residual.py", "stage5b_soap_direction.py"),
    ("stage6_final_analysis.py",),
)


def source_files(root: Path) -> dict[str, Path]:
    project = root / PROJECT
    files: dict[str, Path] = {}
    for stage in range(1, 7):
        matches = sorted(project.glob(f"Stage{stage:02d}_*/02_脚本"))
        if len(matches) != 1:
            raise ValueError(f"Stage{stage:02d}: expected one script directory, found {len(matches)}")
        for path in sorted(matches[0].iterdir()):
            if path.is_file() and path.suffix in {".py", ".cfg"}:
                if path.name in files:
                    raise ValueError(f"Duplicate filename: {path.name}")
                files[path.name] = path
    required = {name for stage in STEPS for name in stage}
    required.update({"cfg_reader.py", "ni_stable.cfg", "ni_metastable.cfg"})
    missing = required - files.keys()
    if missing:
        raise ValueError("Missing inputs/scripts: " + ", ".join(sorted(missing)))
    requirements = project / "06_环境与实现/requirements_stages.txt"
    if not requirements.is_file():
        raise ValueError(f"Missing dependency list: {requirements}")
    files["requirements.txt"] = requirements
    return files


def prepare(root: Path, destination: Path) -> Path:
    root = root.resolve()
    destination = destination.resolve()
    if destination.exists():
        raise ValueError(f"Destination already exists: {destination}")
    # Never add runtime copies inside archived scientific source directories.
    if destination.is_relative_to(root) and not destination.is_relative_to(root / "runs"):
        raise ValueError("Inside the repository, choose a destination under runs/")
    if destination == root / "runs":
        raise ValueError("Choose a named run below runs/, such as runs/ni-soap-001")
    sources = source_files(root)  # Validate the entire plan before creating files.
    try:
        commit = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = None
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".soap-prepare-", dir=destination.parent) as temp:
        staging = Path(temp) / "workspace"
        staging.mkdir()
        records = []
        for name, source in sources.items():
            target = staging / name
            shutil.copyfile(source, target)
            records.append({
                "source": source.relative_to(root).as_posix(),
                "destination": name,
                "sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
            })
        (staging / "outputs").mkdir()
        manifest = {
            "source_head": commit,
            "content_basis": "working-tree bytes; per-file SHA256 identifies copied content",
            "files": records,
            "steps": STEPS,
        }
        (staging / "SOURCE_MANIFEST.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8",
        )
        lines = [
            "# Ni SOAP run order", "",
            "Prepared scripts and inputs only. No scientific computation has run.",
            "Run from this directory using a separate Python environment:", "",
            "```bash", "python -m pip install -r requirements.txt", "```", "",
            "Check Stage01 results before proceeding. Each command must succeed.",
            "Outputs are shared under outputs/. Never mix results from other runs.", "",
        ]
        for number, stage in enumerate(STEPS, 1):
            lines.extend([f"## Stage{number:02d}", "", "```bash"])
            lines.extend(f"python {name}" for name in stage)
            lines.extend(["```", ""])
        (staging / "RUN_ORDER.md").write_text("\n".join(lines), encoding="utf-8")
        staging.rename(destination)
    return destination


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dest", type=Path, required=True, help="New run directory")
    args = parser.parse_args()
    try:
        result = prepare(REPO_ROOT, args.dest)
    except (ValueError, OSError) as error:
        parser.exit(1, f"error: {error}\n")
    print(f"Prepared: {result}")
    print("No computation started. Read RUN_ORDER.md before running scripts.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
