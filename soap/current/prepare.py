"""Prepare one isolated workspace for the recorded Ni SOAP parameter workflow.

Uses standard-library copying only. Does not import DScribe or run calculations.
The archived scientific scripts remain byte-for-byte unchanged.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
CODE_DIR = Path("work/ni_soap_sweep_20261101")


def sources(root: Path) -> dict[str, Path]:
    audit = root / "soap/experiments/2026-11-01_audit"
    files = {str(CODE_DIR / p.name): p for p in sorted((audit / "source").glob("*.py"))}
    files[str(CODE_DIR / "radial_followup.py")] = (
        root / "soap/experiments/2026-11-01_radial_followup/source/radial_followup.py"
    )
    files["tmp/soap_timing/cfg_reader.py"] = (
        root / "soap/experiments/2026-11-01_midterm/source/cfg_reader.py"
    )
    inputs = root / "soap/history/stages/Stage01_数据准备与SOAP初测/02_脚本"
    for name in ("ni_stable.cfg", "ni_metastable.cfg"):
        files[f"tmp/soap_timing/{name}"] = inputs / name
    required_scripts = {"run_sweep.py", "validate_candidates.py", "build_report.py",
                        "extend_radial.py", "metric_diagnostics.py", "verify_results.py"}
    if not required_scripts.issubset({p.name for p in files.values()}):
        raise ValueError("Incomplete archived parameter-workflow source set")
    for path in files.values():
        if not path.is_file():
            raise ValueError(f"Missing source: {path}")
    protocol = json.loads((audit / "protocol.json").read_text(encoding="utf-8"))
    for item in protocol["input_evidence"]:
        source = inputs / item["file"]
        if hashlib.sha256(source.read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError(f"Input differs from recorded audit: {source}")
    return files


def prepare(root: Path, destination: Path) -> Path:
    root, destination = root.resolve(), destination.resolve()
    if destination.exists():
        raise ValueError(f"Destination already exists: {destination}")
    if destination.is_relative_to(root) and not destination.is_relative_to(root / "runs"):
        raise ValueError("Inside the repository, choose a named destination under runs/")
    if destination == root / "runs":
        raise ValueError("Choose a named run below runs/")
    files = sources(root)
    protocol = json.loads((root / "soap/experiments/2026-11-01_audit/protocol.json").read_text(encoding="utf-8"))
    revision = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".soap-prepare-", dir=destination.parent) as temp:
        staging = Path(temp) / "workspace"
        staging.mkdir()
        records, input_records = [], []
        for relative, source in files.items():
            target = staging / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
            digest = hashlib.sha256(target.read_bytes()).hexdigest()
            records.append({"source": source.relative_to(root).as_posix(),
                            "destination": Path(relative).as_posix(), "sha256": digest})
            if source.suffix == ".cfg":
                input_records.append({"file": source.name, "sha256": digest, "bytes": target.stat().st_size})
        (staging / "tmp/soap_timing/inputs.json").write_text(
            json.dumps(input_records, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        (staging / "SOURCE_MANIFEST.json").write_text(json.dumps({
            "source_revision": revision, "files": records,
            "scope": "Prepared original parameter-workflow scripts and verified inputs only; no computation run.",
            "statistics": "20000 stratified pairs; distinct from the PPT 3600-pair sweep",
        }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        (staging / "requirements-recorded.txt").write_text(
            "\n".join(f"{name}=={version}" for name, version in protocol["versions"].items()) + "\n",
            encoding="utf-8")
        (staging / "RUN_ORDER.md").write_text("""# Ni SOAP parameter workflow

This is an isolated copy of the recorded 2026-11-01 audit workflow. No calculation has run.
It uses 20000 stratified pairs and is distinct from the PPT's 3600-pair statistics.
The original scientific scripts are unchanged; the nested layout preserves their imports and output paths.

The recorded Python environment is listed in requirements-recorded.txt. Use an isolated environment.
From work/ni_soap_sweep_20261101, run and review each step before continuing:

```bash
python run_sweep.py prepare
python run_sweep.py sweep
python run_sweep.py convergence
python validate_candidates.py perturb
python validate_candidates.py review
python build_report.py
python verify_results.py
```

Results are written under output/ni_soap_sweep_20261101 inside this workspace.
No existing results are copied. This reproduces an unresolved method, not a validated final parameter choice.
Radial extensions and polynomial diagnostics are optional recorded investigations; consult their reports
before using extend_radial.py, metric_diagnostics.py or radial_followup.py. The tested polynomial l=10
branch failed rotation-invariance acceptance; do not treat it as a validated replacement.
""", encoding="utf-8")
        staging.rename(destination)
    return destination


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dest", type=Path, required=True, help="New isolated workspace")
    args = parser.parse_args()
    try:
        print(f"Prepared: {prepare(REPO_ROOT, args.dest)}")
    except (OSError, ValueError) as error:
        parser.exit(1, f"error: {error}\n")
    print("No computation started. Read RUN_ORDER.md in the destination.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
