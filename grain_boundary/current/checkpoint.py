"""Inspect an aniso candidate and promote it only with an explicit review record."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil
import statistics
import tempfile

PROTOCOL = "li-sigma5-2x2-aniso-10k"
REPOSITORY = Path(__file__).resolve().parents[2]
FIELDS = ("Temp", "Pxx", "Pyy", "Pzz", "Pxy", "Pxz", "Pyz", "Lx", "Ly", "Lz")


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def inside(root: Path, name: str) -> Path:
    if Path(name).is_absolute():
        raise ValueError("Evidence paths must be relative to the candidate directory")
    path = (root / name).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError(f"Evidence path leaves candidate directory: {name}")
    return path


def thermo_summary(log: Path) -> dict:
    text = log.read_text(encoding="utf-8", errors="replace")
    if "CURRENT_RELAX_ANISO_DONE" not in text:
        raise ValueError("Current aniso completion marker is missing")
    rows = []
    header = None
    for line in text.splitlines():
        tokens = line.split()
        if tokens and tokens[0] == "Step" and all(key in tokens for key in FIELDS):
            header = tokens
            continue
        if header and len(tokens) == len(header):
            try:
                values = list(map(float, tokens))
            except ValueError:
                continue
            if not all(math.isfinite(value) for value in values):
                raise ValueError("Non-finite thermo value")
            row = dict(zip(header, values))
            if "Atoms" not in row or row["Atoms"] != 14400:
                raise ValueError("Expected 14400 atoms in the current model")
            if rows and row["Step"] <= rows[-1]["Step"]:
                raise ValueError("Thermo steps must increase; concatenated or resumed logs require separate review")
            rows.append(row)
    if len(rows) < 3:
        raise ValueError("At least three thermo samples are required")
    count = max(2, len(rows) // 2)
    tail = rows[-count:]
    summary = {"protocol_id": PROTOCOL, "status": "pending_review",
               "total_samples": len(rows), "window_samples": count,
               "window_steps": [tail[0]["Step"], tail[-1]["Step"]],
               "window_duration_ps": (tail[-1]["Step"] - tail[0]["Step"]) * 0.001,
               "note": "Standard deviation describes sampled fluctuations, not uncertainty of the mean.",
               "metrics": {}}
    split = max(1, len(tail) // 2)
    for field in FIELDS:
        values = [row[field] for row in tail]
        metrics = {"mean": statistics.mean(values), "sample_std": statistics.stdev(values),
                   "first_half_mean": statistics.mean(values[:split]),
                   "second_half_mean": statistics.mean(values[split:])}
        metrics["block_mean_difference"] = metrics["second_half_mean"] - metrics["first_half_mean"]
        if field in ("Lx", "Ly", "Lz"):
            if metrics["mean"] <= 0:
                raise ValueError("Box lengths must be positive")
            metrics["relative_block_drift"] = metrics["block_mean_difference"] / metrics["mean"]
        summary["metrics"][field] = metrics
    return summary


def inspect(run: Path) -> None:
    run = run.resolve()
    metadata = json.loads((run / "run.json").read_text(encoding="utf-8"))
    if metadata.get("protocol_id") != PROTOCOL or metadata.get("status") != "candidate_pending_review":
        raise ValueError("Run is not a completed candidate for the current protocol")
    review_directory = run / "review"
    if review_directory.exists():
        raise ValueError("review/ already exists; preserve its record before creating a new inspection")
    summary = thermo_summary(run / "relax_aniso/log.lammps")
    names = ["run.json", "relax_aniso/candidate.data", "relax_aniso/candidate.restart",
             "relax_aniso/log.lammps", "relax_aniso/trajectory.dump"]
    for name in names:
        path = inside(run, name)
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f"Missing candidate evidence: {name}")
    # Also bind the protocol and source snapshot used for this run.
    names.extend(path.relative_to(run).as_posix() for path in sorted((run / "source").rglob("*")) if path.is_file())
    required_sources = {"source/model.json", "source/PROTOCOL.md", "source/lammps/relax_aniso.in"}
    if not required_sources.issubset(names):
        raise ValueError("Candidate source/protocol snapshot is incomplete")
    manifest = {"protocol_id": PROTOCOL, "files": {name: digest(inside(run, name)) for name in names}}
    review_directory.mkdir()
    write_json(review_directory / "thermo_summary.json", summary)
    manifest["files"]["review/thermo_summary.json"] = digest(review_directory / "thermo_summary.json")
    write_json(review_directory / "candidate_manifest.json", manifest)
    review = {"protocol_id": PROTOCOL, "approved": False, "reviewed_by": "",
              "candidate_manifest_sha256": digest(review_directory / "candidate_manifest.json"),
              "checks": {key: {"accepted": False, "criterion_and_result": "", "evidence": []}
                         for key in ("stress", "box", "structure")}}
    write_json(review_directory / "review.json", review)


def promote(run: Path, review_file: Path, destination: Path) -> None:
    run, destination = run.resolve(), destination.resolve()
    if destination.exists() or destination.is_relative_to(run):
        raise ValueError("Standard-start destination must be new and outside the candidate run")
    if destination.is_relative_to(REPOSITORY) and not destination.is_relative_to(REPOSITORY / "runs"):
        raise ValueError("Saved outputs inside the repository must be under runs/")
    manifest_path = run / "review/candidate_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    review = json.loads(review_file.read_text(encoding="utf-8"))
    if (manifest.get("protocol_id") != PROTOCOL or review.get("protocol_id") != PROTOCOL
            or review.get("approved") is not True or not str(review.get("reviewed_by", "")).strip()
            or review.get("candidate_manifest_sha256") != digest(manifest_path)):
        raise ValueError("Explicit approval for this candidate manifest is required")
    for name, checksum in manifest["files"].items():
        path = inside(run, name)
        if not path.is_file() or digest(path) != checksum:
            raise ValueError(f"Candidate evidence changed after inspection: {name}")
    evidence = set(manifest["files"])
    for key in ("stress", "box", "structure"):
        check = review.get("checks", {}).get(key, {})
        if (check.get("accepted") is not True or not str(check.get("criterion_and_result", "")).strip()
                or not isinstance(check.get("evidence"), list) or not check["evidence"]):
            raise ValueError(f"Incomplete {key} review")
        for name in check["evidence"]:
            if not isinstance(name, str) or not inside(run, name).is_file():
                raise ValueError(f"Missing {key} evidence: {name}")
            evidence.add(name)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".standard-start-", dir=destination.parent) as temp:
        staging = Path(temp) / "bundle"
        staging.mkdir()
        for suffix in ("data", "restart"):
            shutil.copyfile(run / f"relax_aniso/candidate.{suffix}", staging / f"standard_start.{suffix}")
        shutil.copyfile(review_file, staging / "review.json")
        shutil.copyfile(manifest_path, staging / "candidate_manifest.json")
        for name in sorted(evidence):
            target = staging / "evidence" / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(inside(run, name), target)
        checksums = {path.relative_to(staging).as_posix(): digest(path)
                     for path in sorted(staging.rglob("*")) if path.is_file()}
        write_json(staging / "STANDARD_START.json", {"protocol_id": PROTOCOL, "status": "accepted",
                                                    "reviewed_by": review["reviewed_by"], "files": checksums})
        staging.rename(destination)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    inspection = commands.add_parser("inspect")
    inspection.add_argument("--run", type=Path, required=True)
    saving = commands.add_parser("save")
    saving.add_argument("--run", type=Path, required=True)
    saving.add_argument("--review", type=Path, required=True)
    saving.add_argument("--dest", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.action == "inspect":
            inspect(args.run)
            print("Inspection recorded. Review stress, box and structure in review/review.json.")
        else:
            promote(args.run, args.review, args.dest)
            print(f"Saved reviewed standard start: {args.dest.resolve()}")
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.exit(1, f"error: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
