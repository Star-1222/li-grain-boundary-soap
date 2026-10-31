"""Regression checks for repository organization, without running MD or SOAP."""

import ast
import csv
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "maintenance/tools"))
import build_artifact_catalog as catalog
import prepare_soap_workspace as soap


class RepositoryToolsTests(unittest.TestCase):
    def test_root_contains_two_research_entries_and_one_support_directory(self):
        raw = subprocess.check_output([
            "git", "-C", str(ROOT), "ls-files", "--cached", "--others",
            "--exclude-standard", "-z",
        ])
        paths = {p for p in raw.decode("utf-8").split("\0") if p and (ROOT / p).is_file()}
        directories = {p.split("/", 1)[0] for p in paths if "/" in p}
        self.assertEqual(directories, {"grain_boundary", "soap", "maintenance"})
        self.assertEqual({p for p in paths if "/" not in p}, {"README.md", ".gitignore", ".gitattributes"})

    def test_missing_data_manifest_uses_current_project_paths(self):
        with (ROOT / "maintenance/records/DATA_MANIFEST.csv").open(encoding="utf-8-sig", newline="") as stream:
            records = list(csv.DictReader(stream))
        self.assertTrue(records)
        for record in records:
            path = record["相对路径"].replace("\\", "/")
            self.assertIn(path.split("/", 1)[0], {"grain_boundary", "soap", "maintenance"}, path)
            self.assertRegex(record["SHA256"], r"^[0-9a-f]{64}$")

    def test_every_git_file_has_exactly_one_catalog_group(self):
        raw = subprocess.check_output([
            "git", "-C", str(ROOT), "ls-files", "--cached", "--others",
            "--exclude-standard", "-z",
        ])
        expected = {p for p in raw.decode("utf-8").split("\0") if p and (ROOT / p).is_file()}
        actual = []
        for project, _ in catalog.PROJECTS:
            actual.extend(p.relative_to(ROOT).as_posix() for p in catalog.iter_project_files(project))
        for path in expected:
            self.assertEqual(actual.count(path), 1, path)

    def test_future_experiment_and_docs_are_not_silently_omitted(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            for name in ("new_experiment/log.txt", "docs/guide.md", ".github/workflows/check.yml"):
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("test", encoding="utf-8")
            with patch.object(catalog, "ROOT", root):
                found = {p.relative_to(root).as_posix() for p in catalog.iter_project_files("")}
            self.assertEqual(found, {"new_experiment/log.txt", "docs/guide.md", ".github/workflows/check.yml"})

    def test_workspace_has_exact_inputs_and_resolves_local_imports(self):
        with tempfile.TemporaryDirectory() as temp:
            destination = soap.prepare(ROOT, Path(temp) / "run")
            manifest = json.loads((destination / "SOURCE_MANIFEST.json").read_text(encoding="utf-8"))
            self.assertEqual(list((destination / "outputs").iterdir()), [])
            for record in manifest["files"]:
                content = (destination / record["destination"]).read_bytes()
                self.assertEqual(content, (ROOT / record["source"]).read_bytes())
                self.assertEqual(hashlib.sha256(content).hexdigest(), record["sha256"])
            for script in destination.glob("*.py"):
                tree = ast.parse(script.read_text(encoding="utf-8-sig"))
                for node in ast.walk(tree):
                    if isinstance(node, ast.ImportFrom) and node.module:
                        module = node.module.split(".")[0]
                        if module.startswith("stage") or module == "cfg_reader" or module == "soap_test":
                            self.assertTrue((destination / f"{module}.py").is_file(), f"{script.name}: {module}")
            self.assertEqual(len(manifest["steps"]), 6)

    def test_existing_destination_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as temp:
            destination = Path(temp) / "run"
            destination.mkdir()
            sentinel = destination / "result.txt"
            sentinel.write_text("keep", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "already exists"):
                soap.prepare(ROOT, destination)
            self.assertEqual(sentinel.read_text(), "keep")
            self.assertEqual(list(destination.iterdir()), [sentinel])

    def test_archived_source_directories_cannot_be_runtime_targets(self):
        with self.assertRaisesRegex(ValueError, "under runs"):
            soap.prepare(ROOT, ROOT / soap.PROJECT / "new-run")

    def test_missing_input_fails_before_any_destination_is_created(self):
        with tempfile.TemporaryDirectory() as temp:
            fake_root = Path(temp) / "source"
            fake_root.mkdir()
            destination = Path(temp) / "run"
            with self.assertRaisesRegex(ValueError, "Stage01"):
                soap.prepare(fake_root, destination)
            self.assertFalse(destination.exists())

    def test_new_navigation_links_resolve(self):
        files = [ROOT / "README.md", ROOT / "maintenance/README.md",
                 ROOT / "maintenance/tools/README.md",
                 *sorted((ROOT / "maintenance/docs").glob("*.md"))]
        files += [ROOT / name / "README.md" for name in (
            "grain_boundary", "grain_boundary/experiments/phase1_2026-09",
            "soap", "grain_boundary/archive",
        )]
        for file in files:
            for target in re.findall(r"\[[^\]]*\]\(([^)]+)\)", file.read_text(encoding="utf-8")):
                parsed = urlsplit(target.strip("<>"))
                if not parsed.scheme and parsed.path:
                    self.assertTrue((file.parent / unquote(parsed.path)).exists(), f"{file}: {target}")

    def test_ledger_covers_all_git_files_and_is_stable(self):
        # Redirect both outputs; this test never edits the committed ledger.
        with tempfile.TemporaryDirectory() as temp:
            env = dict(os.environ, RESEARCH_ROOT=str(ROOT), LEDGER_MIRROR="",
                       LEDGER_MD_NAME=str(Path(temp) / "ledger.md"),
                       LEDGER_CSV_NAME=str(Path(temp) / "ledger.csv"))
            command = [sys.executable, str(ROOT / "maintenance/tools/artifact_ledger.py"), "--quiet"]
            subprocess.run(command, env=env, check=True, capture_output=True)
            first = [(Path(temp) / name).read_bytes() for name in ("ledger.md", "ledger.csv")]
            subprocess.run(command, env=env, check=True, capture_output=True)
            second = [(Path(temp) / name).read_bytes() for name in ("ledger.md", "ledger.csv")]
            self.assertEqual(first, second)
            with (Path(temp) / "ledger.csv").open(encoding="utf-8-sig", newline="") as stream:
                paths = {row["相对路径"].replace("\\", "/") for row in csv.DictReader(stream)}
            raw = subprocess.check_output([
                "git", "-C", str(ROOT), "ls-files", "--cached", "--others", "--exclude-standard", "-z",
            ])
            expected = {p for p in raw.decode("utf-8").split("\0") if p and (ROOT / p).is_file()}
            self.assertEqual(paths, expected)
            self.assertTrue(any(p.startswith("grain_boundary/experiments/barostat_2026-10/") for p in paths))


if __name__ == "__main__":
    unittest.main()
