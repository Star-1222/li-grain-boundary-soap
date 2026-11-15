"""Verify input preservation and isolation of the current preparation entry."""

import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("current_soap_prepare", ROOT / "soap/current/prepare.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class CurrentSOAPPreparationTests(unittest.TestCase):
    def test_original_sources_and_verified_inputs_are_copied_without_results(self):
        with tempfile.TemporaryDirectory() as temp:
            destination = module.prepare(ROOT, Path(temp) / "workspace")
            manifest = json.loads((destination / "SOURCE_MANIFEST.json").read_text(encoding="utf-8"))
            for item in manifest["files"]:
                data = (destination / item["destination"]).read_bytes()
                self.assertEqual(data, (ROOT / item["source"]).read_bytes())
                self.assertEqual(hashlib.sha256(data).hexdigest(), item["sha256"])
            self.assertEqual(len(list((destination / "tmp/soap_timing").glob("*.cfg"))), 2)
            self.assertTrue((destination / module.CODE_DIR / "run_sweep.py").is_file())
            self.assertTrue((destination / module.CODE_DIR / "radial_followup.py").is_file())
            self.assertFalse((destination / "output").exists())
            inputs = json.loads((destination / "tmp/soap_timing/inputs.json").read_text(encoding="utf-8"))
            self.assertEqual(len(inputs), 2)

    def test_existing_workspace_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as temp:
            destination = Path(temp)
            sentinel = destination / "result.txt"
            sentinel.write_text("keep", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "already exists"):
                module.prepare(ROOT, destination)
            self.assertEqual(sentinel.read_text(), "keep")

    def test_archived_sources_cannot_be_used_as_output(self):
        with self.assertRaisesRegex(ValueError, "under runs"):
            module.prepare(ROOT, ROOT / "soap/history/new-run")


if __name__ == "__main__":
    unittest.main()
