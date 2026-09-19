import csv
import json
import sys
import tempfile
import unittest
import subprocess
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / ".agents/skills/geo-single-cell-loader/scripts"))
from audit_run import AUDIT_FILES, finalize, prepare_bundle, restore_provenance, run_id
from common import write_csv


class AuditRunTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.workflow = self.root / "data/GSE123/.workflow"
        self.workflow.mkdir(parents=True)
        (self.workflow / "inspection.json").write_text(json.dumps({
            "gse": "GSE123", "files": [
                {"url": "https://example.org/counts.mtx", "scope": "GSM1", "format": "10x_mtx"},
                {"url": "https://example.org/normalized.csv", "scope": "GSM1", "format": "text"},
            ], "discovery_method": "geo_mcp", "mcp_tools_used": ["get_geo_info", "list_geo_files"],
        }), encoding="utf-8")
        write_csv(self.workflow / "sample_manifest.csv", [{
            "sample": "GSM1", "local_path": "data/GSE123/raw/counts.mtx", "file_type": "10x_mtx",
            "count_source": "counts", "count_evidence": "GEO counts table",
            "metadata_evidence": "GSM1", "files_json": json.dumps([{
                "url": "https://example.org/counts.mtx", "local_path": "data/GSE123/raw/counts.mtx"
            }]), "cell_map_path": "",
        }])
        write_csv(self.workflow / "file_selection_reasons.csv", [{
            "url": "https://example.org/normalized.csv", "reason": "normalized values, not raw counts"
        }])
        write_csv(self.workflow / "input_routing.csv", [{
            "local_path": "data/GSE123/raw/counts.csv", "file_type": "text", "reader": "data.table::fread",
            "matrix_rows": "2", "matrix_columns": "3", "estimated_dense_bytes": "48",
            "physical_ram_bytes": "32000000000", "reader_selection_reason": "dense estimate below 25% RAM",
        }])
        (self.workflow / "group_confirmation.json").write_text('{"confirmed_by":"user"}', encoding="utf-8")

    def tearDown(self):
        self.temporary.cleanup()

    def test_bundle_has_only_allowed_files_and_observed_exclusions(self):
        target = self.root / "bundle"
        prepare_bundle("GSE123", self.root, "success", "", target)
        self.assertEqual({p.name for p in target.iterdir()}, set(AUDIT_FILES))
        trace = (target / "decision_trace.md").read_text(encoding="utf-8")
        self.assertIn("normalized values, not raw counts", trace)
        self.assertIn("get_geo_info", trace)
        self.assertIn("dense estimate below 25% RAM", trace)
        self.assertIn("estimated_dense_bytes=48", trace)
        self.assertNotIn("chain-of-thought", trace)

    def test_new_run_refuses_unmigrated_historical_workflow(self):
        with self.assertRaisesRegex(ValueError, "Preexisting workflow"):
            run_id(self.workflow, require_clean=True)
        self.assertFalse((self.workflow / "audit_run.json").exists())

    def test_failed_push_keeps_workflow_and_essentials(self):
        checkout = self.root / "checkout"
        checkout.mkdir()
        class Temporary:
            def cleanup(self): pass
        with patch("audit_run.clone_audit", return_value=(Temporary(), checkout)), patch(
            "audit_run.git", side_effect=RuntimeError("push denied")
        ):
            with self.assertRaisesRegex(RuntimeError, "push denied"):
                finalize("GSE123", self.root, "failure", "download stopped")
        self.assertTrue(self.workflow.exists())
        self.assertTrue((self.root / "data/GSE123/group_confirmation.json").exists())

    def test_successful_push_cleans_only_temporary_workflow(self):
        raw = self.root / "data/GSE123/raw"
        raw.mkdir()
        (raw / "counts.mtx").write_bytes(b"source counts")
        final = self.root / "data/GSE123/seurat_raw.rds"
        final.write_bytes(b"fixture rds")
        (self.workflow / "validation.json").write_text('{"status":"success"}', encoding="utf-8")
        checkout = self.root / "checkout"
        checkout.mkdir()
        class Temporary:
            def cleanup(self): pass
        with patch("audit_run.clone_audit", return_value=(Temporary(), checkout)), patch(
            "audit_run.git", return_value=""
        ):
            finalize("GSE123", self.root, "success", "")
        self.assertFalse(self.workflow.exists())
        self.assertEqual((raw / "counts.mtx").read_bytes(), b"source counts")
        self.assertEqual(final.read_bytes(), b"fixture rds")
        self.assertTrue((self.root / "data/GSE123/group_confirmation.json").exists())
        audit = next((checkout / "GSE123").iterdir())
        self.assertEqual({p.name for p in audit.iterdir()}, set(AUDIT_FILES))

    def test_restore_provenance_reads_prior_audit(self):
        checkout = self.root / "checkout"
        prior = checkout / "GSE123/20260916T000000Z-12345678"
        prior.mkdir(parents=True)
        write_csv(prior / "download_profile.csv", [{
            "local_path": "data/GSE123/raw/counts.mtx", "url": "https://example.org/counts.mtx",
            "bytes": "10", "sha256": "abc", "status": "DOWNLOADED",
        }])
        class Temporary:
            def cleanup(self): pass
        with patch("audit_run.clone_audit", return_value=(Temporary(), checkout)):
            records = json.loads(restore_provenance("GSE123", self.root).read_text())
        self.assertEqual(records[0]["bytes"], 10)
        self.assertEqual(records[0]["sha256"], "abc")

    def test_pushes_to_local_git_remote_then_cleans(self):
        remote = self.root / "audit.git"
        subprocess.run(["git", "init", "--bare", "--initial-branch=main", str(remote)],
                       check=True, capture_output=True)
        with patch("audit_run.AUDIT_REPO", str(remote)):
            finalize("GSE123", self.root, "failure", "fixture interruption")
        self.assertFalse(self.workflow.exists())
        paths = subprocess.run(["git", "--git-dir", str(remote), "ls-tree", "-r", "--name-only", "main"],
                               check=True, capture_output=True, text=True).stdout.splitlines()
        self.assertEqual({Path(p).name for p in paths}, set(AUDIT_FILES))


if __name__ == "__main__":
    unittest.main()
