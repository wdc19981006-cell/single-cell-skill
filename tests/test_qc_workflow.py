"""Staged QC routing, recovery and immutable-input orchestration."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("qc_workflow_tests", ROOT / "qc/run_workflow.py")
qc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(qc)


class QCWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.dataset = self.root / "results/GSE123"
        self.raw = self.dataset / "01_download/seurat_raw.rds"
        self.raw.parent.mkdir(parents=True)
        self.raw.write_bytes(b"immutable raw fixture")
        self.area = self.dataset / "02_QC"

    def tearDown(self):
        self.temporary.cleanup()

    def complete(self, script=None, args=(), **kwargs):
        stage = kwargs["stage"].removesuffix("_validation")
        if stage in qc.STEPS:
            directory, filename, _ = qc.STEPS[stage]
            destination = self.area / directory
            destination.mkdir(parents=True, exist_ok=True)
            (destination / filename).write_bytes(b"validated fixture " + stage.encode())
            status = {"qc1": "COMPLETE_QC1", "qc2": "COMPLETE_QC", "qc3": "COMPLETE_QC3"}[stage]
            (destination / "qc_report.csv").write_text("section,metric,value\nRUN,status," + status + "\n")

    def test_default_runs_all_steps_with_applicable_dependency_routes(self):
        with patch.object(qc, "run_r45", side_effect=self.complete) as calls:
            self.assertEqual(qc.run("GSE123", self.root), "COMPLETE_QC3")
        processing = [c.kwargs["stage"] for c in calls.call_args_list if c.args]
        self.assertEqual(processing, ["qc1", "qc2", "qc3"])
        routes = [c.kwargs.get("route") for c in calls.call_args_list if c.kwargs.get("route")]
        self.assertIn("qc1", routes)
        self.assertIn("qc", routes)
        self.assertIn("ambient", routes)
        self.assertEqual(self.raw.read_bytes(), b"immutable raw fixture")

    def test_resume_validates_existing_checkpoint_without_rebuilding(self):
        with patch.object(qc, "run_r45", side_effect=self.complete):
            qc.run("GSE123", self.root, step="qc1")
        with patch.object(qc, "run_r45", side_effect=self.complete) as calls:
            qc.run("GSE123", self.root, step="qc1")
        self.assertIn("qc1_validation", [c.kwargs["stage"] for c in calls.call_args_list])
        self.assertNotIn("qc1", [c.kwargs["stage"] for c in calls.call_args_list])

    def test_successful_validation_clears_current_failure_but_keeps_artifact(self):
        with patch.object(qc, "run_r45", side_effect=self.complete):
            qc.run("GSE123", self.root, step="qc1")
        qc.checkpoint(self.area / ".workflow", qc="QC_ARTIFACT_CREATED_BUT_RUNTIME_FAILED",
                      error={"category":"R_SCRIPT_ERROR"})
        with patch.object(qc, "run_r45", side_effect=self.complete):
            qc.run("GSE123", self.root, step="qc1")
        state = json.loads((self.area / ".workflow/state.json").read_text())
        self.assertEqual(state["qc"], "COMPLETE_QC1")
        self.assertIsNone(state["error"])

    def test_optional_dependency_failure_keeps_finished_qc2(self):
        def fail_ambient(*args, **kwargs):
            if kwargs.get("route") == "ambient":
                raise qc.RRunError("qc3_healthcheck", None, 1, "MISSING_PACKAGE")
            return self.complete(*args, **kwargs)
        with patch.object(qc, "run_r45", side_effect=fail_ambient):
            with self.assertRaises(qc.RRunError):
                qc.run("GSE123", self.root)
        self.assertTrue((self.area / "02_doublet_cellcycle/seurat_QC-2.rds").is_file())
        self.assertFalse((self.area / "03_ambient_RNA/seurat_QC-3.rds").is_file())

    def test_explicit_qc2_input_starts_qc3_without_raw(self):
        self.raw.unlink()
        existing = self.root / "saved_qc2.rds"
        existing.write_bytes(b"existing QC2")
        with patch.object(qc, "run_r45", side_effect=self.complete) as calls:
            self.assertEqual(qc.run("GSE123", self.root, step="qc3", input_path=existing), "COMPLETE_QC3")
        processing = [c.kwargs["stage"] for c in calls.call_args_list if c.args]
        self.assertEqual(processing, ["qc3"])
        self.assertEqual(existing.read_bytes(), b"existing QC2")

    def test_decision_gate_cannot_start_processing_or_invent_consent(self):
        def gate(script=None, args=(), **kwargs):
            if kwargs["stage"] == "qc_precheck":
                self.area.mkdir(parents=True, exist_ok=True)
                (self.area / "qc_report.csv").write_text("section,metric,value\nRUN,status,NEEDS_USER_DECISION\n")
                raise qc.RRunError("qc_precheck", script, 1, "R_SCRIPT_ERROR")
        with patch.object(qc, "run_r45", side_effect=gate):
            with self.assertRaises(qc.RRunError):
                qc.run("GSE123", self.root, precheck_only=True)
        state = json.loads((self.area / ".workflow/state.json").read_text())
        self.assertEqual(state["error"]["category"], "QC_NEEDS_USER_DECISION")
        self.assertFalse((self.area / "qc_decision.json").exists())

    def test_legacy_raw_is_read_without_move(self):
        self.raw.unlink()
        legacy_raw = self.root / "data/GSE123/seurat_raw.rds"
        legacy_raw.parent.mkdir(parents=True)
        legacy_raw.write_bytes(b"legacy raw")
        with patch.object(qc, "run_r45", side_effect=self.complete):
            qc.run("GSE123", self.root, step="qc1")
        self.assertEqual(legacy_raw.read_bytes(), b"legacy raw")
        self.assertTrue((self.area / "01_basic_QC/seurat_QC-1.rds").exists())
        self.assertFalse(self.raw.exists())

    def test_native_crash_remains_failure_even_after_artifact_creation(self):
        def crash(*args, **kwargs):
            self.complete(*args, **kwargs)
            if kwargs["stage"] == "qc1":
                raise qc.RRunError("qc1", args[0], 3221225477, "R_NATIVE_CRASH", "PASS appeared before crash")
        with patch.object(qc, "run_r45", side_effect=crash):
            with self.assertRaises(qc.RRunError):
                qc.run("GSE123", self.root, step="qc1")
        state = json.loads((self.area / ".workflow/state.json").read_text())
        self.assertEqual(state["qc"], "QC_ARTIFACT_CREATED_BUT_RUNTIME_FAILED")
        self.assertTrue((self.area / "01_basic_QC/seurat_QC-1.rds").exists())

    def test_modified_raw_input_is_rejected_even_when_r_returns_success(self):
        def corrupt(*args, **kwargs):
            self.complete(*args, **kwargs)
            if kwargs["stage"] == "qc1":
                self.raw.write_bytes(b"changed raw")
        with patch.object(qc, "run_r45", side_effect=corrupt):
            with self.assertRaisesRegex(RuntimeError, "protected input changed"):
                qc.run("GSE123", self.root, step="qc1")

    def test_existing_legacy_qc_is_only_validated_when_qc2_is_requested(self):
        legacy = self.root / "data/GSE123/qc"
        legacy.mkdir(parents=True)
        final = legacy / "seurat_qc.rds"
        final.write_bytes(b"legacy complete QC")
        (legacy / "qc_report.csv").write_text("section,metric,value\nRUN,status,COMPLETE_QC\n")
        with patch.object(qc, "run_r45") as calls:
            self.assertEqual(qc.run("GSE123", self.root, step="qc2"), "COMPLETE_QC")
        processing = [c.kwargs["stage"] for c in calls.call_args_list if c.args]
        self.assertEqual(processing, ["qc_legacy_validation"])
        self.assertEqual(final.read_bytes(), b"legacy complete QC")
        self.assertFalse((self.area / "01_basic_QC/seurat_QC-1.rds").exists())


if __name__ == "__main__":
    unittest.main()
