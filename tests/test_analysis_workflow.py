"""Independent stage, user-decision and artifact-preservation regressions."""
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("analysis_workflow_tests", ROOT / "analysis/run_workflow.py")
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


class AnalysisWorkflowTests(unittest.TestCase):
    def test_regression_tradeoff_needs_an_explicit_formal_scheme_choice(self):
        def comparison_r(*args,**kwargs):
            self.fake_r(*args,**kwargs)
            if kwargs["stage"]=="analysis_regression-review":
                staging=Path(args[1][2])
                (staging/"regression_review_result.json").write_text(json.dumps(dict(status="USER_REVIEW_REQUIRED")))
        with patch.object(analysis,"run_r45",side_effect=comparison_r):
            analysis.run("GSE123",self.root)
            decision=self.root/"choice.json"
            value=dict(confirmed_by="user",user_statement="TEST: choose .3",resolution=.3)
            decision.write_text(json.dumps(value))
            with self.assertRaisesRegex(ValueError,"formal-scheme choice"):
                analysis.run("GSE123",self.root,step="finalize",resolution_decision=decision)
            value.update(regression_scheme="baseline",regression_statement="TEST: keep baseline after comparison review")
            decision.write_text(json.dumps(value))
            self.assertEqual(analysis.run("GSE123",self.root,step="finalize",resolution_decision=decision),"COMPLETE")

    def setUp(self):
        audit_patch = patch.object(analysis, "audit")
        self.audit = audit_patch.start()
        self.addCleanup(audit_patch.stop)
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.area = self.root / "results/GSE123/03_analysis"
        self.raw = self.root / "results/GSE123/02_QC/02_doublet_cellcycle/seurat_QC-2.rds"
        self.raw.parent.mkdir(parents=True)
        self.raw.write_bytes(b"immutable QC2 fixture")

    def tearDown(self):
        self.temporary.cleanup()

    def fake_r(self, script=None, args=(), **kwargs):
        if script and Path(script).name == "run_analysis.R":
            step, source, staging, request = args
            staging = Path(staging)
            status = {"harmony": "SKIPPED_DESIGN", "cluster-review": "NEEDS_RESOLUTION_SELECTION"}.get(step, "COMPLETE")
            (staging / f"{step.replace('-', '_')}_result.json").write_text(json.dumps({"status": status}))
            name = analysis.STEPS[step][1]
            if step == "harmony":
                name = "seurat_unintegrated.rds"
            (staging / name).write_bytes(f"{step} artifact".encode())
            (staging / f"{step}_statistics.csv").write_text("metric,value\ncells,100\n")

    def test_qc3_preferred_without_requiring_download(self):
        qc3 = self.raw.parents[1] / "03_ambient_RNA/seurat_QC-3.rds"
        qc3.parent.mkdir()
        qc3.write_bytes(b"QC3 fixture")
        self.assertEqual(analysis.qc_input(self.root, "GSE123"), qc3.resolve())

    def test_all_stops_at_review_and_preserves_source(self):
        with patch.object(analysis, "run_r45", side_effect=self.fake_r) as calls:
            self.assertEqual(analysis.run("GSE123", self.root), "NEEDS_RESOLUTION_SELECTION")
        self.assertFalse((self.area / "seurat_analysis.rds").exists())
        self.assertEqual(self.raw.read_bytes(), b"immutable QC2 fixture")
        stages = [c.kwargs["stage"] for c in calls.call_args_list]
        self.assertNotIn("analysis_finalize", stages)
        self.assertTrue((self.area / "seurat_cluster_review.rds").is_file())
        self.audit.assert_called_once_with("GSE123", self.root)
        for path in self.area.rglob("*_result.json"):
            if ".workflow" in path.parts: continue
            result=json.loads(path.read_text())
            self.assertTrue(Path(result['output_object']).is_file())
            self.assertTrue(all(Path(output).exists() for output in result['outputs']))
            self.assertFalse(any(Path(output).name=='request.json' for output in result['outputs']))

    def test_resume_validates_existing_rds_without_recomputation(self):
        with patch.object(analysis, "run_r45", side_effect=self.fake_r):
            analysis.run("GSE123", self.root, step="pca")
        with patch.object(analysis, "run_r45", side_effect=self.fake_r) as calls:
            analysis.run("GSE123", self.root, step="pca")
        self.assertEqual([c.kwargs["stage"] for c in calls.call_args_list], ["analysis_validation"])

    def test_completed_configuration_change_is_rejected(self):
        with patch.object(analysis, "run_r45", side_effect=self.fake_r):
            analysis.run("GSE123", self.root, step="pca")
        config = self.root / "config.json"
        config.write_text('{"assay":"RNA"}')
        with patch.object(analysis, "run_r45", side_effect=self.fake_r):
            with self.assertRaisesRegex(ValueError, "different inputs/configuration"):
                analysis.run("GSE123", self.root, step="pca", config=config)

    def test_harmony_configuration_does_not_invalidate_baseline(self):
        cfg = analysis.configuration()
        before = analysis.stage_configuration(cfg, "pca")
        cfg["harmony"]["technical_fields"] = ["batch"]
        self.assertEqual(before, analysis.stage_configuration(cfg, "pca"))

    def test_modified_artifact_is_preserved_and_rejected(self):
        with patch.object(analysis, "run_r45", side_effect=self.fake_r):
            analysis.run("GSE123", self.root, step="pca")
        artifact = self.area / "seurat_pca.rds"
        artifact.write_bytes(b"changed independently")
        with patch.object(analysis, "run_r45", side_effect=self.fake_r):
            with self.assertRaisesRegex(ValueError, "missing or changed"):
                analysis.run("GSE123", self.root, step="pca")
        self.assertEqual(artifact.read_bytes(), b"changed independently")

    def test_crash_does_not_publish_staged_artifact(self):
        def crash(*args, **kwargs):
            self.fake_r(*args, **kwargs)
            if kwargs["stage"] == "analysis_pca":
                raise RuntimeError("native crash despite PASS in stdout")
        with patch.object(analysis, "run_r45", side_effect=crash):
            with self.assertRaisesRegex(RuntimeError, "native crash"):
                analysis.run("GSE123", self.root, step="pca")
        self.assertFalse((self.area / "seurat_pca.rds").exists())
        self.assertTrue(list((self.area / ".workflow/staging").rglob("seurat_pca.rds")))

    def test_finalize_requires_literal_user_decision(self):
        (self.area).mkdir(parents=True)
        (self.area / "seurat_cluster_review.rds").write_bytes(b"review object")
        with self.assertRaisesRegex(ValueError, "requires --resolution-decision"):
            analysis.run("GSE123", self.root, step="finalize")
        decision = self.root / "resolution.json"
        decision.write_text(json.dumps({"resolution": .6, "confirmed_by": "model", "user_statement": ""}))
        with self.assertRaisesRegex(ValueError, "confirmed_by=user"):
            analysis.run("GSE123", self.root, step="finalize", resolution_decision=decision)
        decision.write_text(json.dumps({"resolution": .6, "confirmed_by": "user",
                                        "user_statement": "Synthetic test choice: resolution 0.6"}))
        with patch.object(analysis, "run_r45", side_effect=self.fake_r):
            self.assertEqual(analysis.run("GSE123", self.root, step="finalize", resolution_decision=decision), "COMPLETE")
        self.assertTrue((self.area / "seurat_analysis.rds").is_file())

    def test_interrupted_publication_resumes_without_recompute(self):
        original_publish = analysis.publish
        def interrupt(path):
            raise KeyboardInterrupt("publication interrupted")
        with patch.object(analysis, "run_r45", side_effect=self.fake_r), patch.object(analysis, "publish", side_effect=interrupt):
            with self.assertRaises(KeyboardInterrupt):
                analysis.run("GSE123", self.root, step="pca")
        with patch.object(analysis, "run_r45", side_effect=self.fake_r) as calls:
            analysis.run("GSE123", self.root, step="pca")
        self.assertEqual([c.kwargs["stage"] for c in calls.call_args_list], ["analysis_validation"])
        state=json.loads((self.area/".workflow/state.json").read_text())
        self.assertEqual(state["pca"],"COMPLETE")
        self.assertIsNone(state["error"])

    def test_audit_failure_keeps_completed_analysis(self):
        self.audit.side_effect = RuntimeError("audit failed")
        with patch.object(analysis, "run_r45", side_effect=self.fake_r):
            with self.assertRaisesRegex(RuntimeError, "audit failed"):
                analysis.run("GSE123", self.root)
        self.assertTrue((self.area / "seurat_cluster_review.rds").is_file())
        state = json.loads((self.area / ".workflow/state.json").read_text())
        self.assertEqual(state["workflow_audit"], "FAILED")

    def test_validation_failure_reuses_completed_computation(self):
        def validation_failure(*args, **kwargs):
            self.fake_r(*args, **kwargs)
            if kwargs["stage"] == "analysis_validation":
                raise ValueError("validation interrupted")
        with patch.object(analysis, "run_r45", side_effect=validation_failure):
            with self.assertRaisesRegex(ValueError, "validation interrupted"):
                analysis.run("GSE123", self.root, step="pca")
        with patch.object(analysis, "run_r45", side_effect=self.fake_r) as calls:
            analysis.run("GSE123", self.root, step="pca")
        self.assertEqual([c.kwargs["stage"] for c in calls.call_args_list], ["analysis_validation"])

    def test_all_finalizes_only_with_explicit_decision(self):
        decision = self.root / "resolution.json"
        decision.write_text(json.dumps({"resolution": .6, "confirmed_by": "user",
                                        "user_statement": "TEST FIXTURE: resolution .6"}))
        with patch.object(analysis, "run_r45", side_effect=self.fake_r):
            self.assertEqual(analysis.run("GSE123", self.root, resolution_decision=decision), "COMPLETE")
        self.assertTrue((self.area / "seurat_analysis.rds").is_file())


if __name__ == "__main__":
    unittest.main()
