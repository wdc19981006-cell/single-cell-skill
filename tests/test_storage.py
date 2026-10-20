"""Lifecycle tests use byte fixtures; real Seurat acceptance has a separate runner."""
import importlib.util
import json
import os
import sys
import tempfile
import unittest
import io
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO / "runtime"))
from storage import Storage, RETAINED, CLEANED, INCOMPLETE, atomic, read, sha256, manifest_valid


def module(name, path):
    spec = importlib.util.spec_from_file_location(name,REPO / path)
    result = importlib.util.module_from_spec(spec); spec.loader.exec_module(result)
    return result


qc = module("storage_qc","qc/run_workflow.py")
analysis = module("storage_analysis","analysis/run_workflow.py")
audit = module("storage_audit","workflow_audit/run_audit.py")


class StorageTests(unittest.TestCase):
    def test_qc_cli_passes_storage_retention_run_and_precheck_options(self):
        with patch.object(qc,"run",return_value="PASS") as run,redirect_stdout(io.StringIO()):
            qc.main(["Synthetic","--storage-mode","full","--retain-checkpoint","qc2",
                     "--run-id","test","--precheck","--output-root",str(self.root)])
        options=run.call_args.kwargs
        self.assertEqual(options["storage_mode"],"full")
        self.assertEqual(options["retain_checkpoints"],["qc2"])
        self.assertEqual(options["run_id"],"test")
        self.assertTrue(options["precheck_only"])

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name); self.dataset = self.root / "results/Synthetic"
        self.raw = self.dataset / "01_download/seurat_raw.rds"
        self.raw.parent.mkdir(parents=True); self.raw.write_bytes(b"original raw")
        self.calls = []
        self.addCleanup(patch.stopall)
        patch.object(qc,"run_r45",side_effect=self.fake_r).start()
        patch.object(analysis,"run_r45",side_effect=self.fake_r).start()
        patch.object(analysis,"audit").start()

    def fake_r(self,script=None,args=(),**kwargs):
        self.calls.append(kwargs.get("stage"))
        if script and Path(script).name=="run_qc.R":
            _,_,stage,*_ = args
            directory,filename,_ = qc.STEPS[stage]
            report_root=Path(os.environ.get("SINGLE_CELL_QC_RUN_ROOT",str(self.dataset/"02_QC")))
            path = (Path(os.environ["SINGLE_CELL_QC_CHECKPOINT_ROOT"]) if stage!="qc3" else
                    report_root / directory) / filename
            path.parent.mkdir(parents=True,exist_ok=True); path.write_bytes(stage.encode())
            report = report_root / directory / "qc_report.csv"
            report.parent.mkdir(parents=True,exist_ok=True)
            report.write_text("section,metric,value\nRUN,status,"+dict(qc1="COMPLETE_QC1",qc2="COMPLETE_QC",qc3="COMPLETE_QC3")[stage]+"\n")
        if script and Path(script).name=="run_analysis.R":
            step,source,staging,request = args; staging=Path(staging)
            status = {"harmony":"SKIPPED_DESIGN","cluster-review":"NEEDS_RESOLUTION_SELECTION"}.get(step,"COMPLETE")
            atomic(staging / (step.replace("-","_")+"_result.json"),dict(status=status))
            name = "seurat_unintegrated.rds" if step=="harmony" else analysis.STEPS[step][1]
            (staging/name).write_bytes((step+str(staging)).encode())
            (staging/(step+".csv")).write_text("metric,value\ncells,160\n")
            if step=="regression-review":
                candidate=staging/"count/seurat_sensitivity.rds"; candidate.parent.mkdir()
                candidate.write_bytes(b"candidate"); (candidate.parent/"comparison_config.json").write_text("{}")

    def choice(self):
        path=self.root/"decision.json"
        atomic(path,dict(confirmed_by="user",user_statement="SYNTHETIC TEST ONLY: select .6",resolution=.6))
        return path

    def prepare(self,mode=None,retain=None):
        qc.run("Synthetic",self.root,storage_mode=mode,retain_checkpoints=retain)
        return Storage(self.dataset)

    def test_qc3_cleanup_preserves_final_raw_and_formal_analysis_can_start(self):
        storage=self.prepare()
        self.assertFalse((storage.work/"qc/seurat_QC-1.rds").exists())
        self.assertFalse((storage.work/"qc/seurat_QC-2.rds").exists())
        self.assertTrue((self.dataset/"02_QC/03_ambient_RNA/seurat_QC-3.rds").is_file())
        self.assertEqual(self.raw.read_bytes(),b"original raw")
        self.assertEqual(analysis.run("Synthetic",self.root),"NEEDS_RESOLUTION_SELECTION")
        self.assertTrue((storage.work/"analysis/baseline/seurat_cluster_review.rds").exists())

    def test_pause_resume_finalize_cleanup_and_restart_without_recompute(self):
        storage=self.prepare()
        analysis.run("Synthetic",self.root)
        review=storage.work/"analysis/baseline/seurat_cluster_review.rds"
        original=sha256(review)
        analysis.run("Synthetic",self.root)
        self.assertEqual(sha256(review),original)
        analysis.run("Synthetic",self.root,step="finalize",resolution_decision=self.choice())
        self.assertFalse(review.exists())
        final=self.dataset/"03_analysis/seurat_analysis.rds"; original=sha256(final)
        self.calls.clear()
        self.assertEqual(analysis.run("Synthetic",self.root),"COMPLETE")
        self.assertEqual(qc.run("Synthetic",self.root),"COMPLETE_QC3")
        self.assertFalse(any(c in ("analysis_pca","analysis_finalize","qc1","qc2","qc3") for c in self.calls))
        self.assertEqual(sha256(final),original)
        for path in (storage.work/"analysis/baseline/state").glob("*.json"):
            value=read(path)
            if "outputs" in value: self.assertTrue(manifest_valid(value))
        self.assertFalse(list((storage.work/"analysis/baseline/state/staging").rglob("*.rds")))
        self.assertFalse(list((storage.work/"analysis/baseline").rglob("seurat_sensitivity.rds")))
        self.assertTrue(list((self.dataset/"03_analysis").rglob("comparison_config.json")))
        _,warnings,_=audit.collect(self.root,"Synthetic")
        self.assertFalse([w for w in warnings if w["category"]=="FILE_STATE_INCONSISTENCY"])
        summary=read(self.dataset/"storage_summary.json")
        self.assertGreater(summary["saved_bytes"],0)

    def test_validation_failure_keeps_review_and_native_crash_is_incomplete(self):
        self.prepare(); analysis.run("Synthetic",self.root)
        def fail(script=None,args=(),**kwargs):
            self.fake_r(script,args,**kwargs)
            if kwargs.get("stage")=="storage_validation": raise RuntimeError("readback failed")
        with patch.object(analysis,"run_r45",side_effect=fail):
            with self.assertRaisesRegex(RuntimeError,"readback failed"):
                analysis.run("Synthetic",self.root,step="finalize",resolution_decision=self.choice())
        self.assertTrue((self.dataset/".work/checkpoints/analysis/baseline/seurat_cluster_review.rds").exists())
        failed=read(self.dataset/".work/checkpoints/analysis/baseline/storage_validation_result.json")
        self.assertEqual(failed["status"],INCOMPLETE)
        self.assertFalse(failed["cleanup_allowed"])
        self.assertEqual(analysis.run("Synthetic",self.root,step="finalize",resolution_decision=self.choice()),"COMPLETE")

    def test_native_analysis_crash_keeps_incomplete_staging_and_publishes_nothing(self):
        self.prepare()
        def crash(script=None,args=(),**kwargs):
            self.fake_r(script,args,**kwargs)
            if kwargs.get("stage")=="analysis_pca": raise RuntimeError("native crash after PASS")
        with patch.object(analysis,"run_r45",side_effect=crash):
            with self.assertRaisesRegex(RuntimeError,"native crash"):
                analysis.run("Synthetic",self.root,step="pca")
        storage=Storage(self.dataset)
        files=[v for v in storage.registry["files"].values() if v["scope"]=="analysis/baseline"]
        self.assertTrue(files)
        self.assertTrue(all(v["lifecycle"]==INCOMPLETE and Path(v["path"]).is_file() for v in files))
        self.assertFalse((storage.work/"analysis/baseline/seurat_pca.rds").exists())

    def test_qc3_validation_failure_preserves_qc2_and_final_save_failure_preserves_review(self):
        def fail_qc(script=None,args=(),**kwargs):
            self.fake_r(script,args,**kwargs)
            if kwargs.get("stage")=="storage_validation": raise RuntimeError("QC final unreadable")
        with patch.object(qc,"run_r45",side_effect=fail_qc):
            with self.assertRaisesRegex(RuntimeError,"unreadable"): qc.run("Synthetic",self.root)
        self.assertTrue((self.dataset/".work/checkpoints/qc/seurat_QC-2.rds").exists())
        qc.run("Synthetic",self.root)
        analysis.run("Synthetic",self.root)
        def fail_final(script=None,args=(),**kwargs):
            self.fake_r(script,args,**kwargs)
            if kwargs.get("stage")=="analysis_finalize": raise RuntimeError("final save failed")
        with patch.object(analysis,"run_r45",side_effect=fail_final):
            with self.assertRaisesRegex(RuntimeError,"save failed"):
                analysis.run("Synthetic",self.root,step="finalize",resolution_decision=self.choice())
        self.assertTrue((self.dataset/".work/checkpoints/analysis/baseline/seurat_cluster_review.rds").exists())
        self.assertFalse((self.dataset/"03_analysis/seurat_analysis.rds").exists())

    def test_different_config_or_changed_final_cannot_reuse_cleaned_completion(self):
        self.prepare(); analysis.run("Synthetic",self.root,resolution_decision=self.choice())
        config=self.root/"config.json"; atomic(config,dict(nfeatures=111))
        with self.assertRaisesRegex(ValueError,"different inputs/configuration"):
            analysis.run("Synthetic",self.root,config=config)
        final=self.dataset/"03_analysis/seurat_analysis.rds"; final.write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError,"changed after verification"):
            analysis.run("Synthetic",self.root)

    def test_cleanup_rejects_formal_successor_from_another_run(self):
        storage=self.prepare("full")
        storage=Storage(self.dataset,"compact")
        with self.assertRaisesRegex(ValueError,"different run"):
            storage.clean("analysis/baseline",storage.proof("qc"))

    def test_full_retains_all_checkpoints_and_selected_retention_works(self):
        storage=self.prepare("full")
        analysis.run("Synthetic",self.root,resolution_decision=self.choice())
        self.assertTrue((storage.work/"qc/seurat_QC-1.rds").exists())
        self.assertTrue((storage.work/"analysis/baseline/count/seurat_sensitivity.rds").exists())
        analysis.run("Synthetic",self.root,storage_mode="compact",retain_checkpoints=["pca","qc2"])
        storage=Storage(self.dataset)
        self.assertTrue((storage.work/"analysis/baseline/seurat_pca.rds").exists())
        self.assertFalse((storage.work/"analysis/baseline/seurat_cluster_review.rds").exists())
        qc.run("Synthetic",self.root)
        self.assertTrue((storage.work/"qc/seurat_QC-2.rds").exists())
        self.assertFalse((storage.work/"qc/seurat_QC-1.rds").exists())

    def test_cleaned_explicit_stage_rebuilds_in_isolated_run(self):
        self.prepare(); analysis.run("Synthetic",self.root,resolution_decision=self.choice())
        final=self.dataset/"03_analysis/seurat_analysis.rds"; original=sha256(final)
        self.calls.clear()
        analysis.run("Synthetic",self.root,step="harmony")
        rebuilt=read(self.dataset/"03_analysis/last_rebuild_result.json")
        self.assertIn("rebuild-",rebuilt["run_id"])
        self.assertIn("analysis_pca",self.calls)
        self.assertIn("analysis_harmony",self.calls)
        self.assertEqual(sha256(final),original)
        self.assertTrue((self.dataset/".work/checkpoints/analysis/runs"/rebuilt["run_id"]/"seurat_pca.rds").exists())

    def test_run_named_baseline_cannot_reuse_default_baseline_checkpoints(self):
        self.prepare(); analysis.run("Synthetic",self.root,step="pca")
        self.calls.clear()
        analysis.run("Synthetic",self.root,step="pca",run_id="baseline")
        self.assertIn("analysis_pca",self.calls)
        self.assertTrue((self.dataset/".work/checkpoints/analysis/runs/baseline/seurat_pca.rds").is_file())
        self.assertTrue((self.dataset/".work/checkpoints/analysis/baseline/seurat_pca.rds").is_file())

    def test_rebuild_uses_nearest_retained_upstream_in_same_origin_run(self):
        storage=self.prepare(retain=["pca","qc1"])
        analysis.run("Synthetic",self.root,resolution_decision=self.choice())
        self.calls.clear(); analysis.run("Synthetic",self.root,step="harmony")
        self.assertNotIn("analysis_pca",self.calls); self.assertIn("analysis_harmony",self.calls)
        self.calls.clear(); qc.run("Synthetic",self.root,step="qc2")
        self.assertNotIn("qc1",self.calls); self.assertIn("qc2",self.calls)

    def test_cleaned_qc_step_rebuild_preserves_final_and_old_reports(self):
        storage=self.prepare()
        final=self.dataset/"02_QC/03_ambient_RNA/seurat_QC-3.rds"; original=sha256(final)
        report=self.dataset/"02_QC/01_basic_QC/qc1_result.json"; before=sha256(report)
        self.assertEqual(qc.run("Synthetic",self.root,step="qc2"),"COMPLETE_QC")
        rebuilt=read(self.dataset/"02_QC/last_rebuild_result.json")
        self.assertEqual(sha256(final),original); self.assertEqual(sha256(report),before)
        checkpoint=storage.work/"qc/rebuilds"/rebuilt["run_id"]/"seurat_QC-2.rds"
        self.assertTrue(checkpoint.exists())
        self.assertEqual(qc.run("Synthetic",self.root,step="qc2",run_id=rebuilt["run_id"]),"COMPLETE_QC")

    def test_old_results_are_never_enrolled_or_cleaned_even_with_mode_option(self):
        old=self.dataset/"03_analysis/seurat_pca.rds"; old.parent.mkdir(parents=True)
        old.write_bytes(b"legacy result")
        storage=Storage(self.dataset,"compact")
        self.assertFalse(storage.active)
        storage.clean("analysis/baseline",{})
        report=storage.report(self.root/"inventory")
        self.assertEqual(old.read_bytes(),b"legacy result")
        self.assertEqual(report["mode"],"protected_legacy")
        self.assertFalse((self.dataset/".work/storage_policy.json").exists())

    def test_user_raw_reference_survives_and_detects_changes(self):
        source=self.root/"user_original.rds"; self.raw.replace(source)
        qc.run("Synthetic",self.root,step="qc1",input_path=source)
        qc.run("Synthetic",self.root)
        self.assertFalse(self.raw.exists()); self.assertTrue(source.exists())
        source.write_bytes(b"changed externally")
        with self.assertRaisesRegex(ValueError,"Referenced input"):
            qc.run("Synthetic",self.root)

    def test_cleanup_write_ahead_recovers_after_interrupted_unlink(self):
        storage=self.prepare("full")
        storage=Storage(self.dataset,"compact")
        original=Path.unlink
        target=storage.work/"qc/seurat_QC-1.rds"
        def interrupt(path,*a,**k):
            if path==target: raise KeyboardInterrupt("unlink interrupted")
            return original(path,*a,**k)
        with patch.object(Path,"unlink",new=interrupt):
            with self.assertRaises(KeyboardInterrupt): storage.clean("qc",storage.proof("qc"))
        record=read(storage.work/"qc/state/qc1.json")
        self.assertTrue(manifest_valid(record)); self.assertEqual(record["lifecycle"],CLEANED)
        storage=Storage(self.dataset); storage.clean("qc",storage.proof("qc"))
        self.assertFalse(target.exists())

    def test_no_cleanup_with_bad_successor_or_incomplete_checkpoint(self):
        storage=self.prepare("full"); storage=Storage(self.dataset,"compact")
        incomplete=storage.work/"qc/crashed.rds"; incomplete.write_bytes(b"crash")
        storage.incomplete(incomplete.parent,"crash","qc")
        storage.clean("qc",storage.proof("qc")); self.assertTrue(incomplete.exists())
        self.assertEqual(storage.registry["files"][str(incomplete)]["lifecycle"],INCOMPLETE)
        final=Path(storage.proof("qc")["path"]); final.write_bytes(b"corrupt")
        with self.assertRaisesRegex(ValueError,"validated"): storage.clean("qc",storage.proof("qc"))


if __name__=="__main__": unittest.main()
