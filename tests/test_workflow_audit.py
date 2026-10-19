"""Audit distinguishes missing evidence and preserves source provenance."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("audit_runner", ROOT / "workflow_audit/run_audit.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class AuditTests(unittest.TestCase):
    def test_qc_checkpoint_checksum_mismatch_is_failed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp).resolve();area=root/"results/GSE123/02_QC/01_basic_QC"
            area.mkdir(parents=True)
            (area/"sample_QC_summary.csv").write_text("sample,input_cells,retained_cells\nS1,100,95\n")
            artifact=area/"seurat_QC-1.rds";artifact.write_bytes(b"synthetic validated fixture")
            state=area.parent/".workflow/state.json";state.parent.mkdir()
            state.write_text(json.dumps(dict(qc1="COMPLETE_QC1",qc1_sha256=audit.digest(artifact))))
            rows,_,_=audit.collect(root,"GSE123")
            self.assertTrue(any(row['metric']=='checkpoint_validation' and row['status']=='PASS' for row in rows))
            artifact.write_bytes(b"changed")
            _,warnings,_=audit.collect(root,"GSE123")
            self.assertTrue(any(w['category']=='FILE_STATE_INCONSISTENCY' and w['status']=='FAILED' for w in warnings))

    def test_missing_numeric_diagnostic_remains_missing_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp).resolve();area=root/"results/GSE123/03_analysis/01_pca_regression"
            area.mkdir(parents=True)
            (area/"regression_diagnostics.csv").write_text("variable,pc,rho,status\nS.Score,1,,MISSING_METADATA\n")
            _,warnings,_=audit.collect(root,"GSE123")
            self.assertTrue(any(w['status']=='MISSING_EVIDENCE' and 'S.Score' in w['message'] for w in warnings))

    def test_verified_finalize_resolves_only_historical_resolution_gate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp).resolve(); area=root/"results/GSE123/03_analysis"
            state=area/".workflow/state.json"; state.parent.mkdir(parents=True)
            state.write_text(json.dumps({"cluster-review":"NEEDS_RESOLUTION_SELECTION","finalize":"COMPLETE"}))
            source=area/"seurat_cluster_review.rds"; source.write_bytes(b"review fixture")
            final=area/"seurat_analysis.rds"; final.write_bytes(b"validated final fixture")
            manifest=area/".workflow/finalize.json"
            manifest.write_text(json.dumps(dict(status="COMPLETE",binding=dict(input=str(source),input_sha256=audit.digest(source),
                decision=dict(resolution=.3,confirmed_by="user",user_statement="TEST FIXTURE: select .3")),
                outputs={str(final):audit.digest(final)})))
            rows,warnings,_=audit.collect(root,"GSE123")
            self.assertTrue(any(r['metric']=='resolution_selection' and r['status']=='RESOLVED' for r in rows))
            self.assertTrue(any(w['status']=='RESOLVED' and 'NEEDS_RESOLUTION_SELECTION' in w['message'] for w in warnings))
            final.write_bytes(b"modified")
            _,warnings,_=audit.collect(root,"GSE123")
            self.assertTrue(any(w['category']=='FILE_STATE_INCONSISTENCY' and w['status']=='FAILED' for w in warnings))
            self.assertFalse(any(w['status']=='RESOLVED' for w in warnings))

    def test_partial_evidence_is_observed_with_missing_stages_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            area = root / "results/GSE123/02_QC/01_basic_QC"
            area.mkdir(parents=True)
            path = area / "sample_QC_summary.csv"
            path.write_text("sample,input_cells,retained_cells\nS1,100,80\n", encoding="utf-8")
            before = path.read_bytes()
            rows, warnings, inputs = audit.collect(root, "GSE123")
            self.assertEqual(path.read_bytes(), before)
            self.assertEqual([r["value"] for r in rows], ["100", "80"])
            self.assertEqual(inputs[str(path)], audit.digest(path))
            self.assertTrue(any(w["stage"] == "qc3" for w in warnings))
            self.assertFalse(any(r["stage"] == "qc3" for r in rows))

    def test_accession_cannot_escape(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                audit.collect(Path(tmp), "GSE1/../../escape")

    def test_header_only_summary_is_reported_as_empty_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            area = root/'results/GSE123/01_download'
            area.mkdir(parents=True)
            (area/'download_summary.csv').write_text('sample,cells\n',encoding='utf-8')
            rows, warnings, inputs = audit.collect(root,'GSE123')
            self.assertEqual(rows,[])
            self.assertTrue(any(w['category']=='EMPTY_EVIDENCE' and w['stage']=='download' for w in warnings))

    def test_cell_metadata_and_cycle_scores_are_not_expanded_into_audit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            for relative, name in (("01_download", "metadata.csv"),
                                   ("02_QC/02_doublet_cellcycle", "cell_cycle_cells.csv"),
                                   ("03_analysis/03_clustering", "cluster_assignments.csv")):
                area = root / "results/GSE123" / relative
                area.mkdir(parents=True)
                (area/name).write_text("cell,sample,nCount_RNA\nCELL1,S1,100\n",encoding="utf-8")
            rows, warnings, inputs = audit.collect(root,"GSE123")
            self.assertEqual(rows,[])
            self.assertEqual(inputs,{})

    def test_empty_warning_table_means_no_warnings(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            area = root / "results/GSE123/03_analysis/01_pca_regression"
            area.mkdir(parents=True)
            (area / "pca_warnings.csv").write_text("category,severity,message,recommendation\n", encoding="utf-8")
            _, warnings, _ = audit.collect(root, "GSE123")
            self.assertFalse(any(w["category"] == "EMPTY_EVIDENCE" for w in warnings))

    def test_audit_bookkeeping_does_not_create_a_new_scientific_snapshot(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            state = root / "results/GSE123/03_analysis/.workflow/state.json"
            state.parent.mkdir(parents=True)
            state.write_text(json.dumps({"pca":"COMPLETE","updated_at":"first"}))
            def render(_script, args, **kwargs):
                (Path(args[0]) / "workflow_audit.pdf").write_bytes(b"%PDF-test-fixture")
            with patch.object(audit, "run_r45", side_effect=render) as calls:
                first = audit.run("GSE123", root)
                state.write_text(json.dumps({"pca":"COMPLETE","updated_at":"second","workflow_audit":"COMPLETE"}))
                self.assertEqual(audit.run("GSE123", root), first)
                self.assertEqual(calls.call_count, 1)

    def test_real_workflow_state_paths_expose_failures_and_user_gates(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            qc = root / "results/GSE123/02_QC/.workflow/state.json"
            qc.parent.mkdir(parents=True)
            qc.write_text(json.dumps(dict(qc1="COMPLETE_QC1",qc2="FAILED",
                error=dict(stage="qc2",message="Actual captured R error"))),encoding="utf-8")
            analysis = root / "results/GSE123/03_analysis/.workflow/state.json"
            analysis.parent.mkdir(parents=True)
            analysis.write_text(json.dumps({"cluster-review":"NEEDS_RESOLUTION_SELECTION","error":None}),encoding="utf-8")
            rows, warnings, inputs = audit.collect(root,"GSE123")
            self.assertEqual(set(inputs),{str(qc),str(analysis)})
            self.assertTrue(any(row['stage']=='clustering' and row['value']=='NEEDS_RESOLUTION_SELECTION' for row in rows))
            self.assertTrue(any(w['category']=='TECHNICAL_ERROR' and 'Actual captured R error' in w['message'] for w in warnings))
            self.assertTrue(any(w['category']=='BIOLOGICAL_DECISION' and w['stage']=='clustering' for w in warnings))

    def test_renderer_changes_version_report_without_changing_scientific_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            def render(_script, args, **kwargs):
                (Path(args[0]) / "workflow_audit.pdf").write_bytes(b"%PDF-test-fixture")
            original_digest = audit.digest
            def updated_digest(path):
                return "updated-renderer" if Path(path).name == "render_audit.R" else original_digest(path)
            with patch.object(audit, "run_r45", side_effect=render) as calls:
                first = audit.run("GSE123", root)
                with patch.object(audit, "digest", side_effect=updated_digest):
                    second = audit.run("GSE123", root)
                    self.assertNotEqual(first, second)
                    self.assertEqual(audit.run("GSE123", root), second)
                self.assertEqual(calls.call_count, 2)

    def test_changed_evidence_versions_are_reused_and_changed_reports_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            download = root/'results/GSE123/01_download'
            download.mkdir(parents=True)
            summary = download/'download_summary.csv'
            summary.write_text('sample,cells\nS1,10\n',encoding='utf-8')
            def fixture_render(_script,args,**_kwargs):
                # Explicitly synthetic renderer fixture; actual PDF rendering
                # is exercised separately through R 4.5.0 in acceptance tests.
                (Path(args[0])/'workflow_audit.pdf').write_bytes(b'%PDF-synthetic-unit-fixture')
            with patch.object(audit,'run_r45',side_effect=fixture_render) as render:
                original = audit.run('GSE123',root)
                protected = {p.name:p.read_bytes() for p in original.iterdir() if p.is_file()}
                summary.write_text('sample,cells\nS1,11\n',encoding='utf-8')
                version = audit.run('GSE123',root)
                self.assertNotEqual(original,version)
                self.assertEqual(audit.run('GSE123',root),version)
                self.assertEqual(render.call_count,2)
                self.assertEqual(protected,{p.name:p.read_bytes() for p in original.iterdir() if p.is_file()})
                report = version/'workflow_audit.pdf'
                report.write_bytes(report.read_bytes()+b'changed')
                with self.assertRaisesRegex(ValueError,'Saved audit output changed'):
                    audit.run('GSE123',root)


if __name__ == "__main__":
    unittest.main()
