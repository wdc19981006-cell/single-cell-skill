"""QC orchestration preserves decisions and never downloads or rebuilds raw."""
import importlib.util
import json
import sys
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
        self.dataset = self.root / "data/GSE123"
        (self.dataset / "qc").mkdir(parents=True)
        (self.dataset / "seurat_raw.rds").write_bytes(b"immutable raw fixture")

    def tearDown(self):
        self.temporary.cleanup()

    def report(self, status):
        (self.dataset / "qc/qc_report.csv").write_text(
            'section,metric,value\nRUN,status,' + status + '\n', encoding='utf-8')

    def test_existing_final_is_validated_without_run_qc(self):
        (self.dataset / "qc/seurat_qc.rds").write_bytes(b"final fixture")
        self.report("COMPLETE_QC")
        with patch.object(qc, "run_r45") as calls:
            self.assertEqual(qc.run("GSE123", self.root), "COMPLETE_QC")
        self.assertEqual([c.kwargs['stage'] for c in calls.call_args_list], ['healthcheck','qc_validation'])
        self.assertEqual((self.dataset / "seurat_raw.rds").read_bytes(), b"immutable raw fixture")

    def test_decision_gate_cannot_start_qc(self):
        def gate(script=None, *args, **kwargs):
            if kwargs['stage']=='qc_precheck':
                self.report('NEEDS_USER_DECISION')
                raise qc.RRunError('qc_precheck',script,1,'R_SCRIPT_ERROR')
        with patch.object(qc, 'run_r45', side_effect=gate) as calls:
            with self.assertRaises(qc.RRunError): qc.run('GSE123', self.root)
        state = json.loads((self.dataset / '.workflow/state.json').read_text())
        self.assertEqual(state['error']['category'], 'QC_NEEDS_USER_DECISION')
        self.assertEqual([c.kwargs['stage'] for c in calls.call_args_list], ['healthcheck','qc_precheck'])
        self.assertFalse((self.dataset / 'qc/qc_decision.json').exists())

    def test_qc_runs_only_after_health_and_precheck(self):
        def complete(script=None, *args, **kwargs):
            if kwargs['stage']=='qc':
                (self.dataset / 'qc/seurat_qc.rds').write_bytes(b'final fixture')
                self.report('COMPLETE_QC')
        with patch.object(qc, 'run_r45', side_effect=complete) as calls:
            self.assertEqual(qc.run('GSE123', self.root), 'COMPLETE_QC')
        self.assertEqual([c.kwargs['stage'] for c in calls.call_args_list], ['healthcheck','qc_precheck','qc'])


if __name__ == '__main__':
    unittest.main()
