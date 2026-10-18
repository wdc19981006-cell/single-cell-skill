"""Canonical download layout and immutable legacy-resume regressions."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / '.agents/skills/geo-single-cell-loader/scripts'))
from common import dataset_paths, dataset_relative, digest, validate_manifest, verify_confirmation, write_csv
from build_manifest import confirm
from audit_run import AUDIT_FILES, finalize
from inspect_geo import propose_input
import run_confirmed as runner

GSE = 'GSE999999999'


class DownloadLayoutTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def test_new_workflow_uses_canonical_stage(self):
        paths = dataset_paths(self.root, GSE)
        self.assertEqual(paths['dataset'], self.root / f'results/{GSE}/01_download')
        urls = [f'https://example.org/{name}' for name in ('matrix.mtx.gz','features.tsv.gz','barcodes.tsv.gz')]
        plan = propose_input(GSE, 'Synthetic1', urls, self.root)
        self.assertEqual(plan['local_path'], f'results/{GSE}/01_download/raw/Synthetic1')

    def test_existing_legacy_paths_and_receipt_bytes_are_preserved(self):
        workflow = self.root / f'data/{GSE}/.workflow'
        workflow.mkdir(parents=True)
        manifest = workflow / 'sample_manifest.csv'
        row = self.row(f'data/{GSE}')
        write_csv(manifest, [row])
        receipt = manifest.with_suffix('.confirmation.json')
        receipt.write_text(json.dumps(dict(manifest_md5=digest(manifest),confirmed_by='user',
            user_statement='SIMULATED FIXTURE ONLY',groups={'Synthetic1':'Fixture'})),encoding='utf-8')
        before = (manifest.read_bytes(), receipt.read_bytes())
        self.assertEqual(dataset_relative(self.root,GSE),f'data/{GSE}')
        validate_manifest([row],self.root)
        verify_confirmation(manifest)
        self.assertEqual(before,(manifest.read_bytes(),receipt.read_bytes()))
        self.assertFalse((self.root/'results').exists())

    def test_downstream_directories_do_not_migrate_legacy_download(self):
        (self.root/f'data/{GSE}').mkdir(parents=True)
        (self.root/f'results/{GSE}/02_QC/01_basic_QC').mkdir(parents=True)
        (self.root/f'results/{GSE}/03_analysis').mkdir(parents=True)
        self.assertEqual(dataset_relative(self.root,GSE),f'data/{GSE}')

    def test_two_download_roots_require_review(self):
        (self.root/f'data/{GSE}').mkdir(parents=True)
        (self.root/f'results/{GSE}/01_download').mkdir(parents=True)
        with self.assertRaisesRegex(ValueError,'Both legacy and canonical'):
            dataset_paths(self.root,GSE)

    def test_canonical_manifest_confines_raw_and_cell_maps(self):
        area = f'results/{GSE}/01_download'
        row = self.row(area)
        validate_manifest([row],self.root)
        for path in (f'data/{GSE}/raw/x',f'results/GSE1/01_download/raw/x',
                     f'{area}/../02_QC/x',f'{area}/raw/../x',f'{area}/.workflow/x'):
            with self.subTest(path=path),self.assertRaises(ValueError):
                validate_manifest([dict(row,local_path=path)],self.root)

    def test_canonical_audit_still_ten_files_and_keeps_user_artifacts(self):
        paths = dataset_paths(self.root,GSE)
        paths['workflow'].mkdir(parents=True)
        (paths['workflow']/'inspection.json').write_text(json.dumps(dict(gse=GSE,files=[])),encoding='utf-8')
        (paths['workflow']/'validation.json').write_text('{"status":"success"}',encoding='utf-8')
        paths['final'].write_bytes(b'synthetic artifact; R validation is tested separately')
        paths['raw'].mkdir()
        for name in ('metadata.csv','download_summary.csv','download_report.pdf','download.log'):
            (paths['dataset']/name).write_bytes(b'required user artifact')
        checkout = self.root/'audit-checkout'
        checkout.mkdir()
        class Temporary:
            def cleanup(self): pass
        with patch('audit_run.clone_audit',return_value=(Temporary(),checkout)),patch('audit_run.git',return_value=''):
            finalize(GSE,self.root,'success','')
        bundles = list((checkout/GSE).iterdir())
        self.assertEqual(len(bundles),1)
        self.assertEqual({p.name for p in bundles[0].iterdir()},set(AUDIT_FILES))
        self.assertFalse(paths['workflow'].exists())
        for name in ('metadata.csv','download_summary.csv','download_report.pdf','download.log','seurat_raw.rds'):
            self.assertTrue((paths['dataset']/name).is_file())

    def test_stage_log_failure_still_finalizes_failure_audit(self):
        paths = dataset_paths(self.root,GSE)
        workflow = paths['workflow']
        workflow.mkdir(parents=True)
        (workflow/'audit_run.json').write_text(json.dumps({'run_id':'20261018T000000Z-12345678'}),encoding='utf-8')
        (workflow/'inspection.json').write_text(json.dumps(dict(gse=GSE,
            title='Synthetic single-cell RNA-seq fixture')),encoding='utf-8')
        paths['final'].write_bytes(b'explicitly synthetic already-validated fixture')
        (paths['dataset']/'download.log').mkdir()
        with patch.object(runner,'execute_r',return_value=0.1),patch.object(runner,'raw_handoff',return_value='Synthetic complete'), \
             patch.object(runner,'finalize',return_value='synthetic-audit-url') as published:
            with self.assertRaisesRegex(RuntimeError,'Could not preserve download.log'):
                runner.run(GSE,self.root)
        self.assertEqual(published.call_args.args[2],'failure')
        self.assertEqual(paths['final'].read_bytes(),b'explicitly synthetic already-validated fixture')
        self.assertEqual(json.loads((workflow/'state.json').read_text())['download_report'],'FAILED')

    @staticmethod
    def row(area):
        return dict(database=GSE,sample='Synthetic1',tissue='Synthetic',disease='Synthetic',
            source_type='Tissue',group='Fixture',local_path=area+'/raw/Synthetic1',
            file_type='10x_mtx',count_source='counts',count_evidence='Synthetic fixture',
            metadata_evidence='Synthetic fixture',files_json='[]')


if __name__ == '__main__':
    unittest.main()
