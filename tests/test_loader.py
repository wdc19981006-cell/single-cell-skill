"""Offline regression tests. No real GEO downloads or biological group inference."""
import copy
import io
import json
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / '.agents/skills/geo-single-cell-loader/scripts'
sys.path.insert(0, str(SCRIPTS))
from common import detect_trio, file_type, choose_filtered, local, read_csv, validate_manifest, verify_confirmation, write_csv
from build_manifest import confirm
from download_processed import extract_member, run, sha256
from inspect_geo import inspect, parse_soft, propose_input
from sample_info import render

def row(sample='GSM2'):
    return dict(database='GSE999999999', sample=sample,tissue='Synthetic',disease='Synthetic',source_type='Tissue',group='Fixture',local_path='data/GSE999999999/raw/'+sample,file_type='10x_mtx',count_source='counts',count_evidence='synthetic fixture',metadata_evidence='synthetic fixture',files_json='[]')

class LoaderTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.root = Path(self.tmp.name)
        self.workflow = self.root/'data/GSE999999999/.workflow'; self.workflow.mkdir(parents=True)
    def tearDown(self): self.tmp.cleanup()
    def test_duplicate_sample(self):
        with self.assertRaisesRegex(ValueError,'Duplicate'): validate_manifest([row(),row()],self.root)
    def test_missing_group(self):
        r = row(); r.pop('group')
        with self.assertRaisesRegex(ValueError,'group'): validate_manifest([r],self.root)
    def test_optional_evidence(self):
        r=row(); r['patient']='P01'
        with self.assertRaisesRegex(ValueError,'patient'): validate_manifest([r],self.root)
    def test_h5ad_routing(self):
        for name in ['x.h5ad','x.h5ad.gz']: self.assertEqual(file_type(name),'h5ad')
        r=row(); r.update(local_path='data/GSE999999999/raw/x.h5ad',file_type='10x_h5')
        with self.assertRaisesRegex(ValueError,'H5AD'): validate_manifest([r],self.root)
    def test_h5_routing(self): self.assertEqual(file_type('x_filtered_feature_bc_matrix.h5'),'10x_h5_candidate')
    def test_filtered_precedence(self):
        raw='sample_raw_feature_bc_matrix.h5'; filtered='sample_filtered_feature_bc_matrix.h5'
        self.assertEqual(choose_filtered([raw,filtered]),filtered)
        with self.assertRaises(ValueError): choose_filtered([filtered,'another_'+filtered])
    def test_trio_plain_gzip_mixed(self):
        for n in ['matrix.mtx.gz','features.tsv','barcodes.tsv.gz']: (self.root/n).write_bytes(b'fixture')
        self.assertEqual(len(detect_trio(self.root)),3)
        (self.root/'matrix.mtx').write_bytes(b'duplicate')
        with self.assertRaises(ValueError): detect_trio(self.root)
    def test_sample_scoped_input_plan(self):
        urls=['https://example.org/GSM1_'+n for n in ['matrix.mtx.gz','features.tsv.gz','barcodes.tsv.gz']]
        plan=propose_input('GSE1','GSM1',urls)
        self.assertEqual(plan['file_type'],'10x_mtx'); self.assertEqual(len(json.loads(plan['files_json'])),3)
    def test_reject_path_escape(self):
        for path in ['../outside','C:/data/x','datasets/../../x','datasets\\x']:
            with self.assertRaises(ValueError): local(self.root,path,'data')
    def test_shared_matrix_requires_cell_map(self):
        a,b=row('GSM2'),row('GSM1'); b['local_path']=a['local_path']
        with self.assertRaisesRegex(ValueError,'Shared'): validate_manifest([a,b],self.root)
    def test_confirmed_cell_map_change_stops(self):
        dataset=self.workflow
        mapping=dataset/'cell_map.csv'; mapping.write_text('cell,sample\nc1,GSM2\n',encoding='utf-8')
        r=row(); r['cell_map_path']='data/GSE999999999/.workflow/cell_map.csv'
        report=self.workflow/'sample_report.csv'; approval=self.workflow/'group_confirmation.json'; output=self.workflow/'sample_manifest.csv'
        write_csv(report,[r]); approval.write_text(json.dumps(dict(confirmed_by='user',user_statement='SIMULATED',groups={'GSM2':'Fixture'})))
        confirm(report,approval,output,self.root); validate_manifest(read_csv(output),self.root)
        mapping.write_text('cell,sample\nc2,GSM2\n',encoding='utf-8')
        with self.assertRaisesRegex(ValueError,'changed'): validate_manifest(read_csv(output),self.root)
    def test_confirmation_integrity_and_order(self):
        rows=[row('GSM2'),row('GSM1')]
        report=self.workflow/'sample_report.csv'; approval=self.workflow/'group_confirmation.json'; output=self.workflow/'sample_manifest.csv'
        write_csv(report,rows)
        approval.write_text(json.dumps(dict(confirmed_by='user',user_statement='SIMULATED TEST: GSM1=A; GSM2=B',groups={'GSM1':'A','GSM2':'B'})))
        result=confirm(report,approval,output,self.root)
        self.assertEqual([r['group'] for r in result],['B','A']); verify_confirmation(output)
        output.write_text(output.read_text()+'\n')
        with self.assertRaisesRegex(ValueError,'changed'): verify_confirmation(output)
    def test_confirmation_exact_coverage(self):
        report=self.workflow/'sample_report.csv'; approval=self.workflow/'group_confirmation.json'
        write_csv(report,[row()]); approval.write_text(json.dumps(dict(confirmed_by='user',user_statement='SIMULATED',groups={})))
        with self.assertRaisesRegex(ValueError,'exactly'): confirm(report,approval,self.workflow/'sample_manifest.csv',self.root)
    def test_safe_selective_archive(self):
        path=self.root/'fixture.tar'
        with tarfile.open(path,'w') as t:
            info=tarfile.TarInfo('sample/matrix.mtx'); data=b'counts'; info.size=len(data); t.addfile(info,io.BytesIO(data))
        dest=self.root/'out.mtx'; extract_member(path,'sample/matrix.mtx',dest); self.assertEqual(dest.read_bytes(),b'counts')
        with self.assertRaises(ValueError): extract_member(path,'../outside',dest)
    def test_soft_preserves_facts(self):
        records=parse_soft('^SAMPLE = GSM1\n!Sample_characteristics_ch1 = disease: PDAC\n!Sample_characteristics_ch1 = tissue: uninvolved pancreas')
        self.assertEqual(len(records[0]['fields']['Sample_characteristics_ch1']),2)

    def confirmed(self, rows):
        report=self.workflow/'sample_report.csv'; approval=self.workflow/'group_confirmation.json'; output=self.workflow/'sample_manifest.csv'
        write_csv(report,rows)
        approval.write_text(json.dumps(dict(confirmed_by='user',user_statement='SIMULATED TEST ONLY',groups={r['sample']:r['group'] for r in rows})),encoding='utf-8')
        confirm(report,approval,output,self.root)
        return output

    def test_all_manifest_paths_confined_to_dataset(self):
        for target in ['data/GSE999999999/raw/../raw/x','/absolute/x','C:/absolute/x', 'data/GSE1/raw/x','data/GSE999999999/.workflow/x','data/GSE999999999/raw/../../../x']:
            r=row(); r['local_path']=target
            with self.subTest(target=target), self.assertRaises(ValueError): validate_manifest([r],self.root)
        r=row(); r['files_json']=json.dumps([dict(url='https://example.org/x',local_path='data/GSE1/raw/x')])
        with self.assertRaises(ValueError): validate_manifest([r],self.root)

    def test_confirmation_requires_user_statement_and_identity(self):
        for approval in [dict(confirmed_by='assistant',user_statement='SIMULATED'),dict(confirmed_by='user',user_statement='')]:
            write_csv(self.workflow/'sample_report.csv',[row()])
            (self.workflow/'group_confirmation.json').write_text(json.dumps(approval),encoding='utf-8')
            with self.assertRaisesRegex(ValueError,'Explicit user'): confirm(self.workflow/'sample_report.csv',self.workflow/'group_confirmation.json',self.workflow/'sample_manifest.csv',self.root)

    def test_group_report_uses_confirmed_mapping_and_omits_absent_optional(self):
        rows=[row('GSM2'),row('GSM1')]; rows[0].update(group='B',author_sample='合成样本'); rows[1]['group']='A'
        self.confirmed(rows)
        info=(self.workflow.parent/'sample_info.txt').read_text(encoding='utf-8')
        self.assertTrue(info.startswith('STATUS: GROUP_CONFIRMED'))
        self.assertIn('合成样本',info); self.assertIn('GROUP SUMMARY',info)
        self.assertIn('B\n'+'-'*60+'\nGSM2',info)
        self.assertNotIn('patient:',info); self.assertNotIn('cohort:',info)

    def test_full_discovery_writes_only_dataset_metadata_without_downloading(self):
        def soft_mock(acc):
            if acc=='GSE999999999': return '^SERIES = '+acc+'\n!Series_title = Synthetic study\n!Series_sample_id = GSM2\n!Series_sample_id = GSM1', 'https://example.org/'+acc
            return '^SAMPLE = '+acc+'\n!Sample_title = '+acc+'\n!Sample_characteristics_ch1 = disease: synthetic\n!Sample_supplementary_file = https://example.org/'+acc+'_matrix.mtx.gz', 'https://example.org/'+acc
        with patch('inspect_geo.soft',side_effect=soft_mock) as soft, patch('inspect_geo.time.sleep'), patch('download_processed.urllib.request.urlopen',side_effect=AssertionError('No bulk download')):
            result=inspect('GSE999999999',self.root)
        self.assertEqual(soft.call_count,3); self.assertTrue(result['complete']); self.assertEqual(result['inspected_samples'],2)
        self.assertTrue((self.workflow/'inspection.json').is_file()); self.assertTrue((self.workflow/'sample_report.csv').is_file())
        info=(self.workflow.parent/'sample_info.txt').read_text(encoding='utf-8')
        self.assertTrue(info.startswith('STATUS: WAITING_FOR_GROUP_CONFIRMATION'))
        self.assertIn('GSM1',info); self.assertIn('GSM2',info); self.assertNotIn('group:',info)
        self.assertFalse((self.workflow.parent/'raw').exists())
        for area in ('datasets','manifests','logs','output'): self.assertFalse((self.root/area).exists())

    def test_incomplete_discovery_cannot_confirm(self):
        (self.workflow/'inspection.json').write_text(json.dumps(dict(complete=False,samples=[dict(accession='GSM2')],total_samples=2)),encoding='utf-8')
        with self.assertRaisesRegex(ValueError,'Incomplete'): self.confirmed([row()])

    def test_development_report_cannot_confirm_without_full_inspection(self):
        folder=self.workflow/'development'; folder.mkdir()
        write_csv(folder/'sample_report.csv',[row()])
        (self.workflow/'group_confirmation.json').write_text(json.dumps(dict(confirmed_by='user',user_statement='SIMULATED',groups={'GSM2':'Fixture'})),encoding='utf-8')
        with self.assertRaisesRegex(ValueError,'canonical'):
            confirm(folder/'sample_report.csv',self.workflow/'group_confirmation.json',self.workflow/'sample_manifest.csv',self.root)

    def test_evidenced_author_sample_ids_cover_all_gsm(self):
        (self.workflow/'inspection.json').write_text(json.dumps(dict(complete=True,samples=[dict(accession='GSM2')],total_samples=1)),encoding='utf-8')
        r=row('AuthorSample'); r['source_gsm']='GSM2'
        manifest=self.confirmed([r]); self.assertEqual(read_csv(manifest)[0]['sample'],'AuthorSample')

    def test_interrupted_download_reuses_completed_files(self):
        r=row(); targets=['data/GSE999999999/raw/GSM2/'+name for name in ('matrix.mtx.gz','features.tsv.gz')]
        r['files_json']=json.dumps([dict(url='https://example.org/'+str(i),local_path=target) for i,target in enumerate(targets)])
        manifest=self.confirmed([r])
        with patch('download_processed.urllib.request.urlopen',side_effect=[io.BytesIO(b'first'),OSError('interrupted')]):
            with self.assertRaisesRegex(OSError,'interrupted'): run(manifest,self.root)
        self.assertEqual(len(json.loads((self.workflow/'download.json').read_text())),1)
        with patch('download_processed.urllib.request.urlopen',return_value=io.BytesIO(b'second')) as request: run(manifest,self.root)
        self.assertEqual(request.call_count,1); self.assertEqual((self.root/targets[0]).read_bytes(),b'first')

    def test_development_probe_does_not_replace_complete_report(self):
        sentinel=self.workflow/'sample_report.csv'; sentinel.write_text('complete report',encoding='utf-8')
        def soft_mock(acc):
            return (f'^SERIES = {acc}\n!Series_sample_id = GSM2\n!Series_sample_id = GSM1' if acc.startswith('GSE') else f'^SAMPLE = {acc}\n!Sample_title = synthetic'), 'https://example.org/'+acc
        with patch('inspect_geo.soft',side_effect=soft_mock),patch('inspect_geo.time.sleep'):
            result=inspect('GSE999999999',self.root,max_samples=1)
        self.assertFalse(result['complete']); self.assertEqual(sentinel.read_text(),'complete report')
        self.assertTrue((self.workflow/'development/inspection.json').exists())

    def test_download_reuse_and_checksum_refusal(self):
        r=row(); target='data/GSE999999999/raw/GSM2/matrix.mtx.gz'
        r['files_json']=json.dumps([dict(url='https://example.org/matrix.mtx.gz',local_path=target)])
        manifest=self.confirmed([r])
        with patch('download_processed.urllib.request.urlopen',return_value=io.BytesIO(b'fixture')) as request:
            run(manifest,self.root)
        self.assertEqual(request.call_count,1)
        with patch('download_processed.urllib.request.urlopen',side_effect=AssertionError('Must reuse')):
            run(manifest,self.root)
            (self.root/target).write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError,'provenance'): run(manifest,self.root)
        log=json.loads((self.workflow/'download.json').read_text())
        self.assertTrue(all(key in log[0] for key in ('url','sha256','bytes','downloaded_at','local_path')))
        self.assertFalse(list((self.workflow.parent/'raw').rglob('*.json')))
        self.assertEqual((self.root/target).read_bytes(),b'changed')

    def test_unverified_preexisting_file_is_not_overwritten(self):
        r=row(); target='data/GSE999999999/raw/GSM2/matrix.mtx.gz'
        r['files_json']=json.dumps([dict(url='https://example.org/x',local_path=target)])
        manifest=self.confirmed([r]); path=self.root/target; path.parent.mkdir(parents=True); path.write_bytes(b'existing')
        with patch('download_processed.urllib.request.urlopen',side_effect=AssertionError('Must stop')):
            with self.assertRaisesRegex(ValueError,'Unverified'): run(manifest,self.root)
        self.assertEqual(path.read_bytes(),b'existing')

    def test_changed_url_and_expected_checksum_stop(self):
        from download_processed import download
        path=self.workflow.parent/'raw/x'; path.parent.mkdir(); path.write_bytes(b'fixture')
        old=dict(url='https://example.org/x',sha256=sha256(path),bytes=7)
        with self.assertRaisesRegex(ValueError,'provenance'): download('https://example.org/y',path,previous=old)
        with self.assertRaisesRegex(ValueError,'Checksum'): download(old['url'],path,'0'*64,old)

    def test_archive_members_are_verified_and_reused(self):
        archive=io.BytesIO()
        with tarfile.open(fileobj=archive,mode='w') as t:
            info=tarfile.TarInfo('sample/counts'); info.size=7; t.addfile(info,io.BytesIO(b'fixture'))
        r=row(); target='data/GSE999999999/raw/GSM2/counts'
        r['files_json']=json.dumps([dict(url='https://example.org/source.tar',local_path=target,member='sample/counts')])
        manifest=self.confirmed([r])
        with patch('download_processed.urllib.request.urlopen',return_value=io.BytesIO(archive.getvalue())) as request: run(manifest,self.root)
        self.assertEqual(request.call_count,1)
        with patch('download_processed.urllib.request.urlopen',side_effect=AssertionError('Must reuse')):
            run(manifest,self.root)
            (self.root/target).write_bytes(b'tampered')
            with self.assertRaisesRegex(ValueError,'provenance'): run(manifest,self.root)
        self.assertEqual((self.root/target).read_bytes(),b'tampered')

if __name__=='__main__': unittest.main(verbosity=2)
