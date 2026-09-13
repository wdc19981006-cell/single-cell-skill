"""Offline regression tests. No real GEO downloads or biological group inference."""
import copy
import io
import json
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / '.agents/skills/geo-single-cell-loader/scripts'
sys.path.insert(0, str(SCRIPTS))
from common import detect_trio, file_type, choose_filtered, local, read_csv, validate_manifest, verify_confirmation, write_csv
from build_manifest import confirm
from download_processed import extract_member
from inspect_geo import parse_soft, propose_input

def row(sample='GSM2'):
    return dict(database='GSE999999999', sample=sample,tissue='Synthetic',disease='Synthetic',source_type='Tissue',group='Fixture',local_path='datasets/GSE999999999/'+sample,file_type='10x_mtx',count_source='counts',count_evidence='synthetic fixture',metadata_evidence='synthetic fixture',files_json='[]')

class LoaderTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.root = Path(self.tmp.name)
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
        r=row(); r.update(local_path='datasets/GSE999999999/x.h5ad',file_type='10x_h5')
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
            with self.assertRaises(ValueError): local(self.root,path,'datasets')
    def test_shared_matrix_requires_cell_map(self):
        a,b=row('GSM2'),row('GSM1'); b['local_path']=a['local_path']
        with self.assertRaisesRegex(ValueError,'Shared'): validate_manifest([a,b],self.root)
    def test_confirmed_cell_map_change_stops(self):
        dataset=self.root/'datasets/GSE999999999'; dataset.mkdir(parents=True)
        mapping=dataset/'cell_map.csv'; mapping.write_text('cell,sample\nc1,GSM2\n',encoding='utf-8')
        r=row(); r['cell_map_path']='datasets/GSE999999999/cell_map.csv'
        report=self.root/'report.csv'; approval=self.root/'approval.json'; output=self.root/'result.csv'
        write_csv(report,[r]); approval.write_text(json.dumps(dict(confirmed_by='user',user_statement='SIMULATED',groups={'GSM2':'Fixture'})))
        confirm(report,approval,output,self.root); validate_manifest(read_csv(output),self.root)
        mapping.write_text('cell,sample\nc2,GSM2\n',encoding='utf-8')
        with self.assertRaisesRegex(ValueError,'changed'): validate_manifest(read_csv(output),self.root)
    def test_confirmation_integrity_and_order(self):
        rows=[row('GSM2'),row('GSM1')]
        report=self.root/'report.csv'; approval=self.root/'approval.json'; output=self.root/'result.csv'
        write_csv(report,rows)
        approval.write_text(json.dumps(dict(confirmed_by='user',user_statement='SIMULATED TEST: GSM1=A; GSM2=B',groups={'GSM1':'A','GSM2':'B'})))
        result=confirm(report,approval,output,self.root)
        self.assertEqual([r['group'] for r in result],['B','A']); verify_confirmation(output)
        output.write_text(output.read_text()+'\n')
        with self.assertRaisesRegex(ValueError,'changed'): verify_confirmation(output)
    def test_confirmation_exact_coverage(self):
        report=self.root/'report.csv'; approval=self.root/'approval.json'
        write_csv(report,[row()]); approval.write_text(json.dumps(dict(confirmed_by='user',user_statement='SIMULATED',groups={})))
        with self.assertRaisesRegex(ValueError,'exactly'): confirm(report,approval,self.root/'result.csv',self.root)
    def test_safe_selective_archive(self):
        path=self.root/'fixture.tar'
        with tarfile.open(path,'w') as t:
            info=tarfile.TarInfo('sample/matrix.mtx'); data=b'counts'; info.size=len(data); t.addfile(info,io.BytesIO(data))
        dest=self.root/'out.mtx'; extract_member(path,'sample/matrix.mtx',dest); self.assertEqual(dest.read_bytes(),b'counts')
        with self.assertRaises(ValueError): extract_member(path,'../outside',dest)
    def test_soft_preserves_facts(self):
        records=parse_soft('^SAMPLE = GSM1\n!Sample_characteristics_ch1 = disease: PDAC\n!Sample_characteristics_ch1 = tissue: uninvolved pancreas')
        self.assertEqual(len(records[0]['fields']['Sample_characteristics_ch1']),2)

if __name__=='__main__': unittest.main(verbosity=2)
