"""Actual R subprocess integration: primary fixtures, pooled mapping and failure/retry.

Run make_fixtures.py first. Never overwrites the primary final RDS.
"""
import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / '.agents/skills/geo-single-cell-loader/scripts'
sys.path.insert(0, str(SCRIPTS))
from build_manifest import confirm
from common import read_csv, write_csv
from download_processed import run, sha256
from sample_info import refresh

GSE = 'GSE999999999'
REL = f'data/{GSE}/.workflow/sample_manifest.csv'

def main(rscript):
    env = dict(os.environ, GEO_SINGLE_CELL_PYTHON=sys.executable, LC_ALL='C')
    log = ROOT/f'data/{GSE}/.workflow/integration-tests.txt'
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text('', encoding='utf-8')
    def rrun(script, root, *args, success=True):
        result = subprocess.run([rscript, str(script), str(root), *args], env=env, capture_output=True)
        output = (result.stdout + result.stderr).decode('utf-8', errors='replace')
        with log.open('a',encoding='utf-8') as f: f.write(f'\n{script.name}: exit={result.returncode}\n{output}\n')
        if (result.returncode == 0) != success:
            raise AssertionError(f'{script.name}: unexpected exit {result.returncode}\n{output}')
        return output
    def build(root, success=True):
        return rrun(SCRIPTS/'build_seurat.R',root,REL,success=success)
    def validate(root):
        rrun(SCRIPTS/'validate_seurat.R',root,REL,f'data/{GSE}/seurat_raw.rds')
    def approve(root, rows):
        workflow=root/f'data/{GSE}/.workflow'; workflow.mkdir(parents=True,exist_ok=True)
        write_csv(workflow/'sample_report.csv',rows)
        (workflow/'group_confirmation.json').write_text(json.dumps(dict(confirmed_by='user',user_statement='SIMULATED INTEGRATION TEST ONLY',groups={r['sample']:r['group'] for r in rows})),encoding='utf-8')
        confirm(workflow/'sample_report.csv',workflow/'group_confirmation.json',workflow/'sample_manifest.csv',root)
    refresh(ROOT,GSE)
    run(ROOT/REL,ROOT)
    build(ROOT); validate(ROOT)
    rrun(ROOT/'tests/test_final_object.R',ROOT)
    info=(ROOT/f'data/{GSE}/sample_info.txt').read_text(encoding='utf-8')
    assert info.startswith('STATUS: COMPLETE')
    for value in ('Total samples:\n4','Total cells:\n16','Total genes:\n250','GROUP SUMMARY','PROCESSING SUMMARY','raw/','seurat_raw.rds','min.cells = 3','min.features = 200','Counts source:','Warnings:'):
        assert value in info, value
    for name in ('FixtureZ','FixtureA','FixtureH','FixtureT'): assert name in info
    assert 'patient:' not in info
    summary=(ROOT/f'data/{GSE}/.workflow/run_summary.txt').read_text(encoding='utf-8')
    for value in ('Input routing:','Performance:','Manifest rows: 4','Unique input signatures: 4','Unique physical expression inputs: 4','Expression matrices actually read: 4','Cell maps actually read: 0','reader: data.table::fread','file_size_bytes:','total_build_seconds:'):
        assert value in summary, value
    final=ROOT/f'data/{GSE}/seurat_raw.rds'; before=sha256(final)
    assert 'Output exists' in build(ROOT,success=False)
    assert sha256(final)==before
    assert (ROOT/f'data/{GSE}/sample_info.txt').read_text(encoding='utf-8')==info
    print('PASS: primary E2E (4 samples, 16 cells, 250 genes), complete TXT and overwrite refusal',flush=True)

    template=read_csv(ROOT/REL)[3]
    with tempfile.TemporaryDirectory(prefix='geo-pooled-') as tmp:
        root=Path(tmp); raw=root/f'data/{GSE}/raw'; raw.mkdir(parents=True)
        with (raw/'pooled.csv').open('w',newline='') as f:
            writer=csv.writer(f); writer.writerow(['gene']+[f'c{i}' for i in range(8)])
            for i in range(250): writer.writerow([f'Gene{i}']+[1]*8)
        rows=[]
        for sample,group in [('PooledB','B'),('PooledA','A')]:
            r=dict(template,sample=sample,group=group,local_path=f'data/{GSE}/raw/pooled.csv',drop_columns='',cell_map_path=f'data/{GSE}/.workflow/cell_map.csv')
            rows.append(r)
        mapping=root/f'data/{GSE}/.workflow/cell_map.csv'
        write_csv(mapping,[dict(cell=f'c{i}',sample='PooledA' if i<4 else 'PooledB') for i in reversed(range(8))])
        approve(root,rows); build(root); validate(root)
        assert 'Total cells:\n8' in (root/f'data/{GSE}/sample_info.txt').read_text(encoding='utf-8')
        pooled_summary=(root/f'data/{GSE}/.workflow/run_summary.txt').read_text(encoding='utf-8')
        for value in ('Manifest rows: 2','Unique input signatures: 1','Unique physical expression inputs: 1','Expression matrices actually read: 1','Cell maps actually read: 1'):
            assert value in pooled_summary, value
        # Corrupt the subordinate mapping; independent manifest reader must reject it.
        mapping.write_text(mapping.read_text()+'\n',encoding='utf-8')
        message=rrun(SCRIPTS/'validate_seurat.R',root,REL,f'data/{GSE}/seurat_raw.rds',success=False)
        assert 'Cell map changed' in message
    print('PASS: pooled matrix, unsorted cell map, 2 samples/8 cells, mapping checksum rejection',flush=True)

    with tempfile.TemporaryDirectory(prefix='geo-failure-') as tmp:
        root=Path(tmp); raw=root/f'data/{GSE}/raw'; raw.mkdir(parents=True)
        source=raw/'counts.h5ad'; shutil.copy2(ROOT/f'data/{GSE}/raw/counts.h5ad',source)
        r=dict(read_csv(ROOT/REL)[2],local_path=f'data/{GSE}/raw/counts.h5ad',count_source='X')
        approve(root,[r]); before=sha256(source)
        build(root,success=False)
        info_path=root/f'data/{GSE}/sample_info.txt'
        failed=info_path.read_text(encoding='utf-8')
        assert failed.startswith('STATUS: BUILD_FAILED') and 'Failure:' in failed
        assert sha256(source)==before and not (root/f'data/{GSE}/seurat_raw.rds').exists()
        workflow=root/f'data/{GSE}/.workflow'
        (workflow/'sample_manifest.csv').rename(workflow/'failed_manifest.csv')
        (workflow/'sample_manifest.confirmation.json').rename(workflow/'failed_manifest.confirmation.json')
        r['count_source']='counts'; approve(root,[r]); build(root); validate(root)
        complete=info_path.read_text(encoding='utf-8')
        assert complete.startswith('STATUS: COMPLETE') and 'Failure:' not in complete
        assert sha256(source)==before
    print('PASS: normalized H5AD build failure recorded, raw preserved, reconfirmed retry succeeds',flush=True)
    for area in ('datasets','manifests','logs','output'): assert not (ROOT/area).exists()
    print('PASS: no legacy runtime directories; full R subprocess logs in '+str(log),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--rscript',required=True)
    main(parser.parse_args().rscript)
