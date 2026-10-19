"""Execute canonical 10x build/report/recovery through the unchanged R 4.5 runtime.

Uses generated sparse counts only; never initiates a real GEO run or download.
Retains the isolated fixture and JSON evidence under .runtime/ for review.
"""
import gzip
import argparse
import hashlib
import json
import sys
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT/'.agents/skills/geo-single-cell-loader/scripts'
sys.path.insert(0,str(SCRIPTS))
sys.path.insert(0,str(ROOT/'runtime/r45'))
from common import dataset_paths, read_csv, write_csv, storage_root, runtime_output
from build_manifest import confirm
from runtime import RRunError, run_r45


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(output_root=None, resume_fixture=False):
    fixture = storage_root(ROOT,output_root) if output_root else ROOT/'.runtime'/('download-layout-'+uuid.uuid4().hex[:10])
    fixture.mkdir(parents=True,exist_ok=True)
    gse = 'GSE999999999'
    paths = dataset_paths(fixture,gse)
    fresh = not paths['dataset'].exists()
    if not fresh and (not resume_fixture or paths['final'].exists()):
        raise ValueError('Existing build fixture preserved; resume only an explicitly selected failed fixture without a final RDS')
    runtime_output(fixture,gse)
    trio = paths['raw']/'Synthetic1'
    trio.mkdir(parents=True,exist_ok=True)
    if fresh:
        with gzip.open(trio/'matrix.mtx.gz','wt') as stream:
            stream.write('%%MatrixMarket matrix coordinate integer general\n250 4 1000\n')
            for cell in range(1,5):
                for gene in range(1,251): stream.write(f'{gene} {cell} 2\n')
        with gzip.open(trio/'features.tsv.gz','wt') as stream:
            for gene in range(1,251): stream.write(f'ENSG{gene}\tGene{gene}\tGene Expression\n')
        with gzip.open(trio/'barcodes.tsv.gz','wt') as stream:
            stream.write('\n'.join('ACGTACGTACGTACG'+str(cell)+'-1' for cell in range(4))+'\n')
    area = paths['dataset'].relative_to(fixture).as_posix()
    row = dict(database=gse,sample='Synthetic1',tissue='Synthetic',disease='Synthetic',
        source_type='Tissue',group='Fixture',local_path=area+'/raw/Synthetic1',file_type='10x_mtx',
        count_source='counts',count_evidence='Generated integer fixture counts',
        metadata_evidence='Synthetic fixture only',files_json='[]')
    workflow = paths['workflow']; workflow.mkdir(exist_ok=True)
    report = workflow/'sample_report.csv'
    if fresh: write_csv(report,[row])
    approval = workflow/'group_confirmation.json'
    if fresh: approval.write_text(json.dumps(dict(confirmed_by='user',user_statement='SIMULATED FIXTURE ONLY',
        groups={'Synthetic1':'Fixture'})),encoding='utf-8')
    manifest = workflow/'sample_manifest.csv'
    if fresh: confirm(report,approval,manifest,fixture)
    else:
        from common import verify_confirmation
        verify_confirmation(manifest)
    raw_sha = {p.name:sha(p) for p in trio.iterdir()}
    receipt_sha = sha(manifest.with_suffix('.confirmation.json'))
    log = paths['dataset']/'regression.log'
    results = []
    def run(name,args):
        result = run_r45(SCRIPTS/name,args,cwd=fixture,stage='download_layout_'+name,log_path=log)
        results.append(dict(step=name,exit_code=result.returncode,duration_seconds=result.record.get('duration_seconds')))
    run_r45(route='base',cwd=fixture,stage='download_layout_health',log_path=log)
    relative = manifest.relative_to(fixture).as_posix()
    run('build_seurat.R',[fixture,relative])
    rds_sha = sha(paths['final'])
    run('validate_seurat.R',[fixture,relative,area+'/seurat_raw.rds','--record-completion'])
    # Preserve original confirmation for the same archived-artifact interface
    # used after successful audit upload; no manifest or receipt is rewritten.
    (paths['dataset']/'group_confirmation.json').write_bytes(approval.read_bytes())
    run('validate_seurat.R',['--existing',fixture,gse])
    run('report_download.R',['--existing',fixture,gse])
    report_sha = sha(paths['dataset']/'download_report.pdf')
    run('report_download.R',['--existing',fixture,gse])
    assert sha(paths['dataset']/'download_report.pdf')==report_sha
    assert sha(paths['final'])==rds_sha
    assert {p.name:sha(p) for p in trio.iterdir()}==raw_sha
    assert sha(manifest.with_suffix('.confirmation.json'))==receipt_sha
    summary = read_csv(paths['dataset']/'download_summary.csv')
    assert len(summary)==1 and summary[0]['cells']=='4' and summary[0]['genes_detected']=='250'
    assert len(read_csv(paths['dataset']/'metadata.csv'))==4
    assert (paths['dataset']/'download_report.pdf').read_bytes().startswith(b'%PDF-')
    # A changed output is refused; no existing CSV or RDS is overwritten.
    protected = paths['dataset']/'metadata.csv'
    protected.write_bytes(protected.read_bytes()+b'changed\n')
    changed = sha(protected)
    try:
        run('report_download.R',['--existing',fixture,gse])
    except RRunError as error:
        assert 'refusing overwrite' in error.output
        assert sha(protected)==changed and sha(paths['final'])==rds_sha
        results.append(dict(step='changed_metadata_refused',status='PASS'))
    else:
        raise AssertionError('Changed user artifact must stop without overwriting')
    evidence = dict(status='PASS',fixture=str(fixture),results=results,raw_sha256=raw_sha,
        rds_sha256=rds_sha,receipt_sha256=receipt_sha,pdf_sha256=report_sha,
        note='metadata.csv intentionally modified by final refusal test; genuine artifacts retained before that boundary')
    (paths['dataset']/'test_results.json').write_text(json.dumps(evidence,indent=2),encoding='utf-8')
    print(json.dumps(evidence,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-root',type=Path)
    parser.add_argument('--resume-fixture',action='store_true')
    main(**vars(parser.parse_args()))
