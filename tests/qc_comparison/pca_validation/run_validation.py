"""Isolated QC using current orchestration, with a test-only diagnostic observer."""
import argparse
import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]


def sha256(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('gse')
    parser.add_argument('source', type=Path)
    parser.add_argument('root', type=Path)
    parser.add_argument('--mode', choices=('latest','legacy_pc_control'), default='latest')
    args = parser.parse_args()
    source, root = args.source.resolve(), args.root.resolve()
    raw = source / 'seurat_raw.rds'
    if not raw.is_file():
        raise SystemExit('未找到 seurat_raw.rds；停止，不使用 QC 后对象替代。')
    expected = json.loads((source/'qc/workflow_state.json').read_text(encoding='utf-8'))
    digest = sha256(raw)
    assert digest == expected['raw_sha256'] == expected['raw_sha256_after']
    dataset = root / 'data' / args.gse
    dataset.mkdir(parents=True, exist_ok=False)
    (dataset/'qc').mkdir()
    shutil.copyfile(raw,dataset/raw.name)
    shutil.copyfile(source/'qc/qc_decision.json',dataset/'qc/qc_decision.json')
    assert sha256(dataset/raw.name)==digest
    protected = [raw,source/'seurat.all_qc.rds',source/'seurat.qc.rds',
                 source/'qc/seurat_qc.rds',source/'qc/qc_report.csv',source/'qc/qc_decision.json',
                 source/'group_confirmation.json']
    provenance = dict(started_at=datetime.now(timezone.utc).isoformat(),mode=args.mode,
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
        input_raw_sha256=digest,source=str(source),
        protected_files={str(p):sha256(p) for p in protected},
        no_discovery_download_build=True,existing_user_decision_reused=True)
    # Audit starts before any real R computation. The final publisher maps this
    # evidence into the project's ten-file independent-QC audit contract.
    (dataset/'validation_provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    spec=importlib.util.spec_from_file_location('diagnostic_qc_workflow',REPO/'qc/run_workflow.py')
    workflow=importlib.util.module_from_spec(spec);spec.loader.exec_module(workflow)
    actual_run=workflow.run_r45
    def observed_run(script=None,argv=(),**kwargs):
        if kwargs.get('stage')=='qc':
            script=Path(__file__).with_name('run_diagnostics.R')
            argv=[REPO,root,args.gse,dataset/'diagnostics',args.mode]
        return actual_run(script,argv,**kwargs)
    workflow.run_r45=observed_run
    try:
        result=workflow.run(args.gse,root)
        actual_run(REPO/'qc/validate_qc.R',[root,args.gse],stage='diagnostic_final_validation',cwd=REPO,
                   log_path=dataset/'qc/execution.log')
        print(result,flush=True)
    finally:
        after={str(p):sha256(p) for p in protected}
        provenance.update(protected_files_after=after,
            historical_files_unchanged=after==provenance['protected_files'],
            raw_copy_sha256_after=sha256(dataset/raw.name))
        (dataset/'validation_provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
        assert provenance['historical_files_unchanged'] and provenance['raw_copy_sha256_after']==digest


if __name__=='__main__':
    main()
