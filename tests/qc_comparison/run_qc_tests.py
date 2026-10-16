"""Run focused QC regressions and existing-artifact validation through R45."""
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'runtime/r45'))
from runtime import run_r45


def main():
    results = []
    env = dict(os.environ, QC_TEST_RETAIN_FIXTURE='1', PYTHONIOENCODING='utf-8')
    destination = ROOT / 'tests/qc_comparison/qc_test_results.json'
    for name in ('test_qc_pca.R', 'test_qc.R', 'test_qc_real_df.R'):
        result = run_r45(ROOT / 'tests' / name, [ROOT], stage=name, cwd=ROOT, env=env,
                         log_path=ROOT / '.runtime/qc-focused-tests.log')
        results.append(dict(test=name,status='PASS',exit_code=result.returncode,
                            duration_seconds=result.record.get('duration_seconds')))
        destination.write_text(json.dumps(results, indent=2), encoding='utf-8')
        print('PASS', name, flush=True)
    fixture = Path((ROOT / 'data/qc_fixture_path.txt').read_text().strip())
    run_r45(ROOT / 'qc/validate_qc.R', [fixture, 'GSE999999997'],
            stage='QC_existing_artifact_validation',cwd=ROOT,env=env,
            log_path=ROOT / '.runtime/qc-focused-tests.log')
    results.append(dict(test='QC_existing_artifact_validation',status='PASS',exit_code=0))
    destination.write_text(json.dumps(results, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
