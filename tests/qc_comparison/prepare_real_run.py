"""Copy the exact historical raw input and persisted decision to a new QC root."""
import argparse
import json
import shutil
from pathlib import Path
from prepare_reference import sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('reference', type=Path)
    parser.add_argument('root', type=Path)
    parser.add_argument('gse')
    args = parser.parse_args()
    state = json.loads((args.reference / args.gse / 'qc/.workflow/state.json').read_text())
    raw = args.source / 'seurat_raw.rds'
    assert raw.is_file(), 'Missing true raw RDS; never substitute final QC'
    digest = sha256(raw)
    assert digest == state['raw_sha256'] == state['raw_sha256_after']
    target = args.root / 'data' / args.gse
    target.mkdir(parents=True, exist_ok=False)
    (target / 'qc').mkdir()
    shutil.copyfile(raw, target / raw.name)
    assert sha256(target / raw.name) == digest
    decision = args.reference / args.gse / 'qc/qc_decision.json'
    shutil.copyfile(decision, target / 'qc/qc_decision.json')
    (target / 'input_provenance.json').write_text(json.dumps(dict(
        source=str(raw), raw_sha256=digest, historical_raw_sha256=state['raw_sha256'],
        decision_sha256=sha256(decision), purpose='Isolated PCA/plot regression; historical files unchanged'),
        indent=2, ensure_ascii=False), encoding='utf-8')
    print('VERIFIED_RAW_COPY', digest)


if __name__ == '__main__':
    main()
