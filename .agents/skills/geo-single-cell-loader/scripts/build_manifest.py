"""Build a manifest only from an inspected report and explicit user grouping."""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from common import ROOT, digest, local, missing, read_csv, validate_manifest, verify_confirmation, write_csv

def confirm(report, confirmation, output, root):
    rows = read_csv(report)
    approval = json.loads(confirmation.read_text(encoding='utf-8-sig'))
    if approval.get('confirmed_by') != 'user' or not approval.get('user_statement', '').strip():
        raise ValueError('Explicit user statement and confirmed_by=user required')
    groups = approval.get('groups', {})
    if set(groups) != {r['sample'] for r in rows}: raise ValueError('Grouping keys must exactly match report samples')
    for row in rows:
        row['group'] = groups[row['sample']]
        if not missing(row.get('cell_map_path')):
            row['cell_map_md5'] = digest(local(root,row['cell_map_path'],'datasets/' + row['database']))
    validate_manifest(rows, root)
    if output.exists(): raise ValueError('Manifest already exists; choose a new path or explicitly archive it first')
    write_csv(output, rows)
    approval.update(manifest_md5=digest(output), confirmed_at=datetime.now(timezone.utc).isoformat())
    output.with_suffix('.confirmation.json').write_text(json.dumps(approval, ensure_ascii=False, indent=2), encoding='utf-8')
    return rows

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=ROOT)
    p.add_argument('--report', type=Path)
    p.add_argument('--confirmation', type=Path)
    p.add_argument('--output', type=Path)
    p.add_argument('--validate', type=Path)
    a = p.parse_args()
    if a.validate:
        validate_manifest(read_csv(a.validate), a.root.resolve()); verify_confirmation(a.validate)
        print('Manifest and confirmation valid')
    elif a.report and a.confirmation and a.output:
        confirm(a.report, a.confirmation, a.output, a.root.resolve()); print(a.output)
    else: p.error('Provide --validate, or --report --confirmation --output')

if __name__ == '__main__': main()
