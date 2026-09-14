"""Build a manifest only from an inspected report and explicit user grouping."""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from common import ROOT, dataset_paths, digest, local, missing, read_csv, validate_manifest, verify_confirmation, workflow_path, write_csv
from sample_info import render

def confirm(report, confirmation, output, root):
    rows = read_csv(report)
    if not rows: raise ValueError('Empty report')
    paths = dataset_paths(root, rows[0]['database'])
    for path in (report, confirmation, output): workflow_path(root, path, rows[0]['database'])
    for path, name in ((report,'sample_report.csv'),(confirmation,'group_confirmation.json'),(output,'sample_manifest.csv')):
        if path.resolve() != paths['workflow'] / name:
            raise ValueError('Use canonical .workflow/' + name + '; archive previous versions before rebuilding')
    inspection_file = paths['workflow'] / 'inspection.json'
    inspection = json.loads(inspection_file.read_text(encoding='utf-8')) if inspection_file.exists() else {}
    if not inspection and rows[0]['database'] != 'GSE999999999':
        raise ValueError('Complete inspection required before confirmation')
    if inspection:
        discovered = {r['accession'] for r in inspection['samples']}
        reported = {r['sample'] for r in rows}
        if reported != discovered:
            # Author matrix IDs are allowed only with explicit GEO ownership.
            reported = {gsm.strip() for r in rows for gsm in r.get('source_gsm', '').split(';') if gsm.strip()}
        if not inspection.get('complete') or reported != discovered:
            raise ValueError('Incomplete inspection or report sample set changed; review complete discovery')
    if paths['final'].exists(): raise ValueError('Final RDS exists; archive explicitly before reconfirming')
    approval = json.loads(confirmation.read_text(encoding='utf-8-sig'))
    if approval.get('confirmed_by') != 'user' or not approval.get('user_statement', '').strip():
        raise ValueError('Explicit user statement and confirmed_by=user required')
    groups = approval.get('groups', {})
    if set(groups) != {r['sample'] for r in rows}: raise ValueError('Grouping keys must exactly match report samples')
    for row in rows:
        row['group'] = groups[row['sample']]
        if not missing(row.get('cell_map_path')):
            row['cell_map_md5'] = digest(local(root,row['cell_map_path'],'data/' + row['database'] + '/.workflow'))
    validate_manifest(rows, root)
    if output.exists(): raise ValueError('Manifest already exists; choose a new path or explicitly archive it first')
    write_csv(output, rows)
    approval.update(manifest_md5=digest(output), confirmed_at=datetime.now(timezone.utc).isoformat())
    output.with_suffix('.confirmation.json').write_text(json.dumps(approval, ensure_ascii=False, indent=2), encoding='utf-8')
    render(root, rows[0]['database'], rows, inspection, confirmed=True)
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
