"""UTF-8 user report, derived from inspection before approval and manifest after it."""
import argparse
import json
from pathlib import Path
from common import ROOT, OPTIONAL, dataset_paths, missing, read_csv, validate_manifest, verify_confirmation

LINE = '=' * 60

def render(root, gse, rows, inspection=None, confirmed=False):
    paths = dataset_paths(root, gse)
    inspection = inspection or {}
    fields = inspection.get('series', {}).get('fields', {})
    def values(key):
        return '; '.join(fields.get(key, [])) or 'Not yet established from public evidence.'
    def distinct(key):
        return '; '.join(dict.fromkeys(str(r[key]) for r in rows if not missing(r.get(key)))) or 'Not yet established from public evidence.'
    complete = inspection.get('complete', True)
    status = 'GROUP_CONFIRMED' if confirmed else ('WAITING_FOR_GROUP_CONFIRMATION' if complete else 'INCOMPLETE_INSPECTION')
    formats = list(dict.fromkeys(f.get('format', '') for f in inspection.get('files', []) if f.get('format') not in ('unknown', 'raw_reads')))
    lines = [f'STATUS: {status}', '', gse, LINE, '', 'Study:', values('Series_title'), '',
             'Study description:', values('Series_summary'), '', 'Sequencing type:',
             inspection.get('sequencing_type') or 'See source sample descriptions; modality requires evidence review.', '',
             'Tissue:', distinct('tissue'), '', 'Processed expression format:',
             distinct('file_type') if confirmed else (', '.join(formats) + ' (candidates; contents not yet verified)' if formats else 'Unresolved'), '',
             'Total samples:', str(inspection.get('total_samples', len(rows))), '',
             'Inspected samples:', str(len(rows)), '', 'Group:',
             'User-confirmed; see each sample and GROUP SUMMARY.' if confirmed else 'Not assigned. Waiting for user confirmation.', '',
             'Source data directory:', 'raw/ contains source expression files downloaded from GEO, unchanged by this Skill.',
             'raw/ does not necessarily mean sequencing FASTQ raw reads.', '',
             'Evidence review:', 'Unresolved metadata must be reviewed against GEO, publications and public sample metadata.', '',
             LINE, 'SAMPLE INFORMATION', LINE, '']
    if not complete:
        lines += ['INCOMPLETE DEVELOPMENT REPORT: cannot be used to confirm groups or download expression.', '']
    optional = [key for key in OPTIONAL if any(not missing(r.get(key)) for r in rows)]
    for i, row in enumerate(rows, 1):
        lines += [f'Sample {i}', '-' * 60]
        for key in ('sample', 'author_sample', 'tissue', 'disease', 'source_type', 'sample_description') + tuple(optional) + (('group',) if confirmed else ()):
            lines += [key + ':', str(row[key]) if not missing(row.get(key)) else 'Not reported', '']
    if confirmed:
        lines += [LINE, 'GROUP SUMMARY', LINE, '']
        for group in dict.fromkeys(r['group'] for r in rows):
            lines += [group, '-' * 60] + [r['sample'] for r in rows if r['group'] == group] + ['']
    paths['dataset'].mkdir(parents=True, exist_ok=True)
    paths['info'].write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return paths['info']

def refresh(root, gse):
    paths = dataset_paths(root, gse)
    inspection_path = paths['workflow'] / 'inspection.json'
    inspection = json.loads(inspection_path.read_text(encoding='utf-8')) if inspection_path.exists() else {}
    manifest = paths['workflow'] / 'sample_manifest.csv'
    if manifest.exists():
        verify_confirmation(manifest)
        rows = validate_manifest(read_csv(manifest), root)
        if paths['final'].exists():
            raise ValueError('Final RDS exists; do not reset a completed sample_info report')
        return render(root, gse, rows, inspection, confirmed=True)
    return render(root, gse, read_csv(paths['workflow'] / 'sample_report.csv'), inspection)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('gse'); parser.add_argument('--root', type=Path, default=ROOT)
    args = parser.parse_args()
    print(refresh(args.root.resolve(), args.gse))
