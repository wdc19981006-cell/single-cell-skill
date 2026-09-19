"""Complete verified technical input fields after download and before R build."""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from common import ROOT, digest, local, missing, read_csv, validate_manifest, verify_confirmation, workflow_path, write_csv
from pooled_mapping_probe import check_pooled
from tenx_structure_probe import probe_10x_trio
from text_schema_probe import probe_text


def run(manifest, root):
    root = Path(root).resolve()
    rows = read_csv(manifest)
    if not rows:
        raise ValueError('Empty manifest')
    gse = rows[0]['database']
    workflow_path(root, manifest, gse)
    workflow = local(root, f'data/{gse}/.workflow', f'data/{gse}')
    receipt = verify_confirmation(manifest)
    validate_manifest(rows, root, allow_pending_probes=True)
    text_results = []
    text_report = workflow / 'text_schema_probe.json'
    seen = set()
    try:
        for row in rows:
            if row['file_type'] != 'text' or row['local_path'] in seen:
                continue
            seen.add(row['local_path'])
            try:
                result = probe_text(row, root)
                result['status'] = 'VERIFIED_SAMPLE'
                text_results.append(result)
            except Exception as error:
                text_results.append({'local_path': row['local_path'], 'status': 'STOP', 'reason': str(error)})
                raise
            for other in rows:
                if other['local_path'] != row['local_path']:
                    continue
                technical = dict(delimiter=result['delimiter'], orientation=result['orientation'],
                                 feature_column=result['feature_column'],
                                 text_header_missing_id=str(result['header_missing_id']).lower(),
                                 matrix_rows=str(result['matrix_rows']),
                                 matrix_columns=str(result['matrix_columns']),
                                 estimated_dense_bytes=str(result['estimated_dense_bytes']),
                                 physical_ram_bytes=str(result['physical_ram_bytes']),
                                 text_reader=result['text_reader'],
                                 reader_selection_reason=result['reader_selection_reason'])
                for field, observed in technical.items():
                    if not missing(other.get(field)) and other[field] != observed:
                        raise ValueError(f'Shared text input has conflicting explicit {field}')
                    if missing(other.get(field)):
                        other[field] = observed
    finally:
        text_report.write_text(json.dumps(text_results, ensure_ascii=False, indent=2), encoding='utf-8')
    tenx_results = []
    tenx_report = workflow / 'tenx_structure_probe.json'
    seen = set()
    try:
        for row in rows:
            if row['file_type'] != '10x_mtx' or row['local_path'] in seen:
                continue
            seen.add(row['local_path'])
            try:
                result = probe_10x_trio(row, root)
                result['status'] = 'VERIFIED'
                tenx_results.append(result)
            except Exception as error:
                tenx_results.append({'local_path': row['local_path'], 'status': 'STOP', 'reason': str(error)})
                raise
            for other in rows:
                if other['local_path'] != row['local_path']:
                    continue
                observed = result['feature_name_fallback']
                if not missing(other.get('feature_name_fallback')) and other['feature_name_fallback'] != observed:
                    raise ValueError('Shared 10x input has conflicting explicit feature_name_fallback')
                if missing(other.get('feature_name_fallback')):
                    other['feature_name_fallback'] = observed
    finally:
        tenx_report.write_text(json.dumps(tenx_results, ensure_ascii=False, indent=2), encoding='utf-8')
    grouped = {}
    for row in rows:
        grouped.setdefault(row['local_path'], []).append(row)
    for shared in grouped.values():
        if len(shared) <= 1 or not all(missing(row.get('cell_map_path')) for row in shared):
            continue
        cell_map = check_pooled(shared, root, workflow)
        md5 = digest(local(root, cell_map, f'data/{gse}/.workflow'))
        for row in shared:
            row['cell_map_path'] = cell_map
            row['cell_map_md5'] = md5
            row['cell_map_cell_column'] = 'cell'
            row['cell_map_sample_column'] = 'sample'
    validate_manifest(rows, root)
    original = read_csv(manifest)
    if rows != original:
        pending = manifest.with_suffix('.probed.csv')
        write_csv(pending, rows)
        updated = dict(receipt, manifest_md5=digest(pending),
                       technical_probe_at=datetime.now(timezone.utc).isoformat(),
                       technical_probe_fields=('delimiter', 'orientation', 'feature_column', 'text_header_missing_id',
                                               'matrix_rows', 'matrix_columns', 'estimated_dense_bytes',
                                               'physical_ram_bytes', 'text_reader', 'reader_selection_reason',
                                               'feature_name_fallback', 'cell_map_path', 'cell_map_md5'))
        pending_receipt = manifest.with_suffix('.confirmation.pending.json')
        pending_receipt.write_text(json.dumps(updated, ensure_ascii=False, indent=2), encoding='utf-8')
        receipt_path = manifest.with_suffix('.confirmation.json')
        manifest_before, receipt_before = manifest.read_bytes(), receipt_path.read_bytes()
        try:
            pending.replace(manifest)
            pending_receipt.replace(receipt_path)
            verify_confirmation(manifest)
        except Exception:
            manifest.write_bytes(manifest_before)
            receipt_path.write_bytes(receipt_before)
            raise
    print('Pre-build text, 10x structure, and pooled mapping checks passed')
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path)
    parser.add_argument('--root', type=Path, default=ROOT)
    args = parser.parse_args()
    run(args.manifest, args.root)


if __name__ == '__main__':
    main()
