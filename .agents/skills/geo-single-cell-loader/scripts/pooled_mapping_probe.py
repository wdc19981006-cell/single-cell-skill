"""Check pooled cell/sample evidence without reading the expression matrix."""
import csv
import gzip
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from urllib.parse import urlparse

from common import local

SUFFIX = re.compile(r'[-_](\d+)$')
REQUIRED_SOURCES = {'geo_supplementary', 'publication_supplement', 'author_repository'}


def _barcodes(row, root):
    gse = row['database']
    path = local(root, row['local_path'], f'data/{gse}/raw')
    kind = row['file_type']
    if kind == '10x_h5':
        import h5py
        with h5py.File(path, 'r') as handle:
            if 'matrix/barcodes' not in handle:
                raise ValueError('10x H5 lacks matrix/barcodes')
            values = handle['matrix/barcodes'][:]
        return [value.decode('utf-8') if isinstance(value, bytes) else str(value) for value in values]
    if kind == '10x_mtx':
        names = [path / (name + suffix) for name in ('barcodes.tsv',) for suffix in ('', '.gz')]
        found = [candidate for candidate in names if candidate.is_file()]
        if len(found) != 1:
            raise ValueError('Pooled 10x trio has missing or ambiguous barcodes.tsv')
        opener = gzip.open if found[0].name.endswith('.gz') else open
        with opener(found[0], 'rt', encoding='utf-8-sig') as handle:
            return [line.rstrip('\r\n').split('\t')[0] for line in handle]
    if kind == 'text':
        if row.get('orientation') != 'genes_by_cells':
            raise ValueError('Pooled text barcode inspection needs genes_by_cells layout')
        opener = gzip.open if path.name.endswith('.gz') else open
        with opener(path, 'rt', encoding='utf-8-sig', newline='') as handle:
            header = handle.readline(8 * 1024 * 1024 + 1)
        if len(header) > 8 * 1024 * 1024:
            raise ValueError('Pooled text header exceeds barcode probe limit')
        delimiter = {'comma': ',', 'tab': '\t', 'space': ' '}.get(row.get('delimiter'))
        if delimiter is None:
            raise ValueError('Pooled text delimiter is not confirmed')
        cells = next(csv.reader([header], delimiter=delimiter))
        return cells if row.get('text_header_missing_id') == 'true' else cells[1:]
    if kind == 'h5ad':
        import h5py
        with h5py.File(path, 'r') as handle:
            if 'obs/_index' not in handle:
                raise ValueError('H5AD lacks directly readable obs/_index')
            values = handle['obs/_index'][:]
        return [value.decode('utf-8') if isinstance(value, bytes) else str(value) for value in values]
    raise ValueError('Barcode inspection is unavailable for this pooled input type')


def _candidate_rows(path):
    opener = gzip.open if path.name.endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))


def _gsm_ids(rows):
    ids = set()
    for row in rows:
        linked = [value.strip() for value in (row.get('source_gsm') or '').split(';') if value.strip()]
        if linked:
            ids.update(linked)
        elif re.fullmatch(r'GSM\d+', row['sample']):
            ids.add(row['sample'])
    return ids


def check_pooled(rows, root, workflow):
    """Return a verified canonical cell map or stop with an audit of the checks."""
    first = rows[0]
    gse = first['database']
    result = {'local_path': first['local_path'], 'manifest_samples': [r['sample'] for r in rows],
              'source_check_order': ['matrix_internal', 'geo_supplementary',
                                     'publication_supplement', 'author_repository'],
              'status': 'STOP', 'candidate_results': []}
    report = workflow / 'pooled_mapping_probe.json'
    try:
        cells = _barcodes(first, root)
        if not cells or any(not cell for cell in cells) or len(cells) != len(set(cells)):
            raise ValueError('Pooled matrix has blank or duplicate barcode IDs')
        suffixes = Counter(match.group(1) for cell in cells if (match := SUFFIX.search(cell)))
        selected_samples = {row['sample'] for row in rows}
        receipt = json.loads((workflow / 'sample_manifest.confirmation.json').read_text(encoding='utf-8'))
        explicit_subset = receipt.get('selected_samples')
        if explicit_subset is not None:
            if (not isinstance(explicit_subset, list) or len(explicit_subset) != len(set(explicit_subset))
                    or set(explicit_subset) != selected_samples):
                raise ValueError('Confirmed selected samples differ from pooled manifest rows')
            with (workflow / 'sample_report.csv').open(encoding='utf-8-sig', newline='') as handle:
                all_report_rows = list(csv.DictReader(handle))
            source_samples = {row['sample'] for row in all_report_rows}
            if len(source_samples) != len(all_report_rows) or not selected_samples <= source_samples:
                raise ValueError('Selected samples are absent from the complete Stage A report')
            gsm_ids = _gsm_ids(all_report_rows)
        else:
            source_samples = selected_samples
            gsm_ids = _gsm_ids(rows)
        result.update(barcode_count=len(cells), barcode_suffix_counts=dict(sorted(suffixes.items())),
                      gsm_count=len(gsm_ids), gsm_ids=sorted(gsm_ids),
                      explicit_subset=explicit_subset is not None)
        evidence_file = workflow / 'mapping_evidence.json'
        evidence = json.loads(evidence_file.read_text(encoding='utf-8')) if evidence_file.exists() else {}
        checks = evidence.get('checks', [])
        result['source_checks'] = checks
        if suffixes and gsm_ids and len(suffixes) != len(gsm_ids):
            raise ValueError(f'Barcode groups ({len(suffixes)}) and GSM count ({len(gsm_ids)}) differ; no suffix-order mapping is allowed')
        covered = {entry.get('source_type') for entry in checks if entry.get('result') and
                   (entry.get('url') or entry.get('reason'))}
        if not REQUIRED_SOURCES <= covered:
            raise ValueError('Review GEO supplementary metadata, publication supplements, and author repository before pooled mapping')
        for candidate in evidence.get('candidate_maps', []):
            entry = {'source_type': candidate.get('source_type'), 'source_url': candidate.get('source_url'),
                     'status': 'rejected'}
            result['candidate_results'].append(entry)
            try:
                url = urlparse(candidate.get('source_url', ''))
                if url.scheme != 'https' or not url.hostname or url.username or url.password:
                    raise ValueError('Candidate map needs an HTTPS public source URL')
                if not any(check.get('source_type') == candidate.get('source_type') and
                           check.get('url') == candidate['source_url'] for check in checks):
                    raise ValueError('Candidate map source is absent from reviewed source checks')
                relative = candidate['local_path']
                path = local(root, relative, f'data/{gse}/.workflow')
                mapping = _candidate_rows(path)
                cell_key = candidate.get('cell_column', 'cell')
                sample_key = candidate.get('sample_column', 'sample')
                aliases = candidate.get('alias_to_sample', {})
                if aliases and not candidate.get('alias_evidence'):
                    raise ValueError('Alias-to-sample mapping lacks public evidence')
                if not mapping or any(cell_key not in item or sample_key not in item for item in mapping):
                    raise ValueError('Candidate map lacks configured cell/sample columns')
                mapped_cells = [item[cell_key] for item in mapping]
                mapped_samples = [aliases.get(item[sample_key], item[sample_key]) for item in mapping]
                mapped_set = set(mapped_samples)
                if (len(mapped_cells) != len(cells) or len(set(mapped_cells)) != len(cells) or
                        set(mapped_cells) != set(cells) or any(not value for value in mapped_samples) or
                        not selected_samples <= mapped_set or not mapped_set <= source_samples or
                        (explicit_subset is None and mapped_set != selected_samples)):
                    raise ValueError('Candidate map is not exhaustive for matrix cells and manifest samples')
                canonical = workflow / ('cell_map_' + hashlib.sha256(first['local_path'].encode()).hexdigest()[:12] + '.csv')
                with canonical.open('w', encoding='utf-8', newline='') as handle:
                    writer = csv.DictWriter(handle, fieldnames=['cell', 'sample'])
                    writer.writeheader()
                    writer.writerows({'cell': cell, 'sample': sample}
                                     for cell, sample in zip(mapped_cells, mapped_samples))
                result.update(status='VERIFIED', cell_map_path=f'data/{gse}/.workflow/{canonical.name}',
                              mapping_source=candidate['source_url'],
                              excluded_samples=sorted(mapped_set - selected_samples),
                              excluded_cells=sum(sample not in selected_samples for sample in mapped_samples))
                entry['status'] = 'verified'
                return result['cell_map_path']
            except (KeyError, OSError, ValueError) as error:
                entry['reason'] = str(error)
        raise ValueError('No verifiable exhaustive cell-to-sample mapping for pooled matrix')
    except Exception as error:
        result['stop_reason'] = str(error)
        raise
    finally:
        report.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
