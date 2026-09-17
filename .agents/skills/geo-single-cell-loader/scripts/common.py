"""Repository-relative paths, strict manifests and processed-file routing."""
import csv
import ctypes
import hashlib
import json
import math
import os
import re
from pathlib import Path
from urllib.parse import urlparse

ROOT = Path(__file__).resolve().parents[4]
REQUIRED = ('database', 'sample', 'tissue', 'disease', 'source_type', 'group')
OPTIONAL = ('patient', 'specimen', 'cohort', 'treatment')

def local(root, value, area=None):
    p = Path(value)
    if p.is_absolute() or re.match(r'^[A-Za-z]:', value) or '\\' in value or '..' in p.parts:
        raise ValueError('Use repository-relative forward-slash paths: ' + value)
    target = (root / p).resolve()
    base = (root / area).resolve() if area else root.resolve()
    if not base.is_relative_to(root.resolve()) or (area and base != root.resolve() / area):
        raise ValueError('Permitted directory redirects outside its declared location')
    if not target.is_relative_to(base) or target == base:
        raise ValueError('Path escapes permitted directory: ' + value)
    return target

def read_csv(path):
    with open(path, encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))

def write_csv(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(dict.fromkeys(k for row in rows for k in row))
    with open(path, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

def missing(value):
    return value is None or str(value).strip().lower() in ('', 'na', 'nan', 'null', 'none')

def file_type(name):
    name = name.lower().split('?')[0]
    if name.endswith(('.h5ad', '.h5ad.gz')): return 'h5ad'
    if name.endswith(('.h5', '.hdf5')): return '10x_h5_candidate'
    if name.endswith(('.mtx', '.mtx.gz')): return '10x_mtx_candidate'
    if name.endswith(('.rds', '.rds.gz')): return 'rds'
    if name.endswith(('.csv', '.tsv', '.txt', '.csv.gz', '.tsv.gz', '.txt.gz')): return 'text_candidate'
    if name.endswith(('.tar', '.tar.gz', '.tgz', '.zip')): return 'archive'
    if name.endswith(('.fastq', '.fastq.gz', '.sra')): return 'raw_reads'
    return 'unknown'


TRANSFORMED_EXPRESSION = re.compile(r'(^|[_.-])(normalized|normalised|log2?tpm|log1p|log|sct|integrated|scaled|tpm)([_.-]|$)', re.I)
BINARY_COUNTS = {'rds', '10x_h5_candidate', '10x_mtx_candidate', 'h5ad'}


def _candidate_format(item):
    return file_type(urlparse(item['url']).path)


def _transformed(item):
    return bool(TRANSFORMED_EXPRESSION.search(urlparse(item['url']).path.rsplit('/', 1)[-1]))


def _local_binary(item, root):
    path = item.get('local_path')
    if not path: return False
    root = Path(root).resolve()
    target = (root / path).resolve()
    return target.is_relative_to(root) and target.is_file()


def total_physical_memory_bytes():
    """Physical RAM, not available RAM or an R/Python allocation limit."""
    if os.name == 'nt':
        class MemoryStatus(ctypes.Structure):
            _fields_ = [('length', ctypes.c_ulong), ('load', ctypes.c_ulong)] + [
                (name, ctypes.c_ulonglong) for name in ('total_phys', 'avail_phys', 'total_page',
                                                        'avail_page', 'total_virtual', 'avail_virtual', 'extended')]
        status = MemoryStatus()
        status.length = ctypes.sizeof(status)
        if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
            raise OSError('Cannot determine total physical RAM')
        return int(status.total_phys)
    try:
        return int(os.sysconf('SC_PHYS_PAGES') * os.sysconf('SC_PAGE_SIZE'))
    except (AttributeError, OSError, ValueError):
        raise OSError('Cannot determine total physical RAM') from None


def dense_rds_memory_safe(estimated_bytes, total_ram_bytes, fraction=0.25):
    """Require the actual loaded object to fit within the RAM fraction."""
    return (isinstance(estimated_bytes, (int, float)) and math.isfinite(estimated_bytes)
            and isinstance(total_ram_bytes, (int, float)) and math.isfinite(total_ram_bytes)
            and estimated_bytes > 0 and total_ram_bytes > 0
            and estimated_bytes <= total_ram_bytes * fraction)


def rank_raw_count_candidates(candidates, root=ROOT):
    """Probe only evidenced representations of the same raw counts."""
    root = Path(root)
    raw = [item for item in candidates if not _transformed(item) and _candidate_format(item) in BINARY_COUNTS | {'text_candidate'}]
    if len(raw) > 1:
        identities = {item.get('raw_counts_id') for item in raw}
        if len(identities) != 1 or not next(iter(identities)) or any(not item.get('same_counts_evidence') for item in raw):
            raise ValueError('Multiple candidates require one evidenced raw_counts_id before format selection')
    def priority(item):
        kind = _candidate_format(item)
        if kind == 'text_candidate': return 3
        if item.get('public_sparse_raw_evidence'): return 0
        if _local_binary(item, root): return 1
        return 2
    return sorted(raw, key=priority)


def select_verified_raw_counts(candidates, validate, *, root=ROOT, total_ram_bytes=None, audit_path=None):
    """Select verified sparse binary, safe dense RDS, then streaming TXT.

    `validate` inspects content and returns raw_counts, sparse, object_class and
    estimated_memory_bytes. Unknown remote binary is not probed when usable raw
    counts are already verified, unless public evidence explicitly says sparse.
    """
    ranked = rank_raw_count_candidates(candidates, root)
    if total_ram_bytes is None:
        try: total_ram_bytes = total_physical_memory_bytes()
        except OSError: total_ram_bytes = None
    if total_ram_bytes is not None and (not isinstance(total_ram_bytes, (int, float))
                                        or not math.isfinite(total_ram_bytes) or total_ram_bytes <= 0):
        raise ValueError('Invalid total physical RAM')
    attempts = []
    limit = int(total_ram_bytes * 0.25) if total_ram_bytes is not None else None
    for item in candidates:
        if _transformed(item):
            attempts.append(dict(url=item['url'], candidate_format=_candidate_format(item), object_class=None,
                                 estimated_memory_bytes=None, total_ram_bytes=total_ram_bytes, selected_route=None,
                                 status='excluded', reason='normalized/log/SCT/scaled expression is not raw counts',
                                 fallback_reason='normalized/log/SCT/scaled expression is not raw counts'))
    def record(item, status, reason, evidence=None, route=None):
        evidence = evidence or {}
        entry = dict(url=item['url'], candidate_format=_candidate_format(item),
                     object_class=evidence.get('object_class'),
                     estimated_memory_bytes=evidence.get('estimated_memory_bytes'),
                     structure_evidence=item.get('structure_evidence'),
                     total_ram_bytes=total_ram_bytes, selected_route=route, status=status,
                     reason=reason, fallback_reason=reason if status in ('rejected', 'skipped') else None)
        attempts.append(entry)
        return entry
    def finish(item, route, entry, reason):
        entry.update(status='selected', selected_route=route, reason=reason, fallback_reason=None)
        if audit_path is not None:
            path = Path(audit_path)
            path.write_text(json.dumps(dict(selected_url=item['url'], selected_route=route,
                                            total_ram_bytes=total_ram_bytes, dense_limit_bytes=limit,
                                            candidates=attempts), indent=2), encoding='utf-8')
        return item, attempts
    binary = [item for item in ranked if _candidate_format(item) != 'text_candidate']
    text = [item for item in ranked if _candidate_format(item) == 'text_candidate']
    if any(item.get('verified_raw_counts') and not item.get('verification_evidence') for item in text):
        raise ValueError('Previously verified TXT requires verification_evidence')
    verified_txt = any(item.get('verified_raw_counts') for item in text)
    dense_choice = None
    for item in binary:
        kind = _candidate_format(item)
        usable_other = dense_choice is not None or verified_txt
        if not (_local_binary(item, Path(root)) or item.get('public_sparse_raw_evidence') or not usable_other):
            record(item, 'skipped', 'Unknown remote binary; a verified raw-count route already exists')
            continue
        try:
            if item.get('validated_evidence') is not None:
                if not item.get('structure_evidence'):
                    raise ValueError('Cached structure requires structure_evidence')
                evidence = item['validated_evidence']
            else:
                evidence = validate(item)
            if not evidence.get('raw_counts'):
                raise ValueError('Candidate is not verified raw counts')
            size = evidence.get('estimated_memory_bytes')
            if not isinstance(size, (int, float)) or not math.isfinite(size) or size <= 0:
                raise ValueError('Binary candidate has no measured in-memory size')
            if evidence.get('sparse'):
                entry = record(item, 'eligible', evidence.get('reason', 'verified sparse raw counts'), evidence)
                if dense_choice is not None:
                    dense_choice[1].update(status='not_selected', reason='Verified sparse binary has higher priority')
                route = {'rds':'sparse_rds','10x_h5_candidate':'sparse_h5',
                         '10x_mtx_candidate':'sparse_mtx','h5ad':'sparse_h5ad'}[kind]
                return finish(item, route, entry, evidence.get('reason', 'verified sparse raw counts'))
            classes = evidence.get('object_class') or ''
            classes = classes if isinstance(classes, (list, tuple)) else str(classes).split(',')
            if kind == 'rds' and any(name in ('data.frame', 'matrix') for name in classes):
                if not dense_rds_memory_safe(size, total_ram_bytes):
                    reason = ('Total physical RAM is unknown; dense RDS rejected' if limit is None else
                              f'Dense RDS uses {size:g} bytes; 25% RAM limit is {limit} bytes')
                    record(item, 'rejected', reason, evidence)
                    continue
                entry = record(item, 'eligible', 'Dense RDS fits within 25% of total RAM', evidence)
                if dense_choice is None: dense_choice = (item, entry)
                continue
            raise ValueError('Binary candidate is neither sparse counts nor a supported dense RDS')
        except Exception as error:
            record(item, 'rejected', str(error))
    if dense_choice is not None:
        return finish(dense_choice[0], 'dense_rds', dense_choice[1], 'Verified dense raw counts fit within 25% of total RAM')
    for item in text:
        try:
            evidence = ({'raw_counts': True, 'object_class': 'streaming TXT',
                         'reason': item.get('verification_reason', item.get('verification_evidence'))}
                        if item.get('verified_raw_counts') else validate(item))
            if not evidence.get('raw_counts'):
                raise ValueError('TXT candidate is not verified raw counts')
        except Exception as error:
            record(item, 'rejected', str(error))
            continue
        entry = record(item, 'eligible', evidence.get('reason', 'Verified streaming TXT raw counts'), evidence)
        return finish(item, 'streaming_txt', entry, evidence.get('reason', 'Verified streaming TXT raw counts'))
    if audit_path is not None:
        Path(audit_path).write_text(json.dumps(dict(selected_url=None, selected_route=None,
                                                   total_ram_bytes=total_ram_bytes, dense_limit_bytes=limit,
                                                   candidates=attempts), indent=2), encoding='utf-8')
    raise ValueError('No verified raw-count input; attempts: ' + str(attempts))

def choose_filtered(names):
    """Caller must first scope candidates to ONE evidenced sample."""
    filtered = [n for n in names if 'filtered_feature_bc_matrix' in n.lower()]
    choices = filtered or names
    if len(choices) != 1: raise ValueError('Ambiguous input; inspect and explicitly map one sample')
    return choices[0]

def detect_trio(directory):
    directory = Path(directory)
    found = []
    for options in [('matrix.mtx',), ('features.tsv', 'genes.tsv'), ('barcodes.tsv',)]:
        hits = [directory / (n + suffix) for n in options for suffix in ('', '.gz') if (directory / (n + suffix)).is_file()]
        if len(hits) != 1: raise ValueError('Missing or ambiguous 10x member: ' + str(options))
        found.append(hits[0])
    return found

def validate_manifest(rows, root):
    if not rows: raise ValueError('Empty manifest')
    seen, databases, destinations = set(), set(), {}
    for row in rows:
        for key in REQUIRED + ('local_path', 'file_type', 'count_source', 'count_evidence', 'metadata_evidence'):
            if missing(row.get(key)): raise ValueError('Missing ' + key)
        sample = row['sample']
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]*', sample): raise ValueError('Unsafe sample ID')
        if sample in seen: raise ValueError('Duplicate sample: ' + sample)
        seen.add(sample)
        if not re.fullmatch(r'GSE\d+', row['database']): raise ValueError('Invalid GSE')
        databases.add(row['database'])
        area = 'data/' + row['database']
        local(root, row['local_path'], area + '/raw')
        if row['file_type'] not in ('10x_mtx', '10x_h5', 'h5ad', 'text', 'rds'): raise ValueError('Unsupported file_type')
        if row['file_type'] == 'h5ad' and not row['local_path'].lower().endswith(('.h5ad', '.h5ad.gz')): raise ValueError('H5AD extension mismatch')
        if row['file_type'] == '10x_h5' and row['local_path'].lower().endswith(('.h5ad', '.h5ad.gz')): raise ValueError('H5AD is not 10x H5')
        if row['file_type'] == 'text' and row.get('orientation') not in ('genes_by_cells', 'cells_by_genes'): raise ValueError('Text orientation must be inspected')
        for field in OPTIONAL:
            if not missing(row.get(field)) and missing(row.get(field + '_evidence')): raise ValueError('Missing evidence for ' + field)
        if not missing(row.get('cell_map_path')):
            permitted = area + ('/cell_map' if row['cell_map_path'].startswith(area + '/cell_map/') else '/.workflow')
            cellmap = local(root, row['cell_map_path'], permitted)
            if missing(row.get('cell_map_md5')) or not cellmap.is_file() or digest(cellmap) != row['cell_map_md5']:
                raise ValueError('Cell map missing or changed; review and reconfirm')
        files = json.loads(row.get('files_json') or '[]')
        if not isinstance(files, list): raise ValueError('files_json must be a list')
        for item in files:
            local(root, item['local_path'], area + '/raw')
            u = urlparse(item['url'])
            if u.scheme != 'https' or not u.hostname or u.username or u.password: raise ValueError('HTTPS public URL required')
            if u.query and re.search(r'token|signature|credential|key=', u.query, re.I): raise ValueError('Do not store credentials in URLs')
            previous = destinations.setdefault(item['local_path'], item)
            if previous != item: raise ValueError('Conflicting download destination')
    if len(databases) != 1: raise ValueError('One GSE per manifest')
    by_path = {}
    for row in rows: by_path.setdefault(row['local_path'], []).append(row)
    for shared in by_path.values():
        if len(shared) > 1 and (any(missing(r.get('cell_map_path')) for r in shared) or len({r['cell_map_path'] for r in shared}) != 1):
            raise ValueError('Shared matrix requires one explicit cell-to-sample mapping')
    return rows

def digest(path):
    return hashlib.md5(path.read_bytes()).hexdigest()

def verify_confirmation(path):
    receipt = path.with_suffix('.confirmation.json')
    value = json.loads(receipt.read_text(encoding='utf-8'))
    if value.get('manifest_md5') != digest(path) or not value.get('user_statement', '').strip() or value.get('confirmed_by') != 'user':
        raise ValueError('Manifest changed or missing user confirmation; reconfirm before continuing')
    rows = read_csv(path)
    if value.get('groups') != {r['sample']: r['group'] for r in rows}:
        raise ValueError('Confirmed groups differ from manifest')
    selected = value.get('selected_samples')
    if selected is not None and (not isinstance(selected, list) or
                                 any(not isinstance(sample, str) for sample in selected) or
                                 len(selected) != len(set(selected)) or
                                 set(selected) != {r['sample'] for r in rows}):
        raise ValueError('Confirmed selected samples differ from manifest')
    return value

def dataset_paths(root, gse):
    if not re.fullmatch(r'GSE\d+', gse):
        raise ValueError('Expected GSE accession')
    base = local(root, 'data/' + gse, 'data')
    return dict(dataset=base, raw=local(root, f'data/{gse}/raw', f'data/{gse}'),
                workflow=local(root, f'data/{gse}/.workflow', f'data/{gse}'),
                final=local(root, f'data/{gse}/seurat_raw.rds', f'data/{gse}'),
                info=local(root, f'data/{gse}/sample_info.txt', f'data/{gse}'))

def workflow_path(root, path, gse):
    relative = Path(path).resolve().relative_to(root.resolve()).as_posix()
    return local(root, relative, f'data/{gse}/.workflow')
