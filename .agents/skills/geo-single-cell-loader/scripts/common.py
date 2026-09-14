"""Repository-relative paths, strict manifests and processed-file routing."""
import csv
import hashlib
import json
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
    if name.endswith('.rds'): return 'rds'
    if name.endswith(('.csv', '.tsv', '.txt', '.csv.gz', '.tsv.gz', '.txt.gz')): return 'text_candidate'
    if name.endswith(('.tar', '.tar.gz', '.tgz', '.zip')): return 'archive'
    if name.endswith(('.fastq', '.fastq.gz', '.sra')): return 'raw_reads'
    return 'unknown'

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
            cellmap = local(root, row['cell_map_path'], area + '/.workflow')
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
