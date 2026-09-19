"""Lightweight pre-build validation for a physical 10x Matrix Market trio."""
from __future__ import annotations

import gzip
from collections import Counter

from common import local


def _open_text(path):
    return gzip.open(path, 'rt', encoding='utf-8-sig', newline='') if path.name.lower().endswith('.gz') \
        else path.open('rt', encoding='utf-8-sig', newline='')


def _one(directory, names, label):
    candidates = [directory / suffix for name in names for suffix in (name, name + '.gz')]
    hits = [path for path in candidates if path.is_file()]
    if len(hits) != 1:
        raise ValueError(f'Missing or ambiguous 10x {label}: expected exactly one of {", ".join(names)} (plain or .gz)')
    return hits[0]


def _matrix_dimensions(path):
    with _open_text(path) as handle:
        banner = handle.readline().strip()
        if not banner.lower().startswith('%%matrixmarket matrix coordinate'):
            raise ValueError('10x matrix.mtx is not Matrix Market coordinate format')
        for line in handle:
            line = line.strip()
            if not line or line.startswith('%'):
                continue
            fields = line.split()
            if len(fields) != 3:
                break
            try:
                dimensions = tuple(int(value) for value in fields)
            except ValueError:
                break
            if any(value < 0 for value in dimensions):
                break
            return dimensions
    raise ValueError('10x matrix.mtx has no valid dimension line')


def _blank(value):
    return value.strip().lower() in ('', 'na', 'nan', 'null', 'none')


def probe_10x_trio(row, root):
    directory = local(root, row['local_path'], f"data/{row['database']}/raw")
    if not directory.is_dir():
        raise ValueError('10x trio local_path must be a directory')
    matrix_path = _one(directory, ('matrix.mtx',), 'matrix')
    features_path = _one(directory, ('features.tsv', 'genes.tsv'), 'features')
    barcodes_path = _one(directory, ('barcodes.tsv',), 'barcodes')
    matrix_rows, matrix_columns, matrix_nonzero = _matrix_dimensions(matrix_path)

    features = []
    with _open_text(features_path) as handle:
        for line in handle:
            if line.strip():
                fields = line.rstrip('\r\n').split('\t')
                if len(fields) < 1 or len(fields) > 3:
                    raise ValueError('10x features.tsv must contain one to three tab-delimited columns')
                features.append(fields)
    with _open_text(barcodes_path) as handle:
        barcodes = [line.rstrip('\r\n').split('\t')[0] for line in handle if line.strip()]

    if matrix_rows != len(features) or matrix_columns != len(barcodes):
        raise ValueError(f'10x dimensions mismatch: matrix={matrix_rows}x{matrix_columns}, '
                         f'features={len(features)}, barcodes={len(barcodes)}')
    if any(_blank(value) for value in barcodes) or len(set(barcodes)) != len(barcodes):
        raise ValueError('10x barcodes must be nonblank, non-NA, and unique')
    feature_ids = [fields[0] for fields in features]
    if any(_blank(value) for value in feature_ids) or len(set(feature_ids)) != len(feature_ids):
        raise ValueError('10x feature IDs must be nonblank, non-NA, and unique')

    gene_names = [fields[1] if len(fields) >= 2 else '' for fields in features]
    missing_gene_names = sum(_blank(value) for value in gene_names)
    present_gene_names = [value for value in gene_names if not _blank(value)]
    duplicate_gene_names = len(present_gene_names) - len(set(present_gene_names))
    fallback = 'feature_id_for_missing_gene_name' if missing_gene_names else 'none'

    feature_type_counts = {}
    if any(len(fields) >= 3 for fields in features):
        if any(len(fields) < 3 or _blank(fields[2]) for fields in features):
            raise ValueError('10x feature type column is partially missing or NA')
        feature_type_counts = dict(Counter(fields[2] for fields in features))
        if 'Gene Expression' not in feature_type_counts:
            raise ValueError('10x feature types do not contain Gene Expression')

    return {
        'local_path': row['local_path'],
        'matrix_path': matrix_path.name,
        'features_path': features_path.name,
        'barcodes_path': barcodes_path.name,
        'matrix_rows': matrix_rows,
        'matrix_columns': matrix_columns,
        'matrix_nonzero': matrix_nonzero,
        'feature_rows': len(features),
        'barcode_rows': len(barcodes),
        'barcode_unique': True,
        'feature_id_unique': True,
        'missing_gene_names': missing_gene_names,
        'duplicate_gene_names': duplicate_gene_names,
        'feature_name_fallback': fallback,
        'duplicate_gene_name_strategy': 'Seurat::Read10X make.unique' if duplicate_gene_names else 'none',
        'feature_type_counts': feature_type_counts,
        'gene_expression_features': feature_type_counts.get('Gene Expression', len(features)),
    }
