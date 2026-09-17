"""Small pre-build sample of a text count matrix; never loads the full matrix."""
import csv
import gzip
import re
from decimal import Decimal, InvalidOperation

from common import TRANSFORMED_EXPRESSION, local, missing

MAX_LINE = 8 * 1024 * 1024
GENE_COLUMNS = {'gene', 'genes', 'geneid', 'gene_id', 'gene_symbol', 'symbol',
                'feature', 'features', 'feature_id', 'ensembl_id'}
CELL_COLUMNS = {'cell', 'cell_id', 'barcode', 'barcodes'}
BARCODE = re.compile(r'^[ACGTN]{12,}(?:[-_]\d+)?$', re.I)
GENE_ID = re.compile(r'^(?:ENSG|ENSMUSG|ENSMUST)\d+', re.I)
DELIMITERS = {'comma': ',', 'tab': '\t', 'space': ' '}


def _read_lines(path, count=9):
    opener = gzip.open if path.name.lower().endswith('.gz') else open
    lines = []
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as handle:
        for _ in range(count):
            line = handle.readline(MAX_LINE + 1)
            if not line:
                break
            if len(line) > MAX_LINE:
                raise ValueError('Text schema probe line exceeds the small-read limit')
            if line.strip():
                lines.append(line)
    if len(lines) < 2:
        raise ValueError('Text schema probe needs a header and an expression row')
    return lines


def _parse(lines, name):
    delimiter = DELIMITERS[name]
    if name == 'space':
        rows = [re.split(r' +', line.rstrip('\r\n')) for line in lines]
    else:
        try:
            rows = [next(csv.reader([line], delimiter=delimiter)) for line in lines]
        except csv.Error as error:
            raise ValueError(f'Text schema probe cannot parse {name}: {error}') from error
    width = len(rows[0])
    body_widths = {len(row) for row in rows[1:]}
    if width < 2 or len(body_widths) != 1 or next(iter(body_widths)) not in (width, width + 1):
        raise ValueError(f'Text schema probe has inconsistent columns for {name}')
    return rows, next(iter(body_widths)) == width + 1


def _choose_delimiter(lines, configured):
    if not missing(configured):
        if configured not in DELIMITERS:
            raise ValueError('Unsupported explicit text delimiter')
        rows, header_missing_id = _parse(lines, configured)
        return configured, rows, header_missing_id
    candidates = []
    for name in DELIMITERS:
        try:
            rows, header_missing_id = _parse(lines, name)
            candidates.append((len(rows[0]), name, rows, header_missing_id))
        except ValueError:
            continue
    if not candidates:
        raise ValueError('Text delimiter cannot be identified from sampled rows')
    candidates.sort(reverse=True)
    if len(candidates) > 1 and candidates[0][0] == candidates[1][0]:
        raise ValueError('Ambiguous text delimiter; configure it explicitly')
    return candidates[0][1], candidates[0][2], candidates[0][3]


def _orientation(header, rows, header_missing_id):
    if header_missing_id:
        return 'genes_by_cells'
    first = header[0].strip().lower()
    row_barcodes = sum(bool(BARCODE.fullmatch(row[0])) for row in rows) / len(rows)
    gene_headers = sum(bool(GENE_ID.match(value)) for value in header[1:]) / max(len(header) - 1, 1)
    if row_barcodes >= 0.8 and gene_headers >= 0.5:
        return 'cells_by_genes'
    if first == '' or first in GENE_COLUMNS:
        return 'genes_by_cells'
    if first in CELL_COLUMNS:
        return 'cells_by_genes'
    cell_headers = sum(bool(BARCODE.fullmatch(value)) for value in header[1:]) / max(len(header) - 1, 1)
    if cell_headers >= 0.8 and row_barcodes < 0.5:
        return 'genes_by_cells'
    return None


def probe_text(row, root):
    """Return only reliable technical fields and observed evidence."""
    area = f"data/{row['database']}/raw"
    path = local(root, row['local_path'], area)
    if not path.is_file():
        raise ValueError('Text schema probe input is missing: ' + row['local_path'])
    if TRANSFORMED_EXPRESSION.search(path.name):
        raise ValueError('Normalized/processed text is not raw counts: ' + path.name)
    lines = _read_lines(path)
    delimiter, parsed, header_missing_id = _choose_delimiter(lines, row.get('delimiter'))
    header, values = parsed[0], parsed[1:]
    inferred = _orientation(header, values, header_missing_id)
    configured_orientation = row.get('orientation')
    if not missing(configured_orientation) and inferred and configured_orientation != inferred:
        raise ValueError(f'Explicit text orientation conflicts with sampled {inferred} layout')
    orientation = configured_orientation if not missing(configured_orientation) else inferred
    if orientation not in ('genes_by_cells', 'cells_by_genes'):
        raise ValueError('Text orientation is ambiguous; inspect and configure it explicitly')
    observed_id = '__row_names__' if header_missing_id or not header[0] else header[0]
    configured_id = row.get('feature_column')
    if (not missing(configured_id) and configured_id != observed_id and
            not (configured_id == '__row_names__' and orientation == 'genes_by_cells')):
        raise ValueError(f'Explicit feature_column {configured_id!r} differs from sampled ID column {observed_id!r}')
    column_headers = header if header_missing_id else header[1:]
    if any(not value for value in column_headers) or len(set(column_headers)) != len(column_headers):
        raise ValueError('Text cell/feature header is blank or duplicated')
    if any(not line[0] for line in values):
        raise ValueError('Text row ID is blank')
    if orientation == 'cells_by_genes' and len({line[0] for line in values}) != len(values):
        raise ValueError('Sampled cell IDs are duplicated')
    dropped = [value for value in (row.get('drop_columns') or '').split(';') if value]
    if any(value not in column_headers for value in dropped):
        raise ValueError('Explicit drop_columns are absent from text header')
    expression_indexes = [i for i in range(1, len(values[0]))
                          if column_headers[i - 1] not in dropped]
    if not expression_indexes:
        raise ValueError('No expression columns remain after drop_columns')
    for line in values:
        for index in expression_indexes:
            value = line[index].strip()
            try:
                number = Decimal(value)
            except InvalidOperation:
                raise ValueError(f'Nonexpression text value in column {column_headers[index - 1]!r}') from None
            if not number.is_finite() or number < 0 or number != number.to_integral_value():
                raise ValueError('Normalized/processed or invalid text values; sampled expression is not nonnegative integer counts')
    source_bytes = path.stat().st_size
    streaming_threshold = 128 * 1024**2 if delimiter == 'comma' else 256 * 1024**2
    inferred_id = ('__row_names__' if (orientation == 'genes_by_cells' and
                    not dropped and delimiter in ('comma', 'tab') and
                    source_bytes >= streaming_threshold) else observed_id)
    return {
        'local_path': row['local_path'], 'delimiter': delimiter, 'orientation': orientation,
        'feature_column': configured_id if not missing(configured_id) else inferred_id,
        'observed_id_header': observed_id, 'cell_header_unique': True,
        'header_missing_id': header_missing_id,
        'sampled_rows': len(values), 'sampled_columns': len(expression_indexes),
        'sampled_nonnegative_integer_counts': True,
        'full_matrix_validation_required': True, 'source_bytes': source_bytes,
    }
