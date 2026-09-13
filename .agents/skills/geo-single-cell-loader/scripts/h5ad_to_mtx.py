"""Export one explicitly selected AnnData count matrix for R ingestion."""
import argparse
import gzip
import json
import shutil
import tempfile
from pathlib import Path

def assert_ids(values, label):
    values = [str(x) for x in values]
    if not values or len(set(values)) != len(values) or any(not x or '\t' in x or '\n' in x or '\r' in x for x in values):
        raise ValueError(f'Missing, duplicate, or unsafe {label} IDs')
    return values

def export(source, count_source, output):
    try:
        import anndata
        import numpy as np
        import scipy.io
        import scipy.sparse
    except ImportError as e:
        raise RuntimeError('H5AD fallback requires anndata, numpy, and scipy in GEO_SINGLE_CELL_PYTHON') from e
    temporary = None
    path = source
    if str(source).lower().endswith('.gz'):
        handle = tempfile.NamedTemporaryFile(suffix='.h5ad', delete=False)
        temporary = Path(handle.name); handle.close()
        with gzip.open(source, 'rb') as src, temporary.open('wb') as dst: shutil.copyfileobj(src, dst)
        path = temporary
    try:
        adata = anndata.read_h5ad(path)
        if count_source == 'X': matrix = adata.X
        elif count_source.startswith('raw:'):
            requested = count_source[4:]
            if requested != 'X' or adata.raw is None: raise ValueError('Requested AnnData raw:X is absent')
            matrix = adata.raw.X
            genes = assert_ids(adata.raw.var_names, 'raw feature')
        else:
            if count_source not in adata.layers: raise ValueError('Requested AnnData layer absent: ' + count_source)
            matrix = adata.layers[count_source]
        if not count_source.startswith('raw:'): genes = assert_ids(adata.var_names, 'feature')
        cells = assert_ids(adata.obs_names, 'cell')
        matrix = scipy.sparse.coo_matrix(matrix)
        values = matrix.data
        if not np.issubdtype(values.dtype, np.number) or not np.all(np.isfinite(values)) or np.any(values < 0) or np.any(np.abs(values - np.rint(values)) > 1e-8):
            raise ValueError('Selected H5AD matrix is not finite nonnegative integer counts')
        output.mkdir(parents=True, exist_ok=False)
        scipy.io.mmwrite(output / 'matrix.mtx', matrix.transpose().tocsr(), field='integer')
        (output / 'features.tsv').write_text('\n'.join(genes) + '\n', encoding='utf-8')
        (output / 'barcodes.tsv').write_text('\n'.join(cells) + '\n', encoding='utf-8')
        (output / 'selection.json').write_text(json.dumps({'reader': 'anndata', 'count_source': count_source, 'genes': len(genes), 'cells': len(cells)}, indent=2), encoding='utf-8')
    finally:
        if temporary is not None: temporary.unlink(missing_ok=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path); p.add_argument('count_source'); p.add_argument('output', type=Path)
    a = p.parse_args(); export(a.source, a.count_source, a.output)
