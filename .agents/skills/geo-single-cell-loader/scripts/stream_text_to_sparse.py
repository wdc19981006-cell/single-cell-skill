"""Convert a large GEO genes-by-cells tabular gzip into sparse triplets.

The source is read once. Only nonzero entries are written to temporary binary
files; the input remains untouched. The header has cell IDs from column one,
while each following row starts with a gene ID.
"""

import gzip
import itertools
import json
import sys
import warnings
from pathlib import Path

import numpy as np


def convert(source: Path, output: Path) -> dict:
    output.mkdir(parents=True, exist_ok=False)
    opener = gzip.open if source.name.lower().endswith(".gz") else open
    with opener(source, "rb") as src, (output / "i.bin").open("wb") as rows, \
            (output / "j.bin").open("wb") as cols, \
            (output / "x.bin").open("wb") as counts, \
            (output / "genes.tsv").open("w", encoding="utf-8", newline="\n") as genes_file:
        header = src.readline().rstrip(b"\r\n")
        first_row = src.readline()
        if not first_row:
            raise ValueError("Empty matrix")
        _, separator, first_payload = first_row.rstrip(b"\r\n").partition(b"\t")
        if not separator:
            raise ValueError("Matrix row 1 lacks expression fields")
        header_fields = [item.decode("utf-8-sig") for item in header.split(b"\t")]
        value_columns = first_payload.count(b"\t") + 1
        if len(header_fields) == value_columns + 1:
            cells = header_fields[1:]
        elif len(header_fields) == value_columns:
            cells = header_fields
        else:
            raise ValueError("Header and expression column counts differ")
        if not cells or any(not cell for cell in cells) or len(set(cells)) != len(cells):
            raise ValueError("Blank or duplicate matrix cell IDs")
        if len(cells) >= 2**31:
            raise ValueError("Matrix has too many cells for dgCMatrix")
        (output / "cells.tsv").write_text("\n".join(cells) + "\n", encoding="utf-8")
        n_cells = len(cells)
        seen_genes = set()
        n_genes = 0
        nnz = 0
        warnings.simplefilter("error", DeprecationWarning)
        for line in itertools.chain((first_row,), src):
            line = line.rstrip(b"\r\n")
            gene_bytes, sep, payload = line.partition(b"\t")
            if not sep:
                raise ValueError(f"Matrix row {n_genes + 1} lacks expression fields")
            gene = gene_bytes.decode("utf-8")
            if not gene or gene in seen_genes:
                raise ValueError(f"Blank or duplicate gene ID at row {n_genes + 1}")
            seen_genes.add(gene)
            try:
                values = np.fromstring(payload, dtype=np.float64, sep="\t")
            except (ValueError, DeprecationWarning) as exc:
                raise ValueError(f"Malformed count at row {n_genes + 1}: {exc}") from exc
            if values.size != n_cells:
                raise ValueError(f"Matrix row {n_genes + 1} has {values.size} values; expected {n_cells}")
            if not np.all(np.isfinite(values)) or np.any(values < 0) or np.any(values != np.floor(values)):
                raise ValueError(f"Noninteger, negative, or nonfinite count at row {n_genes + 1}")
            positions = np.flatnonzero(values)
            nonzero = positions.size
            if nnz + nonzero >= 2**31:
                raise ValueError("Too many nonzero entries for dgCMatrix")
            np.full(nonzero, n_genes + 1, dtype="<i4").tofile(rows)
            (positions + 1).astype("<i4").tofile(cols)
            values[positions].astype("<f8").tofile(counts)
            genes_file.write(gene + "\n")
            n_genes += 1
            nnz += nonzero
    if not n_genes or not nnz:
        raise ValueError("Empty matrix")
    result = {"genes": n_genes, "cells": n_cells, "nonzero": nnz,
              "reader": "NumPy streaming tabular counts"}
    (output / "dimensions.json").write_text(json.dumps(result), encoding="utf-8")
    print(json.dumps(result), flush=True)
    return result


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("Usage: stream_text_to_sparse.py INPUT.txt[.gz] OUTPUT_DIRECTORY")
    convert(Path(sys.argv[1]), Path(sys.argv[2]))
