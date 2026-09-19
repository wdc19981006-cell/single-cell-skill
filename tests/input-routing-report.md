# Stage B input-routing validation report

Validation date: 2026-09-19 (original real-data metrics retained from 2026-09-16). Runtime: Windows 11, R 4.5.3, Seurat 5.5.1, SeuratObject 5.4.0, Python 3.12. Tests used generated fixtures or previously recorded real-run evidence; the 2026-09-19 change did not redownload or rebuild GSE166504, GSE189175, or GSE289173. Runtime data remains under git-ignored `data/`.

## Reader coverage

| Route / case | Reader or behavior | Fixture / validation | Result |
|---|---|---|---|
| 10x MTX, plain/gzip trio | `Seurat::Read10X()` | 250 × 4 synthetic trio plus GSE231993 12-sample rebuild | PASS |
| 10x H5 | schema validation + `Seurat::Read10X_h5()` | Valid synthetic H5; ordinary non-10x H5 rejected | PASS |
| H5AD raw counts | explicit `counts` and `raw:X`; Python anndata fallback in this runtime | `.h5ad` and `.h5ad.gz`, author obs/PCA/UMAP present | PASS |
| H5AD normalized-only selection | selected normalized `X` | Noninteger/log-normalized matrix rejected | PASS |
| Seurat RDS | `readRDS()` + exact RNA counts extraction | Author object contains normalized data, PCA, UMAP and author cluster | PASS |
| Seurat RDS split layers | exact counts-layer rule | `counts.A`/`counts.B` fixture rejected | PASS |
| CSV | `data.table::fread()` | genes × cells with one dropped annotation column | PASS |
| TSV | `data.table::fread()` | genes × cells | PASS |
| TXT / TXT.GZ | `data.table::fread()` | genes × cells, plain and gzip | PASS |
| cells × genes | `data.table::fread()` + explicit transpose | transposed TSV | PASS |
| pooled/shared input | one reader result + one exhaustive cell map | Three samples share six cells | PASS: expression calls 1, cell-map calls 1 |
| expression + separate metadata | underlying reader + configured metadata keys | `barcode` → `author_sample`, extra author column ignored | PASS |
| archive exact member | transport extraction only | synthetic TAR exact member and traversal rejection | PASS |
| FASTQ/SRA | unsupported V1 route | FASTQ, FASTQ.GZ and SRA candidates | PASS: stopped before quantification |
| global `min.cells` | one merged `CreateSeuratObject()` | gene in one cell of each of three samples retained; two-cell gene removed | PASS |
| clean final object | one RNA counts layer; no author analysis state | serialized synthetic E2E plus independent validator | PASS |
| small CSV reader selection | prebuild dimensions + RAM decision, then `data.table::fread()` | 2 × 2 genes-by-cells fixture | PASS |
| small compressed file with unsafe dense dimensions | Python streaming selected from `genes * cells * 8`, independent of gzip bytes | GSE166504-shaped 25,127 × 82,168 decision fixture, 32 GiB RAM | PASS |
| 10x missing gene name | prebuild records feature-ID fallback | unique feature IDs plus one `NA` gene name | PASS |
| 10x dimension mismatch | prebuild STOP before `Read10X()` | matrix 3 × 2 versus two feature rows | PASS |
| non-single-cell modality gate | STOP before download/build | GSE116504-shaped `Expression profiling by array` Stage A fixture | PASS |
| multi-input memory lifecycle | release prepared/aligned lists after combine and merged reference after construction | exact full-count comparison in global min.cells regression and synthetic E2E | PASS |

The synthetic pooled test increments injected reader functions directly. It does not infer call counts from log text.

## GSE181919 pooled performance rebuild

Global GEO MCP verification returned `complete=true`, 37 samples and two GSE-level text files: `GSE181919_UMI_counts.txt.gz` plus `GSE181919_Barcode_metadata.txt.gz`. The existing confirmed manifest and local expression file were reused; no expression download occurred. The previous final RDS was archived before rebuilding.

| Metric | Observed |
|---|---:|
| Manifest rows | 37 |
| Unique input signatures | 1 |
| Shared expression file size | 127,878,601 bytes |
| Reader | `data.table::fread` |
| Expression matrices actually read | 1 |
| Cell maps actually read | 1 |
| Input features | 20,000 |
| Input/final cells | 54,239 |
| Read seconds | 44.280 |
| Mapping seconds | 0.220 |
| Feature alignment seconds | 3.920 |
| Matrix combine seconds | 3.870 |
| CreateSeuratObject seconds | 9.030 |
| Metadata mapping seconds | 0.040 |
| Validation seconds | 7.780 |
| Serialization seconds | 41.060 |
| Total build seconds | 112.630 |

Result: PASS. The former manifest-row execution model would have requested 37 reads of the pooled matrix; the router performed one physical expression read and one cell-map read while keeping the pooled sparse matrix intact.

## GSE231993 10x non-regression rebuild

The local project initially lacked GSE231993. The published v2 Release asset was downloaded (378,435,397 bytes; SHA256 `cb5da8d09d706607a629d119629a182411376b4d18b811456de7f7fabe0e6493`), its 12 existing 10x trios and confirmed workflow were reused, and its prior final RDS was archived before rebuilding.

| Metric | Expected | Observed | Result |
|---|---:|---:|---|
| Cells | 60,665 | 60,665 | PASS |
| Genes | 25,953 | 25,953 | PASS |
| Total UMI | 305,606,810 | 305,606,810 | PASS |
| Nonzero entries | 58,767,679 | 58,767,679 | PASS |
| Expression matrices read | 12 | 12 | PASS |
| Cell maps read | 0 | 0 | PASS |

Per-sample cells: GSM7307094 4,481; GSM7307095 5,993; GSM7307096 3,866; GSM7307097 6,421; GSM7307098 776; GSM7307099 6,529; GSM7307100 4,596; GSM7307101 5,305; GSM7307102 5,270; GSM7307103 4,633; GSM7307104 6,578; GSM7307105 6,217. All match the required regression values. Total build time was 128.580 seconds; `CreateSeuratObject()` took 6.090 seconds and serialization took 30.690 seconds.

## 2026-09-19 regression commands

- Python unit tests: 67/67 PASS.
- R reader regression (`test_seurat.R`): PASS.
- R pooled call-count regression (`test_pooled_inputs.R`): PASS.
- R global min.cells regression (`test_global_min_cells.R`): PASS.
- R text reader-selection regression (`test_stream_text.R`): PASS.
- Synthetic subprocess E2E, serialized validation, pooled E2E and H5AD failure/retry: PASS.
- No 2026-09-19 real large-data rebuild was performed, as requested.
