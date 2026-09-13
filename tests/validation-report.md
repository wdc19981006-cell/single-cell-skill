# Validation report

Environment and date: Windows 11, 2026-09-13. Python 3.12.10. Final R tests use an isolated R 4.5.3 library with Seurat 5.5.1. Exact dependency versions are captured by `check_dependencies.R` and the synthetic run summary.

## Executed local checks

- `skill-creator/scripts/quick_validate.py` with `PYTHONUTF8=1`: passed. The variable is necessary because the validator otherwise uses this Chinese Windows session's GBK default when reading the UTF-8 Skill file.
- Python `compileall`: passed for Skill and test scripts.
- Python unittest: 15 tests passed. These cover duplicate/missing groups, exact confirmation coverage and change detection, unsorted sample/group mapping, optional-field evidence, path escape, shared-matrix cell map requirement and hash invalidation, H5AD-vs-H5 routing, filtered-over-raw preference, trio detection, safe archive member extraction and SOFT fact preservation.
- R parse: all five R files parsed successfully.
- Reader regression test: passed for actual Read10X trio, Read10X_h5, H5AD.gz selected `counts`, H5AD `raw:X`, genes-by-cells CSV, cells-by-genes TSV and Seurat RDS reconstruction. A deliberately normalized H5AD X was rejected. H5AD sent to the 10x H5 route, a non-Seurat RDS, wrong metadata, changed sample keys and metadata row reordering were rejected.
- Synthetic end-to-end build: four samples with colliding original barcodes were prefixed, independently created, merged, JoinLayers applied, mapped from an unsorted mock manifest, serialized and reloaded. Independent validation reported 250 genes, 16 cells, four samples, two groups and exact metadata/Cells alignment. Reader tests, build, independent validation and dependency check all exited with status 0 under R 4.5.3.

The synthetic fixture uses the reserved, non-GEO identifier `GSE999999999` and a confirmation statement explicitly marked simulated. Binary fixtures, generated manifests, logs and outputs are ignored by Git.

## H5AD compatibility finding

Bioconductor 3.22 `zellkonverter` 1.20.1 was attempted on installed R 4.5.0 and 4.5.3. Its Windows `configure.win` nested Rscript process segfaulted on this host; skipping configure then failed package lazy loading, so it is not reported as installed or tested successfully. Existing dependencies `SingleCellExperiment` 1.32.0, `rhdf5` 2.54.1 and `hdf5r` 1.3.12 are present.

The implemented reader still prefers `zellkonverter::readH5AD(reader="R")` when it is loadable. The executed fallback used explicit `GEO_SINGLE_CELL_PYTHON=C:/Python312/python.exe` with anndata 0.13.3.post0, numpy 2.3.2 and scipy 1.16.1. It exports only the explicitly selected layer/raw matrix and IDs to a temporary Matrix Market file, validates integer counts, and records the reader in `run_summary.txt`. It never routes H5AD through Read10X_h5.

An independent Windows exit issue was diagnosed during validation. This host process omitted `PROCESSOR_ARCHITECTURE`; cli 3.6.6 dereferences that variable during unload and R exited with `0xC0000005` after otherwise successful commands. This matches [r-lib/cli #375](https://github.com/r-lib/cli/issues/375) and the unchecked environment access described in [PR #838](https://github.com/r-lib/cli/pull/838). The scripts now set this normally-present value only inside the current R process when absent, based on `R.version$arch`. A direct cli/Seurat probe then exited 0. Final tests were rerun with this fix and are not accepted based only on printed success text.

## GEO discovery checks

NCBI GEO SOFT metadata and PubMed abstract retrieval were executed without bulk expression downloads. `--max-samples 2` intentionally sampled only two GSM records from each study, so these are format/discovery checks rather than complete biological sample reports.

| Dataset | Observed result | Status |
|---|---|---|
| GSE231993 | 12 GSM; GSM7307094/95 have per-sample gzipped matrix, features and barcodes; Sample data processing says Cell Ranger/10x | Passed 10x trio discovery |
| GSE202051 | 74 GSM; Series publishes three `.h5ad.gz` files including nucleus and organoid collections | Passed H5AD routing discovery; counts layer remains an inspection requirement |
| GSE211644 | 50 GSM; fresh and grown matrix/genes/barcodes are accompanied by separate metadata CSVs | Passed split expression/metadata discovery; explicit cell map required |
| GSE229413 | 32 GSM; GSM7162998/99 each publish raw and filtered feature_bc_matrix H5 files | Passed 10x H5 and same-GSM filtered precedence discovery |

No full large dataset was downloaded. A one-sample real GSE231993 end-to-end run remains gated by explicit user group confirmation, as required by the Skill. It must not be replaced with an inferred or test-only group without that confirmation.
