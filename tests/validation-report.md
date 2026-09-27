# Validation report — 2026-09-27

This report records checks actually run for the targeted QC and raw-handoff fixes. No real GSE was rebuilt or downloaded.

## Environment

| Component | Verified value |
|---|---|
| R | 4.5.0 (2025-04-11 ucrt) |
| R library | `D:/R/R-4.5.0/library` |
| Seurat | 5.4.0 |
| DoubletFinder | 2.0.6 |
| Rscript exit-139 workaround | `qc/r450_rscript.sh` (used for every R check) |
| Python | `C:/Python312/python.exe` |

The R launcher check passed: normal exit remained successful, an intentional R error returned exit 1, and `--vanilla` was rejected because it disables the required exit profile. Startup locale warnings and warnings about packages built under R 4.5.3 were observed; no package was installed, upgraded, or removed.

## Results

| Check | Actual result |
|---|---|
| Python `unittest discover -s tests -p 'test_*.py' -q` | **72/72 passed**, no skips. Includes raw-handoff failure preserving `COMPLETE_RAW` and successful validation. The streaming-route test now uses the R 4.5.0 launcher rather than looking for R 4.5.3. |
| `tests/test_qc.R` | **PASS**: strict QC boundaries, precheck/formal-filter equivalence, overall and per-sample >15% guard, persisted user acceptance, precheck cell-cycle gene coverage, PC rule, 57-sample legend, UMAP/three-panel cell-cycle plots, clean final counts-only object, and existing QC decision/DoubletFinder mock flow. |
| `tests/test_qc_real_df.R` | **PASS**: actual DoubletFinder `paramSweep`, pK selection and both classification calls ran on a synthetic 130-cell object. This verifies real DoubletFinder computation, not a real GEO cohort. |
| Other R tests | **7/7 passed**: `test_10x_literal_na.R`, `test_final_object.R`, `test_global_min_cells.R`, `test_pooled_inputs.R`, `test_seurat.R`, `test_stream_text.R`, `test_utf8_manifest.R`. H5AD fallback was explicitly set to `C:/Python312/python.exe` for `test_seurat.R`. |
| R launcher shell test | **PASS**: `tests/test_r450_launcher.sh`. |
| Synthetic E2E | **PASS** in an isolated ignored-data directory: four-input raw build + independent validation (16 cells/250 genes), pooled matrix/cell map, and normalized-H5AD failure/retry. Both temporary E2E directories were removed after testing; existing GSE data was untouched. |
| `git diff --check` | Run before commit; see commit verification. |

The first `test_seurat.R` run without `GEO_SINGLE_CELL_PYTHON` stopped because the optional zellkonverter route was unavailable. After explicitly selecting the existing Python fallback, an existing synthetic `split_input.rds` could not be overwritten on this host. The test now writes that negative-case fixture to a unique temporary filename; the rerun passed. The first isolated E2E attempt lacked the source file needed by `test_final_object.R`; adding it to a fresh isolated fixture allowed the full E2E rerun to pass. These were test setup issues, not passes claimed from failed runs.

No real GSE QC or real-data DoubletFinder run was performed in this revision.
