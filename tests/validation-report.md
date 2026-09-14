# Validation report

Executed on 2026-09-14, Windows 11, for the dataset layout and natural-language Skill workflow refactor. Only the checks explicitly described below were executed in this revision.

## Environment actually used

| Component | Version/result |
|---|---|
| Python | 3.12.0, 64-bit, `C:/Python312/python.exe` |
| R | 4.5.3 (2026-03-11 ucrt) |
| Seurat / SeuratObject | 5.5.1 / 5.4.0 |
| Matrix / jsonlite / hdf5r | 1.7.4 / 2.0.0 / 1.3.16 |
| anndata / numpy / scipy / h5py | 0.13.3.post0 / 2.3.2 / 1.16.1 / 3.16.0 |
| zellkonverter / SingleCellExperiment / SummarizedExperiment | Not loadable in the R environment used for this run |

H5AD checks used the explicit `GEO_SINGLE_CELL_PYTHON` anndata fallback. The native zellkonverter route was preserved but was not successfully exercised here. R subprocesses used `LC_ALL=C`; the existing process-local Windows architecture workaround remains in place. Dependencies were inspected, not installed or changed.

## Python and Skill checks

- `unittest discover -s tests -p 'test_*.py' -v`: **28 tests passed**, exit 0.
- `compileall` for Skill Python scripts and tests: passed. The first attempt could not write bytecode under the protected `.agents` directory; it was rerun successfully with `PYTHONPYCACHEPREFIX` pointing into ignored `data/GSE999999999/.workflow/pycache`.
- `skill-creator/scripts/quick_validate.py`: **Skill is valid**, exit 0, with `PYTHONUTF8=1`.
- Git whitespace/error check: `git diff --check` passed.

Behavioral coverage includes: all required metadata/group fields; duplicate samples; optional evidence; path confinement and rejection of absolute/parent/cross-GSE paths; H5AD/H5 routing; filtered preference/ambiguity; gzip/plain trio detection; exact user group coverage; user identity marker and literal statement; manifest hash invalidation; cell-map hash invalidation; explicit author sample-to-GSM ownership; archive member path safety; full GSM traversal without bulk downloads; waiting TXT; confirmed groups and UTF-8 TXT; omitted absent optional fields; development probes isolated from complete reports and rejected by confirmation; download to raw; centralized provenance; direct and archive-member reuse; changed URL/checksum rejection; unverified-file refusal; interrupted downloads persisting and reusing completed files.

Implicit invocation remains enabled. Chinese and English natural-language entry examples are present in Skill discovery metadata and instructions. This validates configuration and workflow implementation, not a live independent Codex skill-selection session.

## R regression checks

- Parsed all five Skill R files plus both R test files: **7 files parsed**.
- `tests/test_seurat.R`: passed, process exit 0, after the final R path changes.
- Actual readers exercised: 10x trio, 10x H5, uncompressed H5AD, H5AD.gz counts layer, H5AD raw:X, genes-by-cells CSV, cells-by-genes TSV, and Seurat RDS counts.
- Expected failures verified: normalized H5AD X, H5AD on the H5 route, non-Seurat RDS, missing required metadata, wrong group mapping, changed sample keys, missing barcode prefixes and reordered metadata rows. Expected error output from these negative cases is not a test failure.

## Synthetic end-to-end

`tests/run_end_to_end.py --rscript <verified Rscript>` completed with exit 0 after final code changes. All group statements are marked SIMULATED; no real clinical grouping is implied.

| Case | Observed result |
|---|---|
| Main `GSE999999999` fixture | Four independent inputs, **4 samples, 16 cells, 250 genes**. Real R build, serialization, reload, separate validate script and independent object assertions passed. |
| Stable metadata | Unsorted sample/group manifest matched by keys. Four cells in Fixture_A, twelve in Fixture_B; original colliding barcodes were sample-prefixed; metadata row names exactly equal Cells. |
| Final visible TXT | COMPLETE; sample/group sections; totals; output paths; reader/count sources; versions; construction thresholds and warnings present. Entirely absent optional fields omitted. |
| Existing RDS | Second build refused. RDS SHA256 and completed TXT remained unchanged. |
| Pooled matrix | Explicit reversed-order cell map, **2 samples, 8 cells**; build and independent validation passed. Altered cell-map checksum was rejected. |
| Failure and retry | Deliberately normalized H5AD failed; TXT became BUILD_FAILED with a reason. Raw SHA256 unchanged and no final RDS written. Explicit simulated reconfirmation of counts layer then succeeded and removed stale failure text. |

Retained local result:

```text
data/GSE999999999/
├── raw/
├── seurat_raw.rds
├── sample_info.txt
└── .workflow/
    ├── sample_report.csv
    ├── group_confirmation.json
    ├── sample_manifest.csv
    ├── sample_manifest.confirmation.json
    ├── download.json
    ├── run_summary.txt
    └── integration-tests.txt
```

The synthetic sources were generated locally, so `download.json` is an empty plan. Real downloader IO/reuse/checksum behavior is covered separately by mocked HTTPS responses in the 28 offline tests. Additional pooled/failure fixtures used system temporary directories and were cleaned. No new root datasets/manifests/logs/output directory was created.

## Real GEO metadata-only discovery

[GSE231993](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231993) was rerun without `--max-samples`, using the existing HTTPS SOFT fallback because bio-server GEO MCP tools were not exposed in this session. The initial sandbox network denial was resolved through an authorized network escalation.

- **12/12 GSM records inspected**, complete=true; no partial inventory accepted.
- **37 processed-file candidate URLs**: 36 per-sample trio member links and one Series RAW.tar link; **12 complete candidate 10x input plans** under the new raw layout.
- Linked PubMed abstract **37344477** retrieved. Public GSM characteristics and processing descriptions were reviewed: colon tissue, ulcerative colitis/healthy donor background, explicit individual IDs, and Cell Ranger processing. Original source/inflammation and author sample labels remain in sample descriptions.
- `data/GSE231993/.workflow/inspection.json`, `sample_report.csv` and UTF-8 `data/GSE231993/sample_info.txt` generated. TXT contains all twelve samples, reviewed metadata and **WAITING_FOR_GROUP_CONFIRMATION**.
- Verified no real `raw/`, RDS, group confirmation, manifest or download log was created. No bulk expression matrix was downloaded; counts content and clinical downstream outcomes were not tested.

GSE202051, GSE211644 and GSE229413 were **not rerun in this revision**. Their 2026-09-13 two-GSM development records were preserved during migration; those historical partial checks are not reported as current full-discovery passes. Additional full-text publication supplements were not exhaustively reviewed during this metadata-only regression.

## Existing local data migration and repository contents

The old folders contained **26 actual generated files**, besides `.gitkeep` placeholders. They were copied into the corresponding `data/<GSE>/.workflow/legacy/` archives, checked for equal byte sizes and SHA256, then original copies removed. A final independent pass verified all 26 archived files against their recorded hashes. Per-GSE migration counts: GSE999999999 17, GSE231993 3, and GSE202051/GSE211644/GSE229413 2 each.

The old files consisted of synthetic fixtures/results and historical metadata/manifest records; there were no real downloaded GEO expression files. Original manifest contents and confirmation receipts were archived unchanged. They were not automatically rehashed to approve new paths. Fresh synthetic results occupy the new active paths. Migration audit files remain local at `.workflow/migration.json`. No unattributed file was deleted.

Old root folders and their tracked `.gitkeep` files were removed. `data/.gitkeep` replaces them. Git excludes all real/synthetic runtime outputs, RDS/H5/H5AD/MTX/FASTQ and credential patterns; only Skill code, documentation, test generators and the empty data placeholder are intended for the commit.

## Remaining limits

Public metadata still needs evidence review and explicit user grouping. Download instructions alone do not authorize choosing groups. No real large-matrix end-to-end run was performed. Native zellkonverter remains unverified on this environment. Arbitrary author formats, ambiguous pooled ownership and matrices larger than available memory may need additional work. SHA256 provenance detects local changes; it is not an author signature. Existing RDS and changed confirmations require explicit archiving/reconfirmation; there is no automatic overwrite or legacy-receipt conversion.
