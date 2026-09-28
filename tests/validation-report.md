# R45 runtime and workflow validation — 2026-09-28

The first section records the earlier runtime refactor checks. `GSE225857`
was read only during that stage; no real GEO expression data was downloaded or
rebuilt then. The later GSE231993 acceptance is recorded separately below. No
R package was installed, upgraded, or removed.

## Environment observed

| Item | Actual result |
|---|---|
| R executable | `D:/R/R-4.5.0/bin/Rscript.exe` |
| R version | `4.5.0` (`R version 4.5.0 (2025-04-11 ucrt)`) |
| `R.home()` | `D:/R/R-4.5.0` |
| `.libPaths()` | `D:/R/R-4.5.0/library` |
| Seurat / SeuratObject | `5.4.0` / `5.4.0` |
| DoubletFinder | `2.0.6` |
| Other core packages | Matrix `1.7.3`, jsonlite `2.0.0`, yaml `2.3.12`, ggplot2 `4.0.2`, patchwork `1.3.2` |
| Python used | `C:/Python312/python.exe` (3.12.0) |
| CLI exit workaround | Active in child process; `library(cli)` exited 0 and intentional `stop()` exited 1 |

The QC-route healthcheck exited 0 and wrote `.runtime/r45/health.json`. No
`C.UTF-8` startup warnings were observed in the runtime acceptance checks.
The installed `cli` emitted a package-built-under-R-4.5.3 warning in the
compatibility shell test; the process still used R 4.5.0 and exited 0.

## Executed results

| Check | Observed result |
|---|---|
| Python `unittest discover -s tests -p 'test_*.py'` | **PASS, 81 tests**, no skips. Covers Loader routes, checkpoint/artifact recovery, audit resume, launcher native-crash classification, interruption versus validation failure, and QC orchestration/decision gate. |
| `tests/runtime/run_tests.py --real-fixture` | **PASS**: `1+1`, `library(cli)`, ggplot2, SeuratObject, Seurat and DoubletFinder each exited 0. Intentional `stop("RUNTIME_TEST_ERROR")` exited 1. RDS and Chinese UTF-8 TXT/JSON round trips passed. |
| `tests/test_r450_launcher.sh` | **PASS**: deprecated compatibility command invokes the unified runtime, keeps normal success and R error exit 1, and rejects `--vanilla`. |
| `tests/run_regressions.py` synthetic Loader | **PASS**: four-input raw build and independent validation (16 cells, 250 genes), pooled matrix/cell map, normalized H5AD failure/retry, raw retention, and refusal to overwrite an existing final RDS. Used an isolated ignored synthetic workspace. |
| Existing raw artifact recovery | **PASS**: synthetic `sample_info.txt` was set to `BUILD_FAILED`; independent validation of the existing RDS restored `COMPLETE_RAW`. The RDS SHA256 was unchanged and no build/download ran. |
| Other Loader R tests | **7/7 PASS**: `test_10x_literal_na.R`, `test_final_object.R`, `test_global_min_cells.R`, `test_pooled_inputs.R`, `test_seurat.R`, `test_stream_text.R`, `test_utf8_manifest.R`. Python H5AD fallback was explicitly set to `C:/Python312/python.exe`. |
| QC synthetic `test_qc.R` | **PASS**: precheck, user-decision gates, strict thresholds, mock DoubletFinder flow, cell cycle, plots, final object and raw hash preservation. |
| Real DoubletFinder `test_qc_real_df.R` | **PASS**: installed DoubletFinder performed `paramSweep`, pK selection and both classification calls on a synthetic 130-cell object; exit 0. This was not a real GEO cohort. |
| Existing QC artifact recovery | **PASS**: independent read-only QC validation through `qc/run_workflow.py`; no new QC build. |
| `GSE225857` read-only raw validation | **PASS, exit 0**: 196,473 cells and 17,066 genes. SHA256 before and after: `6691e028ea8ecfda86711956c3b1b5c1d91b184fc87550d94822e3922bd3a7ab`. No download or build was started. |
| Native crash classification | **PASS (simulated exit codes)**: 139, −1073741819 and 3221225477 remain `R_NATIVE_CRASH` even after stdout contains `PASS`. No real native crash was induced. |
| Real GSE QC / real-cohort DoubletFinder | **NOT RUN at this earlier checkpoint**; see the real acceptance section below. |

For this previously cleaned dataset, `validate_seurat.R --existing` checked
the raw object structure, retained user-confirmed groups and dimensions from
`sample_info.txt`. Its original manifest was no longer local, so this check
could not repeat the full manifest/receipt comparison performed at build time.

The first regression attempt stopped at `test_stream_text.R` because the new
test runner omitted that test's required Python executable argument. The runner
was corrected; the complete rerun passed. An earlier runtime attempt stopped
at a syntax error in the new archived-artifact validation branch. That syntax
was corrected, and `GSE225857` then passed read-only validation. Failed
attempts are not counted as passes.

`git diff --check` and the tracked-code R invocation search passed before
commit. Commit, push and final worktree status are recorded in the delivery
summary after they occur.

## REAL GSE231993 FULL ACCEPTANCE — 2026-09-28

The final tested workflow code is commit
`a77ec822fe02069edc8c9a521080d8b30b796480`. A post-build change in that
commit only adds `route=...` to missing-package error messages; its success
path was unchanged and unit-tested. The run used an isolated
`.runtime/acceptance-workspace` checkout whose `data/GSE231993/` directory did
not exist at the start. Existing project data and Release artifacts were not
used. The previous user-confirmed groups were preserved exactly:
GSM7307094–GSM7307101 = UC (8 samples) and GSM7307102–GSM7307105 = HC
(4 samples).

| Check | Observed result |
|---|---|
| R runtime | R `4.5.0`, `R.home()` `D:/R/R-4.5.0`, `.libPaths()` `D:/R/R-4.5.0/library`; Seurat `5.4.0`, DoubletFinder `2.0.6`. All 16 acceptance R records used the fixed R45 launcher. Native crashes: **0**. |
| Code tests | Python unittest **86/86 PASS**; runtime acceptance and compatibility launcher **PASS**; Loader and QC synthetic regression, including installed DoubletFinder and direct production QC entrypoints, **PASS**. Actual route checks passed for `base,text,10x_h5` once each, and H5AD Python fallback passed with `anndata/numpy/scipy` while native `zellkonverter` was unavailable. |
| GEO download | **36/36 DOWNLOADED**, 238,349,484 bytes from the 12 GSM-owned 10x trios. The duplicate series RAW.tar was excluded with a recorded reason. No source file was reused. |
| Loader | 12 samples, **60,665 cells × 25,953 genes**; reader `Seurat::Read10X` for 12 unique inputs; build `135.88 s`; independent raw validation exited 0. Historical expected 60,665 × 25,953: **MATCH**. UC: 37,967 cells; HC: 22,698. |
| Raw object and handoff | One sparse integer RNA counts layer, required metadata and user groups, no normalized or downstream state: **PASS**. `STATUS: COMPLETE_RAW` and the QC handoff question were produced; Loader did not start QC. Read-only raw revalidation after QC exited 0. |
| QC precheck | Initially `NEEDS_USER_DECISION`: predicted 13,814 / 60,665 cells removed (**22.77%**), with 11 samples individually over 15%. User replied **“继续使用当前阈值”**. Decision type `accept_qc_removal`; affected samples GSM7307094–GSM7307105. The decision listed real sample IDs because `ALL` is a report summary row. Rerun precheck: **PASS**. |
| Precheck biology | Human; pMT/pRP/pHB matched 13/101/11 features; cell-cycle S/G2M matched 42/52 features. Thresholds remained nFeature >200 and <4000, nCount >500 and <30000, pMT <25, pHB <1, pRP <100. |
| Cell QC | 60,665 → **46,851** cells; **13,814 removed (22.77%)**. |
| Per-sample DoubletFinder | `DF` predicted **1,872** doublets; homotypic-adjusted `DF_adj` predicted and removed **1,776**. After doublet removal: **45,075** cells. No sample skipped DoubletFinder. |
| Final QC object | **45,075 cells × 25,953 genes**. Cell cycle: G1 **26,548**, S **10,082**, G2M **8,445**. `qc/validate_qc.R` exited 0. An additional read-only comparison found exact raw counts for every retained cell and identical original `database`, `sample`, `group`, `tissue`, `disease`, and `source_type` metadata. Final object has required QC/doublet/cell-cycle fields and no PCA, UMAP, graphs, neighbors, normalized layer, SCT or integrated assay. |
| Plots | All six PDF/PNG pairs (**12 files**) exist and have nonzero size: `before_QC`, `after_QC`, `doublet_umap`, `doublet_vlnplot`, `cell_cycle_phase`, `cell_cycle_score`. All 12 plot report statuses are PASS; all six PNG metadata values round to **300 DPI**. |
| Raw immutability | SHA256 before and after QC: `8c110ea1813777e715b38141385069aeb0d2b36063c60f9a2255d5fd463f6e03` — **MATCH**. |
| Audit | [Successful immutable audit](https://github.com/wdc19981006-cell/single-cell-skill-audit/tree/main/GSE231993/20260928T120755Z-75bb7029). Contains exactly the ten policy files, including file selection, download profile, input routing, build profile, validation and execution log. |

The first audit attempt stopped before expression download because the sandbox
could not reach the private audit repository. That failure was uploaded as a
[separate immutable audit](https://github.com/wdc19981006-cell/single-cell-skill-audit/tree/main/GSE231993/20260928T120415Z-ef48cfc4).
The resumed acceptance performed all 36 GEO downloads. The first real QC
precheck exposed a top-level `sys.frame(1)` entrypoint error; the direct
entrypoints were fixed and regression-tested. A second precheck stopped as
designed for the >15% biological decision. Neither interruption changed raw
counts or thresholds. The completed QC run took `1211.69 s` and returned
`COMPLETE_QC`.
