# R45 runtime and workflow validation — 2026-09-28

This report records executed checks for the runtime refactor. `GSE225857`
was read only. No real GEO expression data was downloaded or rebuilt; no R
package was installed, upgraded, or removed.

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
| Real GSE QC / real-cohort DoubletFinder | **NOT RUN**; no biological user decision to initiate QC was part of this refactor. |

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
