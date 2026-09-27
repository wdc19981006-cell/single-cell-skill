---
name: single-cell-qc
description: Independently run cell QC, per-sample DoubletFinder and cell-cycle scoring on an existing seurat_raw.rds. Use for “对 GSE156625 进行质控”, “质控 GSE156625”, “继续 QC”, or “继续处理 GSE156625” when a validated raw object already exists. Requires user decisions for flagged samples or thresholds.
---

# Single-cell QC

Work from the `single-cell-skill` repository root. This workflow starts independently whenever `data/<GSE>/seurat_raw.rds` already exists, including days after the Loader finished. Check that file first. If missing, stop with “未找到 seurat_raw.rds，请先完成数据下载和原始对象构建。” Do not restart GEO Stage A, redownload raw expression, or rebuild the raw RDS unless the user explicitly requests that separately. The fixed R scripts in `qc/` perform all calculations; the Skill only invokes them, reads reports and handles user decisions.

Use `bash qc/r450_rscript.sh` for dependency checks, precheck, tests, and QC. It invokes only `D:/R/R-4.5.0/bin/Rscript.exe` and applies a project-local exit workaround for this host's `cli` cleanup crash; it does not change the analysis or the installed package library. First verify that executable exists, then report `.libPaths()`, `R.version.string`, `packageVersion("Seurat")`, and `packageVersion("DoubletFinder")` from that R. Its intended library is `D:/R/R-4.5.0/library`. If that executable is missing, stop; never switch to another R found on `PATH`. Never install missing packages automatically; report each missing package, R version, and `.libPaths()`.

1. Run `bash qc/r450_rscript.sh qc/qc_precheck.R <repository-root> <GSE>` and read `data/<GSE>/qc/qc_report.csv`. Precheck may create the QC directory and report, but never filters the input object.
2. If `RUN/status` is `NEEDS_USER_DECISION`, stop. Report each flagged sample with group and cell count, plus any threshold or species conflict. Ask the user to choose the action. Never remove a sample, change a threshold, group, expected doublet rate, or cell-cycle regression on your own.
3. After an explicit user answer, write `data/<GSE>/qc/qc_decision.json` with `decision_type`, `affected_samples`, `confirmed_by: "user"`, their literal `user_statement`, and `removed_samples`, `kept_low_cell_samples`, or `parameter_overrides` as needed. Use decision types `remove_low_cell_samples`, `keep_without_doubletfinder`, or `override_parameters` for the three low-cell choices. `accept_default_conflicts: true` records an explicit choice to retain conflicting defaults. R validates decisions against the raw object. Rerun precheck and read the new report before continuing; never re-ask for the same persisted choice.
4. When precheck is `PASS`, run `bash qc/r450_rscript.sh qc/run_qc.R <repository-root> <GSE>`. Read the single `qc_report.csv`, then verify `qc/seurat_qc.rds` and all six requested PDF/PNG pairs are present and nonempty. Normal QC status is `COMPLETE_QC`; approved skipped-doublet samples retain `QC_COMPLETE_WITH_UNEVALUATED_DOUBLETS`. On a technical error, report the failure; do not silently change analysis code and rerun.

The fixed R pipeline filters cells only, runs DoubletFinder separately for each sample, records both `DF` and homotypic-adjusted `DF_adj`, saves doublet plots before removal, then scores cell cycle without regression. The final object keeps raw RNA counts and QC metadata, with no temporary PCA, UMAP, or cluster state. The source `seurat_raw.rds` is never overwritten. No integration, clustering interpretation, annotation, or differential expression is part of this Skill.

If a user explicitly keeps a low-cell sample without DoubletFinder, record `DF = "NotEvaluated"`, `DF_adj = "NotEvaluated"`, and `doublet_status = "SKIPPED_LOW_CELL"`; never mislabel it as a predicted singlet. The final run status must be `QC_COMPLETE_WITH_UNEVALUATED_DOUBLETS`, and the report must list the unevaluated samples and cell counts.

The historical `qc/references/01_Seurat_1.R` is reference material only. Never execute it as the production pipeline. A synthetic/mock test verifies control flow; only a test that invokes the installed DoubletFinder functions verifies real doublet computation.

All formal QC figures are generated and saved through `qc/plot_utils.R`; use its `save_sc_plot()` and shared UMAP theme/coordinate helpers for later single-cell plots. Do not scatter separate `ggsave()` size decisions through analysis scripts. The Loader's `STATUS: COMPLETE_RAW` means only the raw object exists and never grants permission to start QC.
