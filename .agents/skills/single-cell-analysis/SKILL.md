---
name: single-cell-analysis
description: Run formal PCA, regression diagnostics, design-aware Harmony, clustering review and saved-result audit on an existing QC Seurat RDS. Use for QC后分析, PCA/Harmony/聚类, or resuming analysis; final resolution requires the user's choice.
---

# Formal single-cell analysis

Work from the repository root. Use fixed scripts through `analysis/run_workflow.py GSE --root ROOT`. Input defaults to completed canonical QC-3/QC-2; `--input PATH` supports an existing QC object. Never restart downloads/QC from an analysis-only request. All R executes through `runtime/r45/run_r45.py` in existing R 4.5.0. Do not change packages or the environment.

Configuration is `analysis/default_config.json`, optionally overridden with `--config JSON`. Output is `results/<GSE>/03_analysis/`. Validate checkpoint/config/input provenance before reuse; never overwrite valid RDS files.

1. `--step pca` freshly runs NormalizeData, FindVariableFeatures, ScaleData with baseline `vars.to.regress="pMT"`, and RunPCA. Never reuse DoubletFinder PCA. RNA is default; DecontX-corrected expression requires an explicit assay choice recorded with input/parameters. Read PC selection, elbow, regressions and warnings in `01_pca_regression/`.
2. Extra regression is never default. When baseline diagnostics show concrete RNA-total or cycle risk, configure only the relevant `sensitivity` comparison (`rna_total` or `cell_cycle`). Results stay in `regression_sensitivity/` without replacing baseline. Sample separation first prompts design/integration review.
3. `--step harmony` examines sample, group, batch, library_id and sequencing_run. Never correct by `group`. Reliable technical fields require explicit designation in config; sample correction requires reviewed design and explicit selection. Fully confounded biology cannot be claimed separable by Harmony. Keep unintegrated UMAP and report limitations. Read batch design, mixing metrics, marker comparisons and warnings in `02_harmony/`; mixing alone does not validate integration.
4. `--step cluster-review` scans resolutions on the selected representation and one UMAP. Review `03_clustering/cluster_resolution_review.pdf`, `clustree.pdf`, cell distributions, stability and regression/cluster associations. Ask the user to choose final resolution from the evidence. Never select solely from UMAP or a fixed default.
5. `--step finalize --resolution-decision JSON` requires `confirmed_by="user"`, literal `user_statement`, and a resolution in the actual review. Only then save `seurat_analysis.rds` and formal labels. `--step all` stops at review without a selection. Simulated selections belong only to explicitly marked regression fixtures.

Run `workflow_audit/run_audit.py GSE --root ROOT` to summarize existing evidence into `03_analysis/99_workflow_audit/`. It does not rerun expression analysis. Missing evidence/failures remain explicit; updated reports are versioned. This scientific report differs from the private ten-file download audit.

Fix technical failures with evidence/tests. Uncertain parameter effects need limited comparisons. Biological tradeoffs require the user's decision: never delete samples/cell populations or erase biological structure on the user's behalf. Report actual RDS/CSV/plot validation and limitations; a skipped method is not a successful computation.
