# Fixed formal analysis and targeted review

Run from the code repository with `/c/Python312/python.exe analysis/run_workflow.py DATASET`.
The shared default in `runtime/output_config.json` writes to
`D:/Desktop/codex-result/single-cell-skill/DATASET/03_analysis/`.
`--output-root PATH` changes the collection root (DATASET is appended once).
`--root OLD_WORKSPACE` retains the explicit legacy `OLD_WORKSPACE/results/DATASET`
layout. Code and R package paths never derive from the result root. Non-GEO
datasets use safe ASCII names, such as `Study_A`; path separators, traversal,
Windows device names and trailing dots are rejected.

Input defaults to saved QC-3. `--input ABSOLUTE_RDS` reuses historical QC without
moving it or restarting downloads/QC. QC-2 fallback is useful only with an
explicit `assay=RNA` configuration: default `decontX` missing or invalid counts
stop the run. RNA counts, fractional corrected counts, original QC metadata,
and aligned gene/cell IDs are checked. Corrected counts are never rounded.

Baseline: LogNormalize → 2,000 HVGs → ScaleData(pMT) → PCA → Harmony(sample)
→ one neighbor graph and UMAP → multiple resolutions → regression review.
The reference standard-deviation PC selection rule is retained. `all` stops for
the user's final resolution; `--step pca|harmony|cluster-review|regression-review|finalize`
supports independent execution. Each stage has a short `*_result.json` containing
status, inputs/outputs, parameters, statistics, unresolved issues and next action.

Harmony defaults to `enabled=true,batch_field=sample,allow_sample=true`.
Sample/group confounding is reported and unintegrated results are always kept;
the approved sample policy does not invent a per-run `confirmed_by=user` decision.
Missing/blank/one-level batch fields skip with an explicit reason.
An override such as `{"harmony":{"batch_field":"batch","technical_fields":["batch"]}}`
requires valid technical metadata. `group` is forbidden. Confounded technical
overrides retain the existing explicit acknowledgement guard. Mixing alone does
not establish preservation of biology.

After baseline clustering, `regression_review.R` diagnoses all PCs (only used
PCs trigger comparisons), cluster eta squared and per-cluster metric profiles.
`risk_correlation=0.5` and `regression_review.cluster_eta_warning=0.5` are diagnostic
triggers, not scientific criteria for adopting regression. With `auto_numeric=true`,
only risky `feature` (pMT+nFeature_RNA) and/or `count` (pMT+nCount_RNA) run.
`sensitivity` can restrict requested comparisons when `auto_numeric=false`;
the historical name `rna_total` maps to `count`. A no-risk target does no extra PCA.
Cycle risk first produces proliferation marker evidence. A `cell_cycle` comparison
also needs `regression_review.cell_cycle_review` with `confirmed_by=user`, their
literal `user_statement`, and `marker_interpretation`. Scores alone never authorize it.

Each actual comparison starts from the exact hashed QC input, same cells, assay,
HVGs, seed, PCA bound, Harmony design, used PC count and comparison resolution.
`comparison_resolution` must occur in the baseline scan (default 0.3). Outputs in
`01_pca_regression/regression_sensitivity/SCHEME/` include config/log/RDS, PCA
subspace overlap, ARI, cluster counts/transitions, marker expression and neighborhood
structure, technical associations, phase and sample/group composition. Expression
does not change with regression: marker comparisons concern cluster/neighborhood
structure. UMAP distance is not a selection metric. Recommendations are
`BASELINE_PREFERRED`, `SIMILAR_RESULTS`, or `USER_REVIEW_REQUIRED`, explicitly
rule-based advice rather than proof of biological correctness.

Results and checkpoints are immutable and input/config/checksum-bound. Use
`--run-id decontx-v1` for an alternative under `03_analysis/runs/decontx-v1/`.
No earlier run is searched for Harmony, clustering or checkpoints. Explicit
`--step regression-review --input REVIEW_RDS --run-id REVIEW_ID` supports a
targeted isolated review of an existing baseline. Completed computations survive
validation/publication interruptions; a native crash remains a failure. State JSON
replacement tolerates brief Windows locks and preserves pending files on failure.

To adopt a comparison, use a new `--run-id` and `--regression-decision FILE`:

```json
{"scheme":"feature","confirmed_by":"user","user_statement":"<literal choice>",
 "comparison_config":"<absolute completed feature/comparison_config.json>"}
```

The QC input hash, assay and selected variables must match. PC selection and
resolution scanning run anew; a comparison never automatically becomes formal.

Finalize with `--resolution-decision FILE`, containing `resolution`,
`confirmed_by=user` and the literal `user_statement`. If regression review requires a biological choice, keeping this
run also requires `regression_scheme="baseline"` and a literal
`regression_statement` in that decision. Candidate adoption uses a separate run.
The resolution must be in the reviewed scan. `03_clustering/selected_resolution_0.3/` (for 0.3) then contains
marker DotPlot/FeaturePlots, QC/cycle distributions, cluster sizes and sample/group
composition. Marker plots support annotation review; no cell types are assigned.

`workflow_audit/run_audit.py DATASET [--run-id ID]` reads saved evidence only.
It distinguishes PASS, WARNING, ACTION_REQUIRED, RESOLVED, FAILED and
MISSING_EVIDENCE, groups repetitive warnings, verifies final manifest checksums,
and resolves historical resolution gates only when validated final evidence agrees.
Missing evidence is never inferred to pass. Changed evidence versions the audit;
the private ten-file GEO download audit remains separate.

Routine use: run fixed commands and read JSON summaries; inspect full logs only
for actual errors. Targeted review: run only the matching diagnostic/comparison.
Development: modify modules and run affected regression/acceptance tests, without
recomputing a whole study. Every R process uses `runtime/r45/run_r45.py`; no
packages or system R settings are changed.

API references: [Seurat ScaleData](https://satijalab.org/seurat/reference/scaledata),
[Seurat cycle regression](https://satijalab.org/seurat/articles/cell_cycle_vignette),
[celda DecontX](https://bioconductor.org/packages/release/bioc/vignettes/celda/inst/doc/decontX.html).
The installed `celda::decontX` route is preserved; no package migration is performed.
