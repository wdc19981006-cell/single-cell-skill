# Formal analysis in the existing R 4.5.0 runtime

`analysis/run_workflow.py GSE123 --root PROJECT_ROOT` starts from the completed
QC-3 RDS, then QC-2, with a read-only fallback to the historical
`data/GSE123/qc/seurat_qc.rds`. `--input PATH` accepts an existing RDS. An explicit
input on `--step pca` must contain an RNA counts layer, sample and finite pMT.
An explicit input on subsequent steps must be the appropriate preceding formal
analysis object. No download or QC computation is launched.

All outputs are under `results/GSE123/03_analysis/`, with stage-specific plots,
tables, configuration, result records and logs. The top-level RDS files are
`seurat_pca.rds`, `seurat_unintegrated.rds`, `seurat_harmony.rds` when integration
actually ran, `seurat_cluster_review.rds`, and `seurat_analysis.rds` only after
a user decision. The audit is refreshed after review/finalization using existing
stage tables. Failed computation keeps its staging directory and runtime log.

Independent commands use `--step pca`, `--step harmony`, `--step cluster-review`
or `--step finalize`. `all` runs through cluster review and then stops for the
user's choice, unless an explicit `--resolution-decision` file was provided.
A later `all` validates an already authorized final result.
Successful stages are checksum-bound to their input and relevant settings;
resuming validates existing files and reloads the RDS instead of rebuilding.
Publication is atomic and exclusive. A changed input, configuration or output
does not authorize replacement. Use a separate project root for alternatives.

Pass a JSON file with `--config`; only keys in `default_config.json` are accepted.
The baseline always reruns NormalizeData, variable-feature selection, ScaleData
with pMT alone, and PCA. It reuses `qc_pc_selection` and its standard-deviation
share rule; squared variance is reported separately. An alternative counts
assay such as `decontX` must be explicitly selected with `{"assay":"decontX"}`.
Original RNA counts stay in the object. Never substitute normalized values for
a missing counts layer.

`regression_diagnostics.csv` records Spearman relationships of selected PCs with
pMT, nCount_RNA, nFeature_RNA, S.Score and G2M.Score, including unavailable or
constant variables. Optional `"sensitivity":["rna_total"]` compares pMT against
pMT+nCount_RNA; `"cell_cycle"` compares against pMT+S.Score+G2M.Score. A requested
comparison only runs if its target has measured absolute correlation at or above
`risk_correlation` (default 0.5). This is a diagnostic trigger, not proof of a
technical cause. Each comparison has its own RDS, settings, plots and subspace
overlap table. It never replaces the baseline or silently changes clustering.

Harmony requires a confirmed technical field, for example:

```json
{"harmony":{"technical_fields":["batch"]}}
```

The presence of a field named batch/library_id/sequencing_run alone does not
prove technical meaning. Without confirmed information the design table and
unintegrated UMAP are saved, and Harmony is explicitly skipped. `group` is
rejected as a correction variable. Sample integration additionally requires
`allow_sample:true`. Rank-confounding with biological group, or unavailable
group metadata, blocks correction unless a literal user acknowledgement is
provided in `harmony.confounding_acknowledgement` using `confirmed_by:"user"`
and `user_statement`. The confounding warning remains in both cases; the result
cannot identify technical versus group effects. No correction-strength parameter
is automatically enlarged. Pre/post kNN label fractions and configured-marker
neighborhood associations supplement UMAPs; absent markers are reported.

Clustering scans all configured resolutions on one graph and one saved UMAP,
records cluster/sample/group sizes, numeric metadata associations, and adjacent
resolution adjusted Rand indices. A clustree and a multi-page resolution review
PDF support biological review. Neither stability nor visual mixing chooses a
resolution. Finalize only after the actual user's statement is recorded:

```json
{"confirmed_by":"user","user_statement":"<literal user choice>","resolution":0.6}
```

Run `--step finalize --resolution-decision DECISION.json`. The selected value
must occur in the reviewed scan. `selected_resolution.yaml` and formal labels
are then created together. An unselected review object has no seurat_clusters.

Validation uses `tests/test_analysis_workflow.py` (orchestration and preservation),
`tests/test_analysis.R` (scientific guards and actual PCA), and
`tests/run_analysis_smoke.py --output NEW_DIRECTORY` (actual synthetic PCA,
Harmony, clustree, finalization, audit and resume). Synthetic choices are test
fixtures and never decisions for a real GEO dataset. Every R process runs through
`runtime/r45/run_r45.py`; no package or R installation is changed.

API references consulted for the installed-route compatibility checks:
[Harmony Seurat method](https://github.com/immunogenomics/harmony/blob/master/R/RunHarmony.R),
[clustree input columns](https://lazappi.github.io/clustree/articles/clustree.html),
[ggraph edge guide](https://ggraph.data-imaginist.com/reference/guide_edge_colourbar.html),
[ggplot2 explicit guides](https://ggplot2.tidyverse.org/reference/guides.html).
The installed clustree/ggraph string guide was rejected by ggplot2 4.0.2; explicit
edge-guide objects repair the plot without changing clustering or packages. The
regression test builds an actual clustree grob to cover that observed failure.
