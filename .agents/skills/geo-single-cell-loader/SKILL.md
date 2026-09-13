---
name: geo-single-cell-loader
description: Search and inspect GEO single-cell datasets, summarize sample and disease metadata for user-confirmed grouping, download processed expression matrices, detect common 10x/H5/H5AD/text formats, and create validated Seurat objects with standardized metadata. Use when the user asks to inspect, download, load, or prepare GEO scRNA-seq or snRNA-seq datasets for Seurat.
---

# GEO single-cell loader

Use the directory containing `.agents/skills/geo-single-cell-loader` as repository root. Run all commands from that root. Store data only under `datasets/<GSE>/`, manifests under `manifests/`, results under `output/<GSE>/`, and logs under `logs/`. Never store downloaded data inside the Skill.

## Stage A — inspect and report

1. **Search first.** Prefer available bio-server MCP `search_geo` and `get_geo_info`. If unavailable, use `scripts/inspect_geo.py GSE...`; keyword search uses `--search`. Do not assume an unavailable MCP server exists.
2. **Do not bulk-download before presenting sample information.** Read Series, every relevant GSM, associated publications and publicly linked supplementary sample/clinical metadata. The script fetches metadata and abstracts, not all publication supplements: follow those links yourself. Its format labels are candidates until contents are verified. A development `--max-samples` report is incomplete and cannot represent all samples.
3. Enrich the generated sample report using public evidence. Show `sample`, `author_sample` when available, `tissue`, `disease`, `source_type`, `sample_description`; optional fields only when reliably mapped. Record sequencing modality (scRNA/snRNA/etc.), study context, raw/filtered availability, separate clinical/cell metadata, and unresolved facts. Do not silently invent values to fill blank fields.
4. **group must always be confirmed by the user.** Say: **group 尚未创建，请确认分组方式。** Show facts before asking. A request to “download and read GSE229413” still starts here. Never infer group from disease or filenames. PDAC uninvolved pancreas still has disease PDAC, not Healthy.

## Stage B — explicit grouping, download, create

Read [metadata-schema.md](references/metadata-schema.md) before building the manifest and [format-routing.md](references/format-routing.md) for the selected format.

1. After the user's reply, expand their rules into an exact `sample -> group` JSON map. Preserve their literal statement, `confirmed_by: user`, and all sample keys. This records user intent; it is not an authentication system. Never synthesize a user statement or fill unconfirmed groups. Synthetic fixtures use explicitly labeled simulated approval and never clinical conclusions.
2. Resolve all required biology and file mapping facts. Adopt `input_plans` only after evidence review. Canonical trio names can be proposed from GSM-owned links; patient/specimen identities cannot. Prefer filtered over raw **within one evidenced sample**. For pooled matrices use an exhaustive explicit cell-to-sample map and record its provenance. Stop on unresolved mapping.
3. Build `manifests/<GSE>_sample_manifest.csv` with `build_manifest.py --report ... --confirmation ... --output ...`. **Use manifest as the source of truth.** Exact confirmation receipts bind to manifest contents; edits require reconfirmation. Run `--validate` before downloading. Do not bypass these checks by manually writing receipt hashes.
4. Run `download_processed.py manifests/<GSE>_sample_manifest.csv`. Small public sample metadata may be fetched during Stage A; bulk count matrices may not. For archives, specify exact members and canonical destinations in manifest; no broad extraction. Log download date and SHA256. Inspect unfamiliar input contents and counts provenance; a filename or integer-looking matrix alone does not prove raw counts. If mapping changes, return to confirmation.
5. Check dependencies using `Rscript scripts/check_dependencies.R` (use the actual script path). Use a verified R installation with Seurat; optional missing format packages should stop only the affected route. Prefer `zellkonverter` native R for H5AD. If it is unavailable and a prevalidated Python with `anndata/numpy/scipy` exists, set `GEO_SINGLE_CELL_PYTHON` explicitly for the fallback. Never silently create or install a Python/Conda environment during conversion.
6. Run `Rscript .agents/skills/geo-single-cell-loader/scripts/build_seurat.R . manifests/<GSE>_sample_manifest.csv`, then the independent `validate_seurat.R . manifests/<GSE>_sample_manifest.csv output/<GSE>/seurat_raw.rds`.

**database/sample/tissue/disease/source_type/group are required in final Seurat objects. patient/specimen/cohort/treatment are optional and must never be guessed.** GSM is preferred; use documented author sample IDs only when a GSM is not one matrix sample. Source names and descriptions remain available for traceability.

**Never silently continue after sample mapping ambiguity.** Reject missing, duplicate, unmatched or anomalous mappings. Prefix barcodes with sample before combining `scelist <- list(...)`. Assign metadata by stable keys, verify cell alignment and retained sample sets. Empty samples stop the run. No final RDS is written until in-memory and serialized validation succeed.

**Stop after creation of seurat_raw.rds; do not run downstream scRNA-seq analysis.** The explicitly requested `CreateSeuratObject(min.cells=3, min.features=200)` construction thresholds are recorded; perform no additional QC, normalization, scaling, PCA/UMAP, clustering, doublet detection, annotation, differential expression or enrichment. Never call NormalizeData, SCTransform, FindVariableFeatures, ScaleData, RunPCA, RunUMAP, FindNeighbors, FindClusters or DoubletFinder.

Report outputs, sample/group cell counts, count source, warnings and limitations. FASTQ/SRA-only datasets: “该数据没有可直接读取的processed expression matrix。需要FASTQ/SRA重新定量，超出当前V1范围。” Do not run Cell Ranger.

See [test-datasets.md](references/test-datasets.md) for real discovery cases and synthetic regression tests; development authorization does not assign real biological groups.
