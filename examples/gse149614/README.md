# GSE149614 pooled scRNA-seq workflow

This runbook records the dataset-specific logic for the repository's `geo-single-cell-loader` Skill. The public GEO Series has 21 GSM records from 10 people with hepatocellular carcinoma (HCC). The four specimen sites are primary liver tumor (10), adjacent non-tumor liver (8), portal vein tumor thrombus (2), and metastatic lymph node (1). Adjacent non-tumor tissue still has HCC as the donor disease background.

## Public inputs and evidence

| Input | Purpose |
|---|---|
| `GSE149614_HCC.scRNAseq.S71915.count.txt.gz` | Cell Ranger raw counts, one pooled genes × cells TSV; 165,349,783 compressed bytes at inspection time. |
| `GSE149614_HCC.metadata.updated.txt.gz` | Public per-cell `Cell`, `sample`, `site`, and `patient` fields. |
| `GSE149614_HCC.scRNAseq.S71915.normalized.txt.gz` | Normalized, log-transformed expression; excluded from raw Seurat construction. |

The metadata has 71,915 unique cell IDs and exactly 21 author sample names. A 512 KiB HTTP range probe of the count file showed 71,915 cell IDs in its header, with an exact set match to the metadata; its data rows begin with gene IDs. The first header field is already a cell ID, so the text reader uses `feature_column=__row_names__`, `orientation=genes_by_cells`, and `delimiter=tab`.

Sources: [GEO GSE149614](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614), [primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC9357016/), and its Supplementary Data 1. Original sequencing reads require controlled EGA access; this workflow downloads the public processed raw-count matrix from GEO.

## Run sequence

1. In Codex, use the global `geo` MCP `get_geo_info("GSE149614")` followed by `list_geo_files("GSE149614")`. Require a complete 21-GSM inventory and file ownership. Save the standard `data/GSE149614/.workflow/inspection.json`, `sample_report.csv`, and `sample_info.txt` as described by the Skill.
2. Inspect the small public metadata file at `data/GSE149614/.workflow/public_metadata/GSE149614_HCC.metadata.updated.txt.gz`. `prepare.py check` verifies its `Cell → author sample → GSM` map against the GEO sample records and report. It does not assign groups or download the expression matrix.
3. Show all 21 samples and ask the user to assign `group`. Save their literal reply and an exact GSM-to-group map in `data/GSE149614/.workflow/group_confirmation.json` with `confirmed_by: "user"`. Do not derive group from `site`, `patient`, diagnosis, or sample suffix. Every GSM needs one nonblank group.
4. Run `prepare.py prepare`. It writes one exhaustive `cell_map.csv` and adds the shared count-file download plan to every report row. Then build and validate the canonical manifest.
5. Run `download_processed.py` on the validated manifest. It downloads the count file once under `raw/` and records SHA256, size, URL, and UTC time in `.workflow/download.json`. Run `prepare.py verify-download` to compare all matrix-header cell IDs with the public metadata and check its full row layout and gzip integrity.
6. Check R dependencies, then build and independently validate `seurat_raw.rds`. The Skill creates only a raw RNA counts layer with `min.cells=3` and `min.features=200`; no normalization or downstream analysis is part of this run.

From the repository root on this Windows Git Bash setup, with the repository path adjusted if needed:

```bash
ROOT=D:/CodexProjects/single-cell-skill
PY=/c/Python312/python.exe
RS=/d/R/R-4.5.3/bin/Rscript.exe

mkdir -p data/GSE149614/.workflow/public_metadata
curl --fail --location \
  --output data/GSE149614/.workflow/public_metadata/GSE149614_HCC.metadata.updated.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.metadata.updated.txt.gz
"$PY" "$ROOT/examples/gse149614/prepare.py" check
# Create group_confirmation.json only after the user confirms every group.
"$PY" "$ROOT/examples/gse149614/prepare.py" prepare
"$PY" "$ROOT/.agents/skills/geo-single-cell-loader/scripts/build_manifest.py" \
  --report data/GSE149614/.workflow/sample_report.csv \
  --confirmation data/GSE149614/.workflow/group_confirmation.json \
  --output data/GSE149614/.workflow/sample_manifest.csv
"$PY" "$ROOT/.agents/skills/geo-single-cell-loader/scripts/build_manifest.py" \
  --validate data/GSE149614/.workflow/sample_manifest.csv
"$PY" "$ROOT/.agents/skills/geo-single-cell-loader/scripts/download_processed.py" \
  data/GSE149614/.workflow/sample_manifest.csv
"$PY" "$ROOT/examples/gse149614/prepare.py" verify-download
"$RS" .agents/skills/geo-single-cell-loader/scripts/check_dependencies.R
"$RS" .agents/skills/geo-single-cell-loader/scripts/build_seurat.R \
  . data/GSE149614/.workflow/sample_manifest.csv
"$RS" .agents/skills/geo-single-cell-loader/scripts/validate_seurat.R \
  . data/GSE149614/.workflow/sample_manifest.csv data/GSE149614/seurat_raw.rds
```

The pooled text matrix has roughly 1.85 billion numeric fields. The current text reader materializes a dense table before sparse conversion, so Seurat construction may need more memory than this machine provides. Verify memory before starting the R build; retain the verified downloaded counts and record `BUILD_FAILED` if the build cannot finish. Do not silently subsample cells or genes. Runtime data under `data/` are ignored by Git; only this run logic belongs in the repository history.
