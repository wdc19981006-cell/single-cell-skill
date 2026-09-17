# Real-world discovery cases

Development policy: inspect metadata and format links, use tiny synthetic fixtures for regressions, and request explicit user grouping before any real sample end-to-end run. Do not download all four studies or hard-code their groups. Discovery results below were retrieved from NCBI SOFT on 2026-09-13; initially two GSM records per study were sampled, not a full biological manifest.

| Dataset | Observed discovery | Test purpose |
|---|---|---|
| [GSE231993](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231993) | 12 GSM; ulcerative colitis study; GSM7307094 and GSM7307095 each publish barcodes/features/matrix gzip files | Multi-sample tissue 10x trio, automatic canonical destination proposals |
| [GSE202051](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202051) | 74 GSM; pancreatic cancer molecular taxonomy; three Series H5AD.gz files including adata_010nuc_10x and organoid-related data | snRNA-seq H5AD routing; inspect study heterogeneity and counts layers before selecting input; never Read10X_h5 |
| [GSE211644](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211644) | 50 GSM; pancreatic tumor lymphocyte study; fresh/grown matrices, genes, barcodes and separate metadata.csv.gz | Complex pooled cell/sample and clinical mappings; distinguish fresh tissue from cultured TIL; never infer patient IDs from names |
| [GSE229413](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229413) | 32 GSM; pancreatic donor/neoplastic-lesion study; GSM7162998 and GSM7162999 each publish raw and filtered feature_bc_matrix.h5 | Filtered preference within a known GSM; verify schema before H5 reading |
| GSE264203 local regression | The downloaded pooled 10x H5 has seven barcode suffix groups while the confirmed Series has six GSM; no exhaustive source-backed cell map is available in the local run | Correct STOP before Seurat construction; never pair suffixes with GSM by order |

GEO Series `RAW.tar` may contain processed files, so do not classify it as FASTQ-only by its filename. H5AD filenames establish routing candidates, not raw-count availability. GSE202051 is not uniformly one specimen/source type; inspect actual GSM and cell metadata.

Run metadata-only discovery:

```text
<python-executable> .agents/skills/geo-single-cell-loader/scripts/inspect_geo.py GSE231993 --fallback-reason "geo MCP real call failed: <actual error>"
```

Repeat for the other accessions when needed. A `--max-samples 2` development probe is incomplete and isolated under `data/<GSE>/.workflow/development/`; it cannot replace the full inventory or authorize a manifest. Complete discovery writes `.workflow/inspection.json`, `.workflow/sample_report.csv` and the visible `sample_info.txt`. The report has no group column; fill biological fields only from evidence and refresh the TXT before presenting the complete report to a user.

See [actual test report](../../../../../tests/validation-report.md) for executed tests and limitations. The synthetic generator creates 250-gene/4-cell count matrices with reused barcodes across samples, a deliberately normalized H5AD X plus integer counts layer/raw, and unsorted simulated sample/group mapping.
