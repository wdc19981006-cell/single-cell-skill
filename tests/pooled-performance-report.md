# MCP-first and pooled-input performance report

Validation date: 2026-09-15/16 (Asia/Shanghai)

## GEO MCP validation

| Check | Result |
|---|---|
| MCP available in the active Codex session | YES |
| `get_geo_info("GSE181919")` real tool call | PASS |
| `list_geo_files("GSE181919")` real tool call | PASS |
| Requested/returned accession | GSE181919 / GSE181919 |
| `sample_count` | 37 |
| Returned GSM records | 37 |
| Public files | 2 |
| `complete` | true |
| `get_geo_sample` calls | 0 (not needed) |
| Fallback used | NO |

The file inventory included the Series-level pooled `GSE181919_UMI_counts.txt.gz` expression matrix and `GSE181919_Barcode_metadata.txt.gz` cell/barcode metadata. Both entries had explicit `owner_accession=GSE181919`, `source_level=GSE`, and `candidate_format=txt`.

The desktop session exposed and successfully called the real global `geo` tools. A separate CLI diagnostic was inconsistent: `codex mcp list` reported no configured servers and `codex mcp get geo` reported not found, while `C:\Users\22241\.codex\config.toml` contained the expected `geo` command/args and both configured files existed. Real tool calls are the functional verification; the global MCP server project was not modified.

## Stage A policy

| Policy | Result |
|---|---|
| Every new Stage A discovery is MCP-first | YES |
| Known GSE call order | `get_geo_info`, then `list_geo_files` |
| Keyword call order | `search_geo`, then `get_geo_info`, then `list_geo_files` |
| Fallback only after an actual MCP failure/incomplete result | YES |
| Duplicate fallback after MCP success | NO |
| Fallback reason required and recorded | YES |

Future `inspection.json` files record `source`, `discovery_method`, `mcp_tools_used`, `mcp_complete`, discovery timestamps, and `discovery_seconds`; `fallback_reason` exists only for fallback. The existing confirmed GSE181919 workflow was not silently regenerated during MCP verification.

## GSE181919 pooled build

Baseline timing: unavailable (the pre-fix workflow had no reliable stage timing log).

| Metric | Result |
|---|---:|
| Manifest rows | 37 |
| Unique physical expression inputs | 1 |
| Shared samples | 37 |
| Expression file | `data/GSE181919/raw/GSE181919_UMI_counts.txt.gz` |
| Expression bytes | 127,878,601 |
| Expression SHA256 | `cc5e83a48afa4ec48f2b2f64d449354124cb54976c7a9755f96a759543e5c6d6` |
| Old-loop theoretical expression reads | 37 |
| Actual expression reads after fix | 1 |
| Actual cell-map reads after fix | 1 |
| Reader | `data.table::fread` |
| Manifest read | 0.080 s |
| Expression read | 46.610 s |
| Cell mapping | 0.260 s |
| Feature alignment | 6.160 s |
| Matrix combine | 4.870 s |
| CreateSeuratObject | 11.490 s |
| Metadata mapping | 0.120 s |
| Validation | 0.310 s |
| Serialization and reload | 64.570 s |
| Total build | 136.700 s |
| Raw redownloaded | NO |
| Final cells | 54,239 |
| Final genes | 20,000 |
| Final RDS bytes | 279,066,168 |
| Final RDS SHA256 | `d0ce5e6fcf4a750618cf750ae484a58e6e0f897d48ae2284e11f39267dc4f603` |

The raw file's mtime, byte count, and SHA256 were unchanged when `download_processed.py` verified and reused it. The final RDS SHA256 is identical to the archived pre-fix RDS SHA256. The builder reads the pooled sparse matrix once, maps all columns in place using the exhaustive cell map, prefixes cells once, and never splits the matrix into 37 sample matrices for reassembly.

### Cells per sample

| Sample | Cells | Sample | Cells |
|---|---:|---|---:|
| GSM5514352 | 608 | GSM5514371 | 1,260 |
| GSM5514353 | 1,731 | GSM5514372 | 1,656 |
| GSM5514354 | 697 | GSM5514373 | 2,612 |
| GSM5514355 | 1,143 | GSM5514374 | 2,676 |
| GSM5514356 | 944 | GSM5514375 | 2,695 |
| GSM5514357 | 282 | GSM5514376 | 2,435 |
| GSM5514358 | 523 | GSM5514377 | 2,399 |
| GSM5514359 | 1,049 | GSM5514378 | 951 |
| GSM5514360 | 1,662 | GSM5514379 | 2,228 |
| GSM5514361 | 1,413 | GSM5514380 | 1,186 |
| GSM5514362 | 1,801 | GSM5514381 | 2,102 |
| GSM5514363 | 1,104 | GSM5514382 | 712 |
| GSM5514364 | 2,350 | GSM5514383 | 1,712 |
| GSM5514365 | 1,327 | GSM5514384 | 2,228 |
| GSM5514366 | 536 | GSM5514385 | 2,238 |
| GSM5514367 | 1,206 | GSM5514386 | 1,238 |
| GSM5514368 | 830 | GSM5514387 | 823 |
| GSM5514369 | 1,084 | GSM5514388 | 1,948 |
| GSM5514370 | 850 |  |  |

Groups: CA 23,088; LN 8,204; LP 6,527; NL 16,420.

## GSE231993 regression

Status: **NOT VERIFIED**. `data/GSE231993/` was not present locally, so no data were downloaded merely to perform this conditional non-regression check. Reference targets remain: 60,665 cells; 25,953 genes; 305,606,810 total UMI; 58,767,679 nonzero entries, with per-sample counts documented in the project history.

## Tests

| Suite | Result |
|---|---|
| Python unit tests | 30/30 PASS |
| R test files | 4/4 PASS |
| Pooled regression | PASS (reader calls 1; cell-map calls 1; all required error cases) |
| Global `min.cells` regression | PASS |
| Synthetic end-to-end | PASS |
| Stage A MCP policy mocks | PASS (success: fallback 0; failure: fallback 1) |

The test totals above are refreshed by the final pre-commit test run; failures must be reflected here before release.
