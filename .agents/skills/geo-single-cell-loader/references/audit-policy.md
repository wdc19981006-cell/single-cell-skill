# Real GSE run audit policy

This is a permanent rule for every real GEO accession run. Keep `.workflow/` only while discovery, grouping, download, build or audit upload is in progress. After a terminal run and successful push, the local `data/<GSE>/` directory may contain only `raw/`, `seurat_raw.rds` if built, `sample_info.txt`, `group_confirmation.json` if confirmed, and `cell_map/` only when needed to reconstruct the object. Preserve these files. Existing historical workflows require a separate verified migration; never delete them merely because this rule changed.

The independent private GitHub repository is `wdc19981006-cell/single-cell-skill-audit`. Each real run gets one immutable `<GSE>/<run-id>/` directory. Start the run before Stage A with `audit_run.py start GSE...`. The confirmed Stage B must use `run_confirmed.py GSE...`, which pushes success, failure or interruption audits automatically. If Stage A itself fails or is interrupted, finalize manually with `audit_run.py finalize GSE... --status failure|interrupted --reason "observed cause"`. A workflow waiting for the user's grouping is active, so keep it until the user responds or cancels.

The audit directory contains exactly:

- `decision_trace.md`
- `stage_a.json`
- `file_selection.csv`
- `download_profile.csv`
- `input_routing.csv`
- `build_profile.json`
- `validation.json`
- `sample_summary.csv`
- `run_summary.txt`
- `execution.log`

`decision_trace.md` records observable decisions only: discovered files, selected expression files, excluded candidates with reasons, raw-count evidence, structure, sample/GSM/cell mapping, reader rationale, actual MCP/fallback behavior, exceptions and stopping reasons. Do not store model hidden chain-of-thought. Add inspected evidence to `.workflow/decision_notes.md` and exact excluded URL reasons to `.workflow/file_selection_reasons.csv` before successful finalization. Do not invent reasons when evidence is missing; stop and review.

`download_profile.csv` records bytes, transfer seconds, MB/s (decimal), retries, REUSED/DOWNLOADED, SHA seconds and archive extraction seconds for each source. `input_routing.csv` records format, reader, unique physical inputs, expression and cell-map reads, dimensions, read seconds and dense conversion. `build_profile.json` records read, mapping, combine, CreateSeuratObject, validation, saveRDS, total time and obtainable memory estimates. A missing stage in a failed run is explicitly marked unavailable.

Only text/JSON/CSV audit evidence is allowed in the audit repository. Never upload expression matrices, RDS, H5, H5AD, source archives, raw cell maps or other large data. Push the audit at the end of every real terminal run without waiting for another user instruction. If push fails, retain local `.workflow/` for retry and report the failure. Once push succeeds, remove only temporary `.workflow/` files; never remove raw data, RDS, group confirmation or necessary cell maps. Existing raw files can be reused only after URL, size and SHA verification against provenance restored from prior audit runs.

During a real run, record newly discovered issues and stop when required; avoid broad restructuring mid-run. Later fixes must target a data-structure type, add a regression test, preserve already validated data behavior, and be pushed to the `single-cell-skill` main repository. Input Router should select an optimized route automatically when the verified structure matches a previously optimized type; never key a route only to a GSE accession.
