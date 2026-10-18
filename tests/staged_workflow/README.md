# Staged workflow acceptance evidence

These are actual text records from R 4.5.0 validation on 2026-10-18.
See `../../docs/staged-workflow-validation.md` for reproduction and limits.

- `GSE337566`: legacy download reporting and nine unchanged protected files.
- `canonical_download`: isolated synthetic new-layout download regression.
- `GSE167297`: full existing 27,073-cell cohort, three QC steps, fresh PCA,
  unintegrated cluster review, 17 object checks, restart and plot checks.
- `analysis-smoke-v1`: explicitly synthetic, actual Harmony and test-only
  resolution finalization; this does not select a real cohort's resolution.
- `analysis-sensitivity-v1`: measured risk gates and corrected-assay PCA.
- `r_regressions.*`: format routing, original QC algorithms, entrypoints and
  artifact recovery. Logs preserve earlier failed attempts; the JSON records
  the final successful run.
- `python_tests.log`: final Python unit suite output.

Text copies use UTF-8/LF and trim log-line trailing padding for stable Git
evidence hashes; hashes inside runtime
records still refer to the original artifacts before this text normalization.
Paths in copied logs refer to the original local acceptance roots. Binary
objects/plots remain local and are represented by hashes and rendering checks.
`GSE167297/test_results.json` describes the successful resume starting at QC-3;
QC-1/QC-2 computation and recovery are recorded in the stage logs and state.
Its artifact hashes are the snapshot at the end of that run. Later audit
versions have their own `final_audit/audit_inputs.json` manifest.

All final audit rows retain source provenance; cell-level expression and
metadata tables are not uploaded here. The final report records missing
capture provenance, absent confirmed batch information, an unavailable marker,
and the user's pending resolution choice. No new real GEO download was run.
