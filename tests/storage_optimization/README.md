# Storage lifecycle acceptance — 2026-10-20

`python_results.json` and `python_tests.log` record focused lifecycle, legacy
compatibility, analysis, QC, audit and output-path tests. `test_results.json`
records an actual synthetic R 4.5.0 / Seurat 5.4.0 end-to-end run. R fixtures use
explicit SYNTHETIC selections only; no real dataset resolution was selected.

Verified actual PCA, crossed-batch Harmony, shared-coordinate resolution review,
targeted count regression, pause/resume, final independent readback, exact RNA
and fractional decontX expression hashes, FeaturePlot, FindMarkers, annotation
field assignment, compaction, checkpoint cleanup, final no-op restart, and audit
resolution of cleaned dependencies. Seven intermediate RDS were cleaned. Read
the JSON for measured bytes; this small fixture does not estimate real-study savings.

The staged QC R regression (`SC_TEST_COMPACT_QC=1` through `run_r45.py`) also
passed, exercising fixture DoubletFinder/DecontX APIs, real cycle/scientific QC
guards and a separate compact QC-3 copy with exact original/corrected counts.
These fixtures do not claim new real DoubletFinder or DecontX method validation.
`r_validation.json` records the successful runtime exits and logs.

Local `inventory/` (gitignored) contains read-only file/size/dependency inventories
of old outputs. Unknown dependencies remain UNKNOWN. No old RDS was deleted,
moved, compacted, or deserialized for inventory. Code does not expose automatic
legacy cleanup; any later deletion requires explicit confirmation of a concrete list.
