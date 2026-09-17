# Pre-build probes for text and pooled inputs

The confirmed manifest may leave `delimiter`, `orientation`, and
`feature_column` blank for a text input. After verified download and before
R, `run_confirmed.py` calls `prebuild_probe.py`. It reads only the header and
first eight expression rows of each unique CSV/TSV/TXT input. Explicit
configuration wins: a reliable contradictory observation stops the run.
The probe fills only blank technical fields, including whether the header
already contains an ID column. A named first gene column may retain an
explicit `__row_names__` streaming route; the R text reader uses the probe's
`text_header_missing_id` field to avoid inserting a second ID column. The
probe checks unique cell headers and
sampled nonnegative integer values, and records `text_schema_probe.json`.
The full R reader still validates all counts.

For a pooled input with no `cell_map_path`, inspect in this order:

1. Matrix barcode/library IDs, including suffix groups.
2. GEO supplementary cell/sample metadata.
3. Linked publication supplements.
4. The author's public repository.

Record what was actually checked in temporary
`.workflow/mapping_evidence.json`. Use one entry per source category, with
`source_type`, public `url` when present, `result`, and a concrete
`reason` if no URL or map exists. Use `geo_supplementary`,
`publication_supplement`, and `author_repository` as source types.
A candidate map is a small, publicly sourced CSV under
`.workflow/mapping_candidates/` with explicit cell and sample columns:

```json
{
  "checks": [
    {"source_type":"geo_supplementary","url":"https://example.org/cells.csv","result":"candidate found"},
    {"source_type":"publication_supplement","url":"https://example.org/supplement","result":"no matching map"},
    {"source_type":"author_repository","reason":"No public repository link found","result":"not_found"}
  ],
  "candidate_maps": [
    {"source_type":"geo_supplementary","source_url":"https://example.org/cells.csv",
     "local_path":"data/GSE123/.workflow/mapping_candidates/cells.csv",
     "cell_column":"cell","sample_column":"sample"}
  ]
}
```

An author sample alias may be converted only with an explicit
`alias_to_sample` map and `alias_evidence` public citation. The candidate
must contain each matrix cell exactly once and cover exactly the manifest
samples. A mismatch between barcode suffix-group count and GSM count stops
the missing-map route; suffix order never establishes sample identity.
Successful evidence creates one canonical temporary cell map and fills the
manifest `cell_map_path` and MD5. The audit trace includes the probe result.

The technical completion script verifies the original confirmation receipt
before editing and preserves the exact user group map and sample rows. It
refreshes the manifest checksum only for verified technical fields. Manual
manifest changes still require new user confirmation.
