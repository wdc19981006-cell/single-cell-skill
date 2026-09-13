# Manifest and metadata contract

One CSV row per stable expression sample; one GSE per manifest. UTF-8, header present, comma quoting supported. No duplicate sample rows, even for three-file inputs. Required final cell metadata:

| Field | Meaning |
|---|---|
| database | GEO Series accession, e.g. GSE231993 |
| sample | GSM preferred; otherwise publicly documented stable matrix/sample ID |
| tissue | Anatomy only: Pancreas, Colon, Blood, Liver, Lung, Bone_Marrow |
| disease | Donor disease background, e.g. Healthy, PDAC, UC; uninvolved PDAC tissue remains PDAC |
| source_type | Tissue, PBMC, Whole_Blood, Bone_Marrow, Ascites, Pleural_Fluid, Organoid, Cell_Line, Cultured_TIL |
| group | User-confirmed analysis label; never derived automatically from disease |

Required workflow columns: `local_path`, `file_type`, `count_source`, `count_evidence`, `metadata_evidence`. Evidence must identify a public URL/table/record and relevant field/statement; review semantic validity, not just whether a URL is present. Paths are relative to repository root, use forward slashes and stay inside `datasets/<GSE>/`.

Optional `patient`, `specimen`, `cohort`, `treatment` require corresponding `<field>_evidence` for each non-NA value. Omit absent columns; partially documented columns may contain NA, reported in run summary. Never infer patient identity from P01, T01, donor3 or filenames. `author_sample` and `sample_description` are report/manifest helper fields, not mandatory cell columns.

## Download plan

`files_json` is a JSON array inside one CSV cell. Example:

```json
[{"url":"https://ftp.ncbi.nlm.nih.gov/.../sample_matrix.mtx.gz","local_path":"datasets/GSE123/GSM123/matrix.mtx.gz"}]
```

For a trio, put three entries in the same array, with canonical matrix/features/barcodes destinations. For an archive member add `member` with the exact archive member name. Optional `sha256` means the downloaded file/archive SHA256, not extracted-member SHA256. One archive URL is downloaded once. `files_json=[]` permits existing local inputs, which are identified as such in the summary. URLs cannot contain credentials. Downloaded hashes are local integrity/provenance records, not author-signed checksums.

## Pooled matrices and separate metadata

When multiple samples share an input, every related manifest row must identify the same `cell_map_path`. This is a CSV with columns `cell,sample`; cell IDs exactly match the input matrix, without prefixes. Duplicate cells, extra cells, missing cells, missing samples or foreign samples stop loading. The manifest references this mapping as an explicit subordinate artifact; expression/clinical columns are never joined using position or filename guessing. Record the cell map source in `metadata_evidence`.

## Group confirmation

Save the user's reply in `manifests/<GSE>_group_confirmation.json`:

```json
{"confirmed_by":"user","user_statement":"Literal user grouping instruction","groups":{"GSM123":"UserChosenLabel"}}
```

`build_manifest.py` requires exact sample key coverage and nonempty groups and writes `<GSE>_sample_manifest.confirmation.json`. The receipt carries the actual statement, group map, UTC time and manifest MD5 (change detection only, not a security signature). A CSV edit invalidates it. Archive the previous manifest and rebuild from an explicitly reconfirmed report. Never claim the receipt proves identity or independent verification of public evidence.
