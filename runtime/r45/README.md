# R 4.5.0 runtime

`run_r45.py` is the only supported R process launcher. It uses the fixed
installation in `config.py`, a process-local startup profile, `CLI_NO_THREAD=1`,
and `LC_ALL=C`. It records stdout, stderr, exit code and environment in ignored
`.runtime/r45/runs/` files. Native crashes remain failures.

Use `/c/Python312/python.exe D:/CodexProjects/single-cell-skill/runtime/r45/run_r45.py --healthcheck` before a raw
build, or add `--route qc`, `text`, `10x_h5`, or `h5ad_native` to check only the
extra dependencies for that route. `--expr "sessionInfo()"` is available for
one-off diagnostics. Python workflow code imports `runtime.py::run_r45`.

The profile checks both the exact version from `config.py` and normalized
`R.home()`. It rejects conflicting libraries or error handlers. If the
architecture environment variable is missing, it obtains the real value from
R; ARM64 is selected only during exit cleanup. Site/user environment files and
the site profile are disabled for this subprocess. No global configuration,
package installation, numerical method or Windows Error Reporting setting is
changed. `--vanilla` and `--no-init-file` are rejected.

`health.json` is written from actual R output after the process exits. Core
packages are required for every R workflow; optional dependencies are checked
only for the requested route. Raw workflows and `qc/run_workflow.py` each run
one health check, rather than trusting a stale PASS from another process.

Each run has a unique JSON record and stdout/stderr files. The initial record
is written before R starts, then updated with its original exit code and
duration. The helper uses a per-call record, so concurrent callers do not
infer status from `latest.json`. A catchable interruption stops this child
process tree and records exit 130; it never turns a crash into success.

Raw and QC checkpoints are separate. A final RDS is never overwritten.
Existing raw artifacts receive independent validation; existing QC artifacts
also require a matching completion report and PDF/PNG outputs. Raw inputs and
user decisions survive failures. Loader audit cleanup archives original
manifest/receipt bytes and checkpoint metadata inside `stage_a.json`; the
ten-file audit layout stays unchanged. Historical completed raw artifacts
without a manifest use `validate_seurat.R --existing ROOT GSE`, which checks
structure, confirmed groups and archived dimensions without writing the RDS.

Error records include category, stage, script, exit code, message and timestamp.
The runtime categories are `R_RUNTIME_MISSING`, `R_PACKAGE_MISSING`,
`R_NATIVE_CRASH`, and `R_SCRIPT_ERROR`. Workflow categories include
`DOWNLOAD_ERROR`, `CHECKSUM_ERROR`, `FORMAT_ERROR`, `MAPPING_ERROR`,
`RAW_BUILD_ERROR`, `RAW_VALIDATION_ERROR`, `QC_NEEDS_USER_DECISION`, and
`QC_ERROR`. `RAW_ARTIFACT_CREATED_BUT_RUNTIME_FAILED` records a preserved
artifact without claiming runtime success.

Run `tests/runtime/run_tests.py --real-fixture` for runtime acceptance and the
read-only GSE225857 fixture; `tests/run_regressions.py` creates a separate
synthetic workspace and runs the existing Loader and QC tests.
