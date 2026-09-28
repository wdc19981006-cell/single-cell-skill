# single-cell-skill repository rules

For every real GSE run, follow
`.agents/skills/geo-single-cell-loader/references/audit-policy.md`.
Start the audit before discovery and use `run_confirmed.py` for confirmed
download/build/validation. Push success, failure and catchable interruption
audits to `wdc19981006-cell/single-cell-skill-audit/<GSE>/<run-id>/`.
The audit repository contains only the ten specified text/JSON/CSV files.
After successful upload, keep locally only raw inputs, the final RDS when
built, sample_info.txt, group_confirmation.json, and cell maps needed to
reconstruct the object. Retain temporary workflow files if upload fails.

During a real run, record issues and avoid broad restructuring. Subsequent
fixes must target a data-structure type, include a regression test, preserve
validated routes, and be pushed to the main single-cell-skill repository.
Input routing decisions must depend on verified structure rather than a
specific GSE accession.

## R runtime

All R computation, including temporary Codex checks, goes through
`runtime/r45/run_r45.py`. Never invoke `Rscript` or `Rscript.exe` directly,
including `Rscript -e`. Use `run_r45.py --expr "sessionInfo()"` instead.
The only environment definition is `runtime/r45/config.py`. Do not select an
R from PATH, change system R configuration or install/update/remove packages.
Run `run_r45.py --healthcheck` once per independent raw workflow, and use
`--route qc` for QC. Missing optional packages stop only the applicable route.

`qc/run_workflow.py` runs independent QC through the same runtime. Preserve
user decisions and the existing biological thresholds. Runtime crashes remain
failures even if stdout contains PASS. Keep completed artifacts; validate them
before considering another build. Never overwrite a final RDS.

An already completed artifact can be checked read-only with
`run_r45.py .agents/skills/geo-single-cell-loader/scripts/validate_seurat.R
--existing <root> <GSE>`. This regression check logs to `.runtime/`; it does
not initiate a new GEO discovery/download/build run.
