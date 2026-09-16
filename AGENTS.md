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
