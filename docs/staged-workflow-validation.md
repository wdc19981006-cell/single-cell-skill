# Three-stage workflow implementation and validation

This change started from local and GitHub `main` commit
`74a9fbaf69914f4b964776482a2276cdd26ac9ed`. The worktree was clean. The initial
inspection found a working download/router/audit system, an independent
two-step QC pipeline under `data/<GSE>/qc`, a recently validated PCA rule,
and no DecontX or formal analysis stage. The old cell filter did not apply
the requested retained-cell gene expression rule.

## Implementation

New downloads use `results/<GSE>/01_download`. Existing legacy datasets,
confirmed manifests and receipts remain in place. A shared layout resolver
updates path validation, probes, build/validation and audit cleanup without
changing format-reader algorithms. Validated raw objects produce metadata,
sample counts, a genuine integrity PDF and a retained download log. The
private download audit still contains exactly its ten allowed text files.

QC outputs are separated into basic QC, DoubletFinder/cell-cycle, and
ambient RNA directories. QC-1 preserves default biological thresholds,
adds the >3 retained-cell gene rule, and recalculates current count metrics.
QC-2 preserves the reference PCA/DoubletFinder algorithm and filters using
DF_adj. Explicit capture-library metadata is honored; legacy sample-only
mapping is identified as unverified. DecontX uses original sparse RNA
counts and saves a separate `decontX` assay; it removes no cells by a
fixed contamination percentage. Each step has its own RDS, tables, plots,
configuration, log and independently validated restart path.

Formal analysis freshly computes PCA with pMT-only baseline regression.
Recorded PC/metadata relationships gate explicitly requested targeted
sensitivity comparisons. Harmony checks technical metadata and confounding;
biological group is forbidden as a correction variable. Resolution scans
share one UMAP. Only a recorded user-selected reviewed resolution creates
formal labels and `seurat_analysis.rds`.

The local scientific audit reads saved summaries and states, reports
missing evidence and biological decisions, and versions changed reports.
It excludes large cell-level tables and does not rerun scientific analysis.
R remains fixed at 4.5.0, with no package installations or updates.

## Reproduction

Run commands from the repository in Git Bash. All R calculations use the
fixed launcher. The Python suite and synthetic R regression are:

```bash
/c/Python312/python.exe -X utf8 -m unittest discover -s tests -p 'test_*.py'
/c/Python312/python.exe D:/CodexProjects/single-cell-skill/tests/run_regressions.py
```

The real downstream acceptance runner copies the complete existing raw
object and the exact historical user decision into an isolated test root:

```bash
/c/Python312/python.exe D:/CodexProjects/single-cell-skill/tests/run_staged_acceptance.py GSE167297 \
  --raw-source D:/CodexProjects/single-cell-skill/.runtime/qc-pca-fixed/data/GSE167297/seurat_raw.rds \
  --decision-source D:/CodexProjects/single-cell-skill/.runtime/qc-pca-fixed/data/GSE167297/qc/qc_decision.json \
  --root D:/CodexProjects/single-cell-skill/.runtime/staged-full-20261018
```

Use a new isolated root for a new computation. Reusing the same root validates
completed artifacts and resumes missing steps. No new GEO discovery or
download is needed for this test. The real cluster resolution remains a user
decision. The separate synthetic analysis test exercises actual Harmony,
test-only finalization, audit and no-recompute resume:

```bash
/c/Python312/python.exe D:/CodexProjects/single-cell-skill/tests/run_analysis_smoke.py \
  --output D:/CodexProjects/single-cell-skill/.runtime/new-analysis-smoke
```

## Evidence and limits

The exact source hashes, executed stages, runtime exits, object checks and
figure checks are saved in [the test evidence directory](../tests/staged_workflow/).
Large RDS/matrix files and generated plot binaries remain local in ignored
data/result directories. Failed attempts are retained in the runtime logs;
they are not counted as successful computations.

GSE337566 separately verifies legacy stage-one reporting: two samples,
14,289 cells and 21,506 features; nine protected source/RDS/confirmation/TXT
files had identical SHA256 before and after. Canonical synthetic downloads
and legacy 10x/H5/H5AD/text/pooled routes are tested separately.

GSE167297 uses the full existing cohort of 27,073 cells and 22,026 genes.
QC-1 retains 23,243 cells and 20,555 genes; QC-2 retains 22,420 cells after
823 DF_adj doublets. DecontX retains all 22,420 cells and original RNA
counts, alongside the corrected assay. It estimates mean contamination
0.085805 and median 0.013430; these are measurements, not removal thresholds.

The historical input lacks independently verified capture-library metadata;
this regression tests its existing per-sample route and records that limit.
It also lacks a confirmed independent technical batch, so the real-cohort
default retains unintegrated PCA. The crossed-batch synthetic test exercises
the installed Harmony implementation. No real-cohort final resolution is
invented for a software acceptance test.

The first full QC-2 process saved its artifact but failed during its closing
code path (`object 'gse' not found`). This attempt remained a failure. A later
independent read-only validation successfully recovered the RDS without
rerunning DoubletFinder. Current entrypoint, recovery and native-crash tests
exercise the corrected behavior. Historical failure logs remain available.

## Final acceptance on 2026-10-18

- Python: 120 tests passed; R regression runner: all 15 recorded checks passed.
- GSE167297: 17 independent object invariants passed. Fresh formal PCA selected
  21 PCs, and nine resolutions (0.1 through 1.2) were reviewed on one UMAP.
- QC plots: 9 actual PDF pages and 9 PNG files passed rendering/DPI checks.
  Formal analysis: 6 PDF plots and 6 PNG files passed. PNGs are 300 dpi.
- The final audit is a new 28-page version, rendered and visually checked;
  prior 32/35-page snapshots remain local for provenance. The latest report
  contains capture-mapping limits and the resolution decision, with no stale
  current runtime failure. Long tokens wrap inside the page, and PC/rho/risk
  records stay together. The complete observations remain in CSV.
- Source raw RDS and literal historical user decision hashes are unchanged.
  GSE167297's raw-only resume fixture intentionally lacks a download summary;
  stage-one acceptance is separately established with GSE337566 and the
  canonical synthetic download test. This is not a new GEO download run.
- PECAM1 is absent from the retained real-cohort matrix; marker diagnostics
  explicitly report this and require tissue-relevant marker review.

See [changed files](../tests/staged_workflow/changed_files.txt),
[evidence hashes](../tests/staged_workflow/evidence_manifest.json), and
[real object checks](../tests/staged_workflow/GSE167297/object_validation.csv).
The software implementation is validated through the review gate; real
biological finalization remains pending the user's resolution selection.
