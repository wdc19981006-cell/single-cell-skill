"""Diagnose saved stage evidence without rerunning any expression analysis."""
import argparse
import csv
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "runtime/r45"))
from runtime import run_r45
from workflow import digest

STEPS = {
    "download": ("01_download", ("download_summary.csv",)),
    "qc1": ("02_QC/01_basic_QC", ("sample_QC_summary.csv", "gene_filter_summary.csv")),
    "qc2": ("02_QC/02_doublet_cellcycle", ("doublet_summary.csv", "cell_cycle_summary.csv")),
    "qc3": ("02_QC/03_ambient_RNA", ("decontX_summary.csv", "contamination_by_sample.csv")),
    "pca": ("03_analysis/01_pca_regression", ("pca_selection.csv", "regression_diagnostics.csv")),
    "harmony": ("03_analysis/02_harmony", ("batch_design_check.csv", "integration_metrics.csv", "marker_comparison.csv")),
    "clustering": ("03_analysis/03_clustering", ("cluster_summary.csv",)),
}
FIELDS = ("stage", "source_file", "source_sha256", "row", "sample", "metric", "value", "status")
WARNING_FIELDS = ("stage", "category", "severity", "message", "recommendation", "source_file")


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as source:
        return list(csv.DictReader(source))


def write_csv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def safe_child(root, relative):
    target = root / relative
    if target.resolve() != target or not target.resolve().is_relative_to(root):
        raise ValueError("Audit source directory redirects or escapes the project")
    return target


def collect(root, gse):
    """Keep absent evidence explicit; never infer scientific success from files."""
    if not re.fullmatch(r"GSE\d+", gse):
        raise ValueError("Expected GSE accession")
    root = Path(root).resolve()
    project = safe_child(root, f"results/{gse}")
    audit, warnings, inputs = [], [], {}
    state_sources = set()
    def record_state(stage, state_path):
        if not state_path.is_file() or state_path in state_sources:
            return
        state_sources.add(state_path)
        state_path = safe_child(state_path.parent, state_path.name)
        state = json.loads(state_path.read_text(encoding="utf-8"))
        sha = inputs[str(state_path)] = digest(state_path)
        for key, value in state.items():
            if not isinstance(value, str) or key == "workflow_audit" or key.endswith(("sha256", "updated_at")):
                continue
            target_stage = {"cluster-review": "clustering", "finalize": "clustering"}.get(key, key)
            target_stage = target_stage if target_stage in STEPS else stage
            audit.append(dict(stage=target_stage, source_file=str(state_path), source_sha256=sha,
                row=1, sample="ALL", metric=key, value=value, status="RECORDED_STATE"))
            if any(word in value.upper() for word in ("NEEDS_", "SKIPPED", "UNEVALUATED", "UNVERIFIED", "PARTIAL", "FAILED", "ERROR")):
                category = ("BIOLOGICAL_DECISION" if "NEEDS_" in value.upper() else
                            "TECHNICAL_ERROR" if any(word in value.upper() for word in ("FAILED", "ERROR")) else
                            "DESIGN_LIMITATION" if "SKIPPED" in value.upper() else "REQUIRES_REVIEW")
                warnings.append(dict(stage=target_stage, category=category, severity="error" if category=="TECHNICAL_ERROR" else "review",
                    message=f"{key}: {value}",
                    recommendation="Review the saved evidence and obtain the recorded user decision where required; preserve completed results.",
                    source_file=str(state_path)))
        if state.get("error"):
            error = state["error"]
            decision = "NEEDS_USER_DECISION" in json.dumps(error).upper()
            warnings.append(dict(stage=stage, category="BIOLOGICAL_DECISION" if decision else "TECHNICAL_ERROR",
                severity="review" if decision else "error", message=json.dumps(error, ensure_ascii=False),
                recommendation="Resolve the reported decision or concrete error and resume only that step; preserve completed RDS files.",
                source_file=str(state_path)))
    for stage, (relative, required) in STEPS.items():
        area = safe_child(project, relative)
        if stage == "download" and not area.exists():
            area = safe_child(root, f"data/{gse}")
        present = list(area.glob("*.csv")) if area.is_dir() else []
        # Cell-level tables belong to their step; audit aggregates already saved
        # summaries instead of duplicating potentially millions of records.
        present = sorted(p for p in present if p.name.lower() != "metadata.csv" and not any(s in p.name.lower() for s in ("cell_qc_all", "cell_level", "per_cell", "cell_contamination", "decontx_cells", "doublet_cells", "cell_cycle_cells", "cluster_assignments")))
        for name in required:
            if not (area / name).is_file():
                warnings.append(dict(stage=stage, category="MISSING_EVIDENCE", severity="review",
                    message=f"{name} is unavailable; this component is not validated by this audit.",
                    recommendation="Run only the missing step or supply its original evidence; do not rerun completed predecessors.",
                    source_file=str(area / name)))
        for path in present:
            path = safe_child(area, path.name)
            sha = digest(path)
            rows = read_csv(path)
            inputs[str(path)] = sha
            if not rows and not any(word in path.name.lower() for word in ("warning", "risk")):
                warnings.append(dict(stage=stage, category="EMPTY_EVIDENCE", severity="review",
                    message=f"{path.name} has no data records; this component is not validated by this audit.",
                    recommendation="Inspect the corresponding step log and resume only its missing output; preserve completed results.",
                    source_file=str(path)))
            if path.name == "gene_filter_summary.csv" and rows and "retained" in rows[0]:
                retained = sum(r["retained"].upper() == "TRUE" for r in rows)
                rows = [dict(input_genes=len(rows), retained_genes=retained,
                             removed_genes=len(rows)-retained,
                             gene_retention=retained/len(rows))]
            for index, row in enumerate(rows, 1):
                sample = row.get("sample", row.get("library_id", "ALL"))
                if any(word in path.name.lower() for word in ("warning", "risk")):
                    warnings.append(dict(stage=stage, category=row.get("category", "PARAMETER_UNCERTAINTY"),
                        severity=row.get("severity", row.get("status", "review")),
                        message=row.get("message", row.get("note", str(row))),
                        recommendation=row.get("recommendation", "Review the cited evidence before choosing biological parameters."),
                        source_file=str(path)))
                elif any(word in row.get("status", "").upper() for word in ("WARNING", "FAILED", "NEEDS_", "NOT_EVALUATED", "UNVERIFIED")):
                    warnings.append(dict(stage=stage, category=("TECHNICAL_ERROR" if "FAILED" in row.get("status", "").upper() else "REQUIRES_REVIEW"),
                        severity="review", message=json.dumps(row,ensure_ascii=False),
                        recommendation="Review the reported limitation; preserve cells and sample identities until the relevant decision is supported.",
                        source_file=str(path)))
                if "metric" in row and "value" in row:
                    values = [(row["metric"], row["value"])]
                else:
                    values = [(k, v) for k, v in row.items() if k not in ("sample", "library_id")]
                for key, value in values:
                    audit.append(dict(stage=stage, source_file=str(path), source_sha256=sha,
                        row=index, sample=sample, metric=key, value=value, status=row.get("status", "OBSERVED")))
        # Carry failures and user-decision gates into the unified report.
        for state_path in (area / "state.json", area / "analysis_state.json", area / ".workflow/state.json"):
            record_state(stage, state_path)
    for stage, relative in (("qc", "02_QC/.workflow/state.json"),
                            ("analysis", "03_analysis/.workflow/state.json")):
        record_state(stage, safe_child(project, relative))
    return audit, warnings, inputs


def run(gse, root=REPO):
    root = Path(root).resolve()
    rows, warnings, inputs = collect(root, gse)
    area = safe_child(root, f"results/{gse}/03_analysis/99_workflow_audit")
    area.mkdir(parents=True, exist_ok=True)
    semantic_inputs = dict(inputs)
    semantic_inputs["report_renderer_sha256"] = digest(REPO / "workflow_audit/render_audit.R")
    for path in inputs:
        if Path(path).name in ("state.json", "analysis_state.json"):
            state = json.loads(Path(path).read_text(encoding="utf-8"))
            state.pop("updated_at", None)
            state.pop("workflow_audit", None)
            semantic_inputs[path] = hashlib.sha256(json.dumps(state, sort_keys=True).encode()).hexdigest()
    identity = hashlib.sha256(json.dumps(semantic_inputs, sort_keys=True).encode()).hexdigest()
    expected = ("workflow_audit.csv", "workflow_warnings.csv", "workflow_audit.pdf", "audit.log")
    state_path = area / "audit_inputs.json"
    # Input changes create immutable report versions. Repeated checks of the
    # same version must reuse that report instead of making another PDF.
    candidates = [area] + sorted((p for p in area.iterdir() if p.is_dir()), reverse=True)
    for candidate in candidates:
        candidate = safe_child(area, candidate.relative_to(area)) if candidate != area else area
        saved = candidate / "audit_inputs.json"
        previous = json.loads(saved.read_text(encoding="utf-8")) if saved.exists() else None
        if previous and previous.get("identity") == identity:
            if not all((candidate / p).is_file() and digest(candidate / p) == previous.get("outputs", {}).get(p) for p in expected):
                raise ValueError("Saved audit output changed; preserve it and use an isolated root for a new report")
            return candidate
    if any((area / name).exists() for name in expected):
        # An updated audit is a new version; never replace a validated report.
        area = area / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        area.mkdir()
        state_path = area / "audit_inputs.json"
    write_csv(area / "workflow_audit.csv", rows, FIELDS)
    write_csv(area / "workflow_warnings.csv", warnings, WARNING_FIELDS)
    with (area / "audit.log").open("w", encoding="utf-8") as log:
        log.write(f"UTC {datetime.now(timezone.utc).isoformat()}\nInput files: {len(inputs)}\nObserved fields: {len(rows)}\nWarnings: {len(warnings)}\n")
        log.write("Read existing stage summaries only. No normalization, filtering, integration or clustering was rerun.\n")
        run_r45(REPO / "workflow_audit/render_audit.R", [area, gse], stage="workflow_audit", cwd=REPO, log_path=log)
    for path, sha in inputs.items():
        if digest(path) != sha:
            raise ValueError(f"Input evidence changed during the audit: {path}")
    outputs = {name: digest(area / name) for name in expected}
    state_path.write_text(json.dumps(dict(identity=identity, inputs=inputs, outputs=outputs,
        created_at=datetime.now(timezone.utc).isoformat()), indent=2), encoding="utf-8")
    return area


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gse")
    parser.add_argument("--root", type=Path, default=REPO)
    args = parser.parse_args()
    print(run(args.gse, args.root))
