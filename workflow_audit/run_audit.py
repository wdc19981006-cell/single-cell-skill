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
sys.path.insert(0, str(REPO / "runtime"))
from project_paths import storage_root, project_path, analysis_path, runtime_output, dataset_name
from storage import manifest_valid, Storage, CLEANED, read
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
WARNING_FIELDS = ("stage", "category", "severity", "message", "recommendation", "source_file", "status")


def issue_status(value):
    value = str(value).upper()
    if value.strip().startswith(("INFO ", "PASS ")): return "PASS"
    if value == "RESOLVED": return "RESOLVED"
    if any(v in value for v in ("FAILED", "ERROR")): return "FAILED"
    if any(v in value for v in ("MISSING_EVIDENCE", "EMPTY_EVIDENCE")): return "MISSING_EVIDENCE"
    if any(v in value for v in ("NEEDS_", "ACTION_REQUIRED", "USER_REVIEW_REQUIRED", "BIOLOGICAL_DECISION")): return "ACTION_REQUIRED"
    return "WARNING"


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


def collect(root, gse, run_id=None):
    """Keep absent evidence explicit; never infer scientific success from files."""
    dataset_name(gse)
    root = Path(root).resolve()
    project = project_path(root, gse)
    analysis_area = analysis_path(root, gse, run_id)
    managed = (project / ".work/storage_policy.json").is_file()
    scope = "analysis/runs/"+run_id if run_id else "analysis/baseline"
    analysis_objects = project / ".work/checkpoints" / scope if managed else analysis_area
    analysis_state = analysis_objects / "state" if managed else analysis_area / ".workflow"
    qc_state = project / ".work/checkpoints/qc/state" if managed else project / "02_QC/.workflow"
    registry_path = project / ".work/checkpoints/registry.json"
    registry = read(registry_path).get("files",{}) if managed and registry_path.exists() else {}
    def cleaned_source(path):
        item = registry.get(str(path),{})
        return item.get("lifecycle")==CLEANED and Storage.valid_proof(item.get("successor",{}))
    def allowed(path):
        return path.resolve().is_relative_to(analysis_area) or (managed and path.resolve().is_relative_to(analysis_objects))
    audit, warnings, inputs = [], [], {}
    resolved = False
    final_manifest = analysis_state / "finalize.json"
    if final_manifest.is_file():
        inputs[str(final_manifest)] = digest(final_manifest)
        saved = json.loads(final_manifest.read_text(encoding="utf-8"))
        decision = saved.get("binding", {}).get("decision") or {}
        outputs = saved.get("outputs", {})
        source = Path(saved.get("binding", {}).get("input", ""))
        resolved = (saved.get("status") == "COMPLETE" and decision.get("confirmed_by") == "user"
                    and bool(decision.get("user_statement")) and bool(outputs)
                    and source == analysis_objects / "seurat_cluster_review.rds"
                    and ((source.is_file() and digest(source) == saved.get("binding", {}).get("input_sha256")) or cleaned_source(source))
                    and str(analysis_area / "seurat_analysis.rds") in outputs)
        for name, checksum in outputs.items():
            path = Path(name)
            valid = path.resolve().is_relative_to(analysis_area) and path.is_file() and digest(path) == checksum
            resolved = resolved and valid
            if path.is_file() and path.resolve().is_relative_to(analysis_area): inputs[str(path)] = digest(path)
        if source.is_file(): inputs[str(source)] = digest(source)
        if not resolved:
            warnings.append(dict(stage="clustering",category="FILE_STATE_INCONSISTENCY",severity="FAILED",
                message="Finalize manifest, input, user decision or published output checksums disagree",
                recommendation="Preserve outputs and inspect finalize validation; selection is not resolved",source_file=str(final_manifest)))
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
            target_stage = {"cluster-review": "clustering", "finalize": "clustering", "regression-review":"pca"}.get(key, key)
            target_stage = target_stage if target_stage in STEPS else stage
            audit.append(dict(stage=target_stage, source_file=str(state_path), source_sha256=sha,
                row=1, sample="ALL", metric=key, value=value, status="RECORDED_STATE"))
            if any(word in value.upper() for word in ("NEEDS_", "REVIEW_REQUIRED", "SKIPPED", "UNEVALUATED", "UNVERIFIED", "PARTIAL", "FAILED", "ERROR")):
                category = ("BIOLOGICAL_DECISION" if any(word in value.upper() for word in ("NEEDS_","REVIEW_REQUIRED")) else
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
        area = safe_child(analysis_area, relative.removeprefix("03_analysis/")) if relative.startswith("03_analysis/") else safe_child(project, relative)
        if stage == "download" and not area.exists():
            area = safe_child(root, f"data/{gse}")
        present = list(area.glob("*.csv")) if area.is_dir() else []
        # Cell-level tables belong to their step; audit aggregates already saved
        # summaries instead of duplicating potentially millions of records.
        present = sorted(p for p in present if p.name.lower() != "metadata.csv" and not any(s in p.name.lower() for s in ("cell_qc_all", "cell_level", "per_cell", "cell_contamination", "decontx_cells", "doublet_cells", "cell_cycle_cells", "cluster_assignments")))
        if stage=="pca":
            sensitivity = area / "regression_sensitivity"
            present.extend(p for p in sensitivity.rglob("*.csv") if p.name in
                ("sensitivity_summary.csv", "regression_review_warnings.csv"))
        if stage=="clustering" and resolved:
            selected=area / ("selected_resolution_"+format(float(decision["resolution"]),"g"))
            for name in ("final_cluster_sizes.csv","marker_availability.csv"):
                if (selected/name).is_file(): present.append(selected/name)
                else: warnings.append(dict(stage=stage,category="MISSING_EVIDENCE",severity="MISSING_EVIDENCE",
                    message=f"Selected resolution diagnostic {name} is missing",recommendation="Preserve the final RDS; generate only the missing selected-resolution review in an isolated output",source_file=str(selected/name)))
        for name in required:
            if not (area / name).is_file():
                warnings.append(dict(stage=stage, category="MISSING_EVIDENCE", severity="review",
                    message=f"{name} is unavailable; this component is not validated by this audit.",
                    recommendation="Run only the missing step or supply its original evidence; do not rerun completed predecessors.",
                    source_file=str(area / name)))
        for path in present:
            path = safe_child(area, path.relative_to(area))
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
                def diagnostic(category, message, severity="WARNING"):
                    warnings.append(dict(stage=stage,category=category,severity=severity,message=message,
                        recommendation="Review the cited diagnostic evidence; no automatic cell/sample deletion or regression adoption",source_file=str(path)))
                def number(key, default=0):
                    try: return float(row.get(key, default))
                    except (ValueError, TypeError): return default
                if path.name=="sample_QC_summary.csv" and number("input_cells")>0 and 1-number("retained_cells")/number("input_cells")>.15:
                    diagnostic("HIGH_QC_REMOVAL",f"{sample}: removed fraction {1-number('retained_cells')/number('input_cells'):.1%}")
                if path.name=="contamination_by_sample.csv" and (number("mean_contamination")>=.5 or number("median_contamination")>=.5 or number("p95_contamination",number("p90_contamination"))>=.8):
                    diagnostic("EXTREME_CONTAMINATION",f"{sample}: mean={row.get('mean_contamination')} median={row.get('median_contamination')} P95/P90={row.get('p95_contamination',row.get('p90_contamination'))}")
                if path.name=="regression_diagnostics.csv" and abs(number("rho"))>=.5 and row.get("variable")!="pMT":
                    diagnostic("REGRESSION_ASSOCIATION",f"{row.get('variable')} PC{row.get('pc')}: rho={row.get('rho')}")
                if path.name=="sensitivity_summary.csv" and row.get("recommendation")=="USER_REVIEW_REQUIRED":
                    diagnostic("REGRESSION_BIOLOGICAL_TRADEOFF",f"{row.get('target')}: {row.get('status')}, ARI={row.get('ari')}","ACTION_REQUIRED")
                if "unverified" in json.dumps(row).lower():
                    diagnostic("CAPTURE_MAPPING_UNVERIFIED",f"{sample}: capture-library evidence remains unverified")
                if path.name=="cluster_summary.csv" and (number("cells",100)<20 or number("fraction")>.8):
                    diagnostic("CLUSTER_SIZE_REVIEW",f"resolution={row.get('resolution')} cluster={row.get('cluster')} cells={row.get('cells')} fraction={row.get('fraction')}")
                if path.name=="final_cluster_sizes.csv" and (number("cells",100)<20 or number("fraction")>.8):
                    diagnostic("FINAL_CLUSTER_SIZE_REVIEW",f"selected resolution={decision['resolution']} cluster={row.get('cluster')} cells={row.get('cells')} fraction={row.get('fraction')}")
                if path.name=="marker_availability.csv" and row.get("status")=="MARKER_ABSENT":
                    diagnostic("MISSING_MARKER_DIAGNOSTIC",f"Selected-resolution marker absent: {row.get('marker')}")
                if any(word in path.name.lower() for word in ("warning", "risk")):
                    warnings.append(dict(stage=stage, category=row.get("category", "PARAMETER_UNCERTAINTY"),
                        severity=row.get("severity", row.get("status", "review")),
                        message=row.get("message", row.get("note", str(row))),
                        recommendation=row.get("recommendation", "Review the cited evidence before choosing biological parameters."),
                        source_file=str(path)))
                elif "MISSING" in row.get("status", "").upper():
                    warnings.append(dict(stage=stage,category="MISSING_EVIDENCE",severity="MISSING_EVIDENCE",
                        message=json.dumps(row,ensure_ascii=False),recommendation="Supply this missing diagnostic evidence; absence is not a pass",source_file=str(path)))
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
        if present and stage in ("pca","harmony","clustering"):
            proof=analysis_state/({"clustering":"cluster-review"}.get(stage,stage)+".json")
            if not proof.is_file(): warnings.append(dict(stage=stage,category="MISSING_EVIDENCE",severity="MISSING_EVIDENCE",
                message="Summary tables exist without a validated computation checkpoint",recommendation="Supply original validation/provenance; do not infer success from tables alone",source_file=str(proof)))
        if present and stage in ("qc1","qc2","qc3"):
            state_path=qc_state/"state.json"
            state=json.loads(state_path.read_text(encoding="utf-8")) if state_path.is_file() else {}
            artifact=area/f"seurat_QC-{stage[-1]}.rds"
            if managed and stage!="qc3": artifact=project/".work/checkpoints/qc"/artifact.name
            checksum=state.get(f"{stage}_sha256")
            if checksum and artifact.is_file():
                actual=inputs.get(str(artifact)) or digest(artifact);inputs[str(artifact)]=actual
                if actual!=checksum: warnings.append(dict(stage=stage,category="FILE_STATE_INCONSISTENCY",severity="FAILED",
                    message="QC artifact checksum differs from its validated checkpoint",recommendation="Preserve the artifact and inspect this stage's validation",source_file=str(artifact)))
                else: audit.append(dict(stage=stage,source_file=str(artifact),source_sha256=actual,row=1,sample="ALL",
                    metric="checkpoint_validation",value=state.get(stage),status="PASS"))
            elif managed and cleaned_source(artifact):
                audit.append(dict(stage=stage,source_file=str(registry_path),source_sha256=digest(registry_path),row=1,
                    sample="ALL",metric="checkpoint_validation",value=CLEANED,status="PASS"))
            else: warnings.append(dict(stage=stage,category="MISSING_EVIDENCE",severity="MISSING_EVIDENCE",
                message="QC summaries lack a matching validated artifact checkpoint",recommendation="Supply original QC validation evidence; preserve existing results",source_file=str(state_path)))
    record_state("qc", qc_state/"state.json")
    record_state("analysis", analysis_state/"state.json")
    for step in ("pca","harmony","cluster-review","regression-review"):
        manifest=analysis_state / f"{step}.json"
        stage="clustering" if step=="cluster-review" else "pca" if step=="regression-review" else step
        if not manifest.is_file():
            continue
        saved=json.loads(manifest.read_text(encoding="utf-8")); inputs[str(manifest)]=digest(manifest)
        outputs=saved.get("outputs",{})
        consistent=bool(outputs) and manifest_valid(saved)
        for name,checksum in outputs.items():
            path=Path(name)
            permitted=allowed(path)
            if not permitted or not path.is_file():
                consistent=False; continue
            actual=inputs.get(str(path)) or digest(path); inputs[str(path)]=actual
            if actual!=checksum: consistent=False
        source=Path(saved.get("binding",{}).get("input",""))
        source_missing=not source.is_file() and not cleaned_source(source)
        if source.is_file():
            actual=inputs.get(str(source)) or digest(source);inputs[str(source)]=actual
            if actual!=saved.get("binding",{}).get("input_sha256"): consistent=False
        status="FAILED" if not consistent else "MISSING_EVIDENCE" if source_missing else (
            "WARNING" if saved.get("status")=="SKIPPED_DESIGN" else "PASS")
        audit.append(dict(stage=stage,source_file=str(manifest),source_sha256=inputs[str(manifest)],
            row=1,sample="ALL",metric=f"{step}_computation_validation",value=saved.get("status"),status=status))
        if status in ("FAILED","MISSING_EVIDENCE"):
            warnings.append(dict(stage=stage,category="FILE_STATE_INCONSISTENCY" if status=="FAILED" else "MISSING_EVIDENCE",
                severity=status,message=f"{step}: published outputs/input and checkpoint disagree or source is missing",
                recommendation="Preserve artifacts; inspect only this stage's provenance and validation log",source_file=str(manifest)))
    for issue in warnings:
        issue["status"] = issue_status(issue["severity"] + " " + issue["category"])
        if resolved and issue["stage"]=="clustering" and issue["status"]=="ACTION_REQUIRED" and (
                "NEEDS_RESOLUTION_SELECTION" in issue["message"] or "resolution" in issue["message"].lower()):
            issue["status"] = "RESOLVED"
        if (resolved and decision.get("regression_scheme")=="baseline" and decision.get("regression_statement","").strip()
                and issue["status"]=="ACTION_REQUIRED" and issue["stage"]=="pca" and
                (issue["category"]=="REGRESSION_BIOLOGICAL_TRADEOFF" or "regression-review: USER_REVIEW_REQUIRED" in issue["message"])):
            issue["status"]="RESOLVED"
    if resolved:
        audit.append(dict(stage="clustering",source_file=str(final_manifest),source_sha256=inputs[str(final_manifest)],
            row=1,sample="ALL",metric="resolution_selection",value=str(decision["resolution"]),status="RESOLVED"))
    return audit, warnings, inputs


def run(gse, root=None, output_root=None, run_id=None):
    root = storage_root(root, output_root)
    runtime_output(root, gse)
    rows, warnings, inputs = collect(root, gse, run_id)
    area = safe_child(analysis_path(root,gse,run_id), "99_workflow_audit")
    area.mkdir(parents=True, exist_ok=True)
    semantic_inputs = dict(inputs)
    semantic_inputs["report_renderer_sha256"] = digest(REPO / "workflow_audit/render_audit.R")
    semantic_inputs["collector_sha256"] = digest(Path(__file__))
    for path in inputs:
        if Path(path).name in ("state.json", "analysis_state.json"):
            state = json.loads(Path(path).read_text(encoding="utf-8"))
            state.pop("updated_at", None)
            state.pop("workflow_audit", None)
            semantic_inputs[path] = hashlib.sha256(json.dumps(state, sort_keys=True).encode()).hexdigest()
    identity = hashlib.sha256(json.dumps(semantic_inputs, sort_keys=True).encode()).hexdigest()
    expected = ("workflow_audit.csv", "workflow_warnings.csv", "workflow_audit.pdf", "audit.log", "audit_summary.json")
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
    current = [w for w in warnings if w["status"] not in ("RESOLVED", "PASS")]
    grouped = []
    for stage, category, status in sorted({(w['stage'],w['category'],w['status']) for w in current}):
        items = [w for w in current if (w['stage'],w['category'],w['status'])==(stage,category,status)]
        grouped.append(dict(stage=stage,category=category,status=status,count=len(items),
            examples=[w['message'] for w in items[:3]],recommendation=items[0]['recommendation']))
    (area / "audit_summary.json").write_text(json.dumps(dict(stage="workflow_audit",status="FAILED" if any(w['status']=='FAILED' for w in current) else
        "ACTION_REQUIRED" if any(w['status']=='ACTION_REQUIRED' for w in current) else "WARNING" if current else "PASS",
        inputs=list(inputs),outputs=[str(area / p) for p in expected],unresolved_issues=grouped,
        needs_user_decision=any(w['status']=='ACTION_REQUIRED' for w in current),resolved_count=sum(w['status']=='RESOLVED' for w in warnings),
        next_command="Resume only the reported unresolved step",limitations="Observed summaries alone do not validate absent computation"),indent=2),encoding="utf-8")
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
    parser.add_argument("--root", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--run-id")
    args = parser.parse_args()
    print(run(args.gse, args.root, args.output_root, args.run_id))
