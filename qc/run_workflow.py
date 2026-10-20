"""Independent/resumable QC-1, QC-2 and QC-3 through the fixed R45 runtime."""
import argparse
import csv
import json
import re
import sys
import os
import uuid
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "runtime"))
from project_paths import storage_root, project_path, runtime_output, dataset_name, flat_storage
from storage import Storage, session, formal_validate, RETAINED, atomic, read, manifest_valid
sys.path.insert(0, str(REPO / "runtime/r45"))
from runtime import RRunError, run_r45
from workflow import checkpoint, digest, failure
from run_r45 import atomic_json

STEPS = {"qc1": ("01_basic_QC", "seurat_QC-1.rds", "qc1"),
         "qc2": ("02_doublet_cellcycle", "seurat_QC-2.rds", "qc"),
         "qc3": ("03_ambient_RNA", "seurat_QC-3.rds", "ambient")}


def report_status(path):
    if not path.is_file():
        return None
    with path.open(encoding="utf-8-sig", newline="") as source:
        rows = list(csv.DictReader(source))
    statuses = [row["value"] for row in rows if row["section"] == "RUN" and row["metric"] == "status"]
    return statuses[-1] if statuses else None


def stage_summary(area, stage, status, source, output, error=None):
    report = area / "qc_report.csv"
    rows = []
    if report.is_file():
        with report.open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
    issues = [row for row in rows if row.get("status") not in ("PASS", "USED", "COMPLETE_QC1", "COMPLETE_QC", "COMPLETE_QC3")]
    needs = status == "NEEDS_USER_DECISION" or any(row.get("status")=="NEEDS_USER_DECISION" for row in issues)
    atomic_json(area / f"{stage}_result.json", dict(stage=stage,status=status,input=str(source) if source else None,
        output=str(output) if output else None,parameters=[r for r in rows if r.get("section")=="CONFIG"],
        statistics=[r for r in rows if r.get("section")=="RUN"],unresolved_issues=issues,error=error,
        needs_user_decision=needs,next_command="Review saved decision gate before resuming" if needs else
        "qc/run_workflow.py DATASET --step " + {"qc1":"qc2", "qc2":"qc3", "qc3":"qc3"}.get(stage,"qc1")))


def resolve_paths(root, gse):
    dataset = project_path(root, gse)
    canonical = dataset / "01_download/seurat_raw.rds"
    legacy = root / "data" / gse
    return dataset, canonical if canonical.is_file() else legacy / "seurat_raw.rds", legacy


def _run(gse, root=None, precheck_only=False, step="all", input_path=None, output_root=None, lifecycle=None,run_id=None):
    dataset_name(gse)
    if step not in (*STEPS, "all"):
        raise ValueError("Expected qc1, qc2, qc3 or all")
    if input_path is not None and step == "all":
        raise ValueError("--input requires one explicit --step")
    root = storage_root(root, output_root)
    runtime_output(root, gse)
    dataset, raw, legacy = resolve_paths(root, gse)
    if lifecycle and lifecycle.active and not raw.is_file():
        raw = lifecycle.referenced("raw") or raw
    explicit = Path(input_path).resolve() if input_path else None
    if explicit and not explicit.is_file():
        raise ValueError(f"Input RDS missing: {explicit}")
    if not raw.is_file() and not explicit:
        raise ValueError("未找到 seurat_raw.rds，请先完成数据下载和原始对象构建。")
    area = dataset / "02_QC"
    scope = "qc"
    if run_id:
        dataset_name(run_id)
        area = area / "rebuilds" / run_id
        scope = "qc/rebuilds/"+run_id
        os.environ["SINGLE_CELL_QC_CHECKPOINT_ROOT"] = lifecycle.checkpoint_root(scope).as_posix()
        os.environ["SINGLE_CELL_QC_RUN_ROOT"] = area.as_posix()
    area.mkdir(parents=True, exist_ok=True)
    if area.resolve() != area:
        raise ValueError("QC directory redirects outside its location")
    state_dir = lifecycle.checkpoint_root(scope+"/state") if lifecycle and lifecycle.active else area / ".workflow"
    def artifact(current):
        directory, filename, _ = STEPS[current]
        if lifecycle and lifecycle.active and current!="qc3":
            return lifecycle.checkpoint_root(scope) / filename
        return area / directory / filename
    protected = {raw: digest(raw)} if raw.is_file() else {}
    if explicit:
        protected[explicit] = digest(explicit)
    stage, script = "healthcheck", None
    current_final = None
    try:
        # Check dependencies at the applicable step so a missing optional DecontX
        # dependency cannot prevent QC-1/QC-2 from being saved.
        with (area / "execution.log").open("a", encoding="utf-8") as log:
            if lifecycle and lifecycle.active:
                if step=="qc1" and explicit: lifecycle.reference("raw",explicit)
                elif raw.is_file(): lifecycle.reference("raw",raw)
                if step in ("all","qc3") and not explicit and artifact("qc3").is_file() and lifecycle.proof(scope):
                    proof = lifecycle.proof(scope)
                    if not lifecycle.valid_proof(proof): raise ValueError("QC final changed after verification")
                    for manifest in state_dir.glob("qc*.json"):
                        if not manifest_valid(read(manifest)):
                            raise ValueError("QC checkpoint missing or changed without valid cleanup evidence")
                    run_r45(REPO / "qc/validate_qc.R",[root,gse,"qc3","--retained"],
                        stage="qc3_validation",log_path=log,cwd=REPO)
                    lifecycle.clean(scope,proof)
                    return report_status(area / STEPS["qc3"][0] / "qc_report.csv")
            run_r45(stage="healthcheck", route="qc1", log_path=log, cwd=REPO)
            checkpoint(state_dir, r_health="PASS", raw_sha256=protected.get(raw))
            if precheck_only:
                stage, script = "qc_precheck", REPO / "qc/qc_precheck.R"
                args = [root, gse, "qc1" if step == "qc1" else "qc2"]
                if explicit:
                    args.append(explicit)
                run_r45(script, args, stage=stage, log_path=log, cwd=REPO)
                status = report_status(area / "qc_report.csv")
                if status != "PASS":
                    raise ValueError(f"QC precheck returned {status}")
                checkpoint(state_dir, qc_precheck="COMPLETE")
                stage_summary(area,"precheck",status,raw,None)
                return status
            selected = list(STEPS) if step == "all" else [step]
            if lifecycle and lifecycle.active and not explicit and step in ("qc2","qc3"):
                selected = [step]
                if step=="qc3" and not artifact("qc2").is_file(): selected.insert(0,"qc2")
                if (step=="qc2" or "qc2" in selected) and not artifact("qc1").is_file(): selected.insert(0,"qc1")
            if step == "qc2" and explicit is None:
                prior = artifact("qc1")
                legacy_final = legacy / "qc/seurat_qc.rds"
                if not prior.is_file() and not legacy_final.is_file():
                    if "qc1" not in selected: selected.insert(0, "qc1")
                elif not prior.is_file() and legacy_final.is_file():
                    stage, script = "qc_legacy_validation", REPO / "qc/validate_qc.R"
                    run_r45(script, [root, gse, "legacy"], stage=stage, log_path=log, cwd=REPO)
                    status = report_status(legacy / "qc/qc_report.csv")
                    checkpoint(state_dir, qc2=status, qc=status, error=None,
                               legacy_qc_sha256=digest(legacy_final))
                    return status
            for current in selected:
                directory, filename, route = STEPS[current]
                destination = area / directory
                destination.mkdir(exist_ok=True)
                current_final = artifact(current)
                stage = current
                if current_final.is_file():
                    stage, script = f"{current}_validation", REPO / "qc/validate_qc.R"
                else:
                    script = REPO / "qc/run_qc.R"
                arguments = [root, gse, current]
                if explicit:
                    arguments.append(explicit)
                with (destination / f"{current}.log").open("a", encoding="utf-8") as step_log:
                    run_r45(stage=f"{current}_healthcheck", route=route, log_path=step_log, cwd=REPO)
                    run_r45(script, arguments, stage=stage, log_path=step_log, cwd=REPO)
                status = report_status(destination / "qc_report.csv")
                allowed = {"qc1": ("COMPLETE_QC1",), "qc2": ("COMPLETE_QC", "QC_COMPLETE_WITH_UNEVALUATED_DOUBLETS"), "qc3": ("COMPLETE_QC3",)}[current]
                if status not in allowed or not current_final.is_file():
                    raise ValueError(f"{current} completion report/artifact missing or incomplete: {status}")
                checkpoint(state_dir, **{current: status, f"{current}_sha256": digest(current_final)},
                           qc=status, error=None)
                source = explicit or (raw if current=="qc1" else
                    artifact("qc1" if current=="qc2" else "qc2"))
                if current=="qc3" and not source.is_file(): source=legacy / "qc/seurat_qc.rds"
                stage_summary(destination,current,status,source,current_final)
                if lifecycle and lifecycle.active:
                    for protected_path,checksum in protected.items():
                        if digest(protected_path)!=checksum:
                            raise RuntimeError(f"Protected input changed before storage commit: {protected_path}")
                    record_path = state_dir / f"{current}.json"
                    atomic(record_path,dict(lifecycle=RETAINED,status=status,binding=dict(input=str(source),
                        input_sha256=digest(source)),outputs={str(current_final):digest(current_final)}))
                    lifecycle.own(current_final,current,scope,formal=current=="qc3",manifest=record_path)
                    if current=="qc3":
                        # Existing QC validator reads the on-disk object in a
                        # separate process before the signature/API validator.
                        run_r45(REPO / "qc/validate_qc.R",[root,gse,current,source],
                            stage="qc3_readback",log_path=log,cwd=REPO)
                        proof = formal_validate(lifecycle,current_final,"qc",scope,run_r45,REPO,log)
                        lifecycle.clean(scope,proof)
            return status
    except BaseException as error:
        if lifecycle and lifecycle.active:
            lifecycle.incomplete(lifecycle.checkpoint_root(scope),stage,scope)
        category = "QC_ERROR"
        needs_decision = report_status(area / "qc_report.csv") == "NEEDS_USER_DECISION"
        if needs_decision and isinstance(error, RRunError) and error.category == "R_SCRIPT_ERROR":
            error.category = category = "QC_NEEDS_USER_DECISION"
        checkpoint(state_dir, **{stage: category}, error=failure(error, stage, script, category),
                   qc="QC_ARTIFACT_CREATED_BUT_RUNTIME_FAILED" if current_final and current_final.is_file() else category)
        stage_summary(area,"precheck" if needs_decision else stage,"NEEDS_USER_DECISION" if needs_decision else "FAILED",
            explicit or raw,current_final,str(error))
        raise
    finally:
        if lifecycle and lifecycle.active: lifecycle.report()
        for path, before in protected.items():
            after = digest(path)
            if path == raw:
                checkpoint(state_dir, raw_sha256_after=after)
            if after != before:
                raise RuntimeError(f"QC_ERROR: protected input changed during QC: {path}")


def run(gse, root=None, precheck_only=False, step="all", input_path=None, output_root=None,
        storage_mode=None, retain_checkpoints=None,run_id=None):
    root = storage_root(root,output_root)
    lifecycle = Storage(project_path(root,gse),storage_mode,retain_checkpoints)
    if run_id and not lifecycle.active: raise ValueError("QC --run-id requires a new managed dataset/output root")
    with session(lifecycle):
        if lifecycle.active and not run_id and step in ("qc1","qc2") and not precheck_only:
            previous = lifecycle.checkpoint_root("qc/state") / f"{step}.json"
            if previous.exists() and (read(previous).get("lifecycle")=="COMPLETE_CLEANED" or
                    not Path(read(previous)["binding"]["input"]).is_file()):
                run_id="rebuild-"+uuid.uuid4().hex[:12]
                if step=="qc2" and input_path is None:
                    upstream=lifecycle.checkpoint_root("qc")/"seurat_QC-1.rds"
                    upstream_manifest=lifecycle.checkpoint_root("qc/state/qc1.json")
                    if upstream.is_file() and upstream_manifest.exists():
                        saved=read(upstream_manifest)
                        if str(upstream) in saved.get("outputs",{}) and manifest_valid(saved): input_path=upstream
                resume=[sys.executable,str(REPO/"qc/run_workflow.py"),gse,
                    "--output-root" if flat_storage(root) else "--root",str(root),"--run-id",run_id,"--step",step]
                if input_path: resume += ["--input",str(input_path)]
                atomic(lifecycle.dataset / "02_QC/last_rebuild_result.json",dict(run_id=run_id,
                    requested_step=step,output=str(lifecycle.dataset / "02_QC/rebuilds" / run_id),
                    resume_command=resume))
        if lifecycle.active and step in ("qc1","qc2","qc3"):
            role=("qc/rebuilds/"+run_id if run_id else "qc")+":"+step+"_input"
            if input_path: lifecycle.reference(role,input_path)
            elif lifecycle.registry["references"].get(role): input_path=lifecycle.referenced(role)
        return _run(gse,root,precheck_only,step,input_path,lifecycle=lifecycle,run_id=run_id)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gse")
    parser.add_argument("--root", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--run-id",help="Isolated QC checkpoint reconstruction; preserves canonical QC-3")
    parser.add_argument("--storage-mode",choices=("compact","full"))
    parser.add_argument("--retain-checkpoint",dest="retain_checkpoints",action="append",
                        help="Keep stage name or RDS basename (repeatable)")
    parser.add_argument("--precheck", dest="precheck_only",action="store_true")
    parser.add_argument("--step", choices=[*STEPS, "all"], default="all")
    parser.add_argument("--input", dest="input_path", type=Path)
    print(run(**vars(parser.parse_args(argv))))


if __name__ == "__main__": main()
