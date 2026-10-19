"""Independent formal PCA, guarded Harmony and user-reviewed clustering in R 4.5.0."""
import argparse
import importlib.util
import json
import os
import re
import sys
import uuid
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "runtime"))
from project_paths import storage_root, project_path, analysis_path, runtime_output, dataset_name, flat_storage
sys.path.insert(0, str(REPO / "runtime/r45"))
from runtime import run_r45
from run_r45 import atomic_json
from workflow import checkpoint, digest, failure

STEPS = {
    "pca": ("01_pca_regression", "seurat_pca.rds", "analysis_pca"),
    "harmony": ("02_harmony", "seurat_harmony.rds", "harmony"),
    "cluster-review": ("03_clustering", "seurat_cluster_review.rds", "clustering"),
    "regression-review": ("01_pca_regression/regression_sensitivity", "regression_review.rds", "analysis_pca"),
    "finalize": ("03_clustering", "seurat_analysis.rds", "analysis_pca"),
}


def configuration(path=None):
    cfg = json.loads((REPO / "analysis/default_config.json").read_text(encoding="utf-8"))
    overrides = json.loads(Path(path).read_text(encoding="utf-8-sig")) if path else {}
    if not isinstance(overrides, dict) or set(overrides) - set(cfg):
        raise ValueError("Unknown analysis configuration keys")
    for key, value in overrides.items():
        if isinstance(cfg[key], dict):
            if not isinstance(value, dict) or set(value) - set(cfg[key]):
                raise ValueError(f"Unknown {key} configuration keys")
            cfg[key].update(value)
        else:
            cfg[key] = value
    return cfg


def qc_input(root, gse, explicit=None):
    if explicit:
        path = Path(explicit).resolve()
        if not path.is_file():
            raise ValueError(f"Analysis input does not exist: {path}")
        return path
    project = project_path(root, gse)
    candidates = [project / "02_QC/03_ambient_RNA/seurat_QC-3.rds",
                  project / "02_QC/02_doublet_cellcycle/seurat_QC-2.rds",
                  root / "data" / gse / "qc/seurat_qc.rds"]
    for path in candidates:
        if path.is_file():
            return path.resolve()
    raise ValueError("No completed QC-3/QC-2 object; provide --input for an existing QC RDS")


def selection(path):
    if not path:
        raise ValueError("Finalization requires --resolution-decision written from the user's choice")
    decision = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    value = decision.get("resolution")
    if (decision.get("confirmed_by") != "user" or not isinstance(decision.get("user_statement"), str)
            or not decision["user_statement"].strip() or isinstance(value, bool)
            or not isinstance(value, (float, int)) or value <= 0):
        raise ValueError("Resolution decision requires confirmed_by=user, user_statement and positive resolution")
    return decision


def stage_configuration(cfg, step):
    keys = {
        "pca": ("assay", "seed", "nfeatures", "max_pcs", "pc_cumulative_variance",
                "pc_individual_variance", "pc_elbow_drop", "risk_correlation", "regression_decision"),
        "harmony": ("seed", "harmony", "markers"),
        "cluster-review": ("seed", "resolutions", "cluster_reduction"),
        "regression-review": tuple(cfg),
        "finalize": ("markers", "proliferation_markers"),
    }[step]
    return {key: cfg[key] for key in keys}


def validate(artifact, kind, log):
    run_r45(REPO / "analysis/validate_analysis.R", [artifact, kind],
            stage="analysis_validation", log_path=log, cwd=REPO)


def audit(gse, root, run_id=None):
    spec = importlib.util.spec_from_file_location("analysis_workflow_audit", REPO / "workflow_audit/run_audit.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.run(gse, root, run_id=run_id) if run_id else module.run(gse, root)


def publish(pending_path):
    """Finish an interrupted publication without recomputing or replacing any result."""
    pending = json.loads(pending_path.read_text(encoding="utf-8"))
    for item in pending["files"]:
        src, dst = Path(item["source"]), Path(item["target"])
        if dst.exists():
            if not dst.is_file() or digest(dst) != item["sha256"]:
                raise ValueError(f"Existing analysis output differs; preserved: {dst}")
        else:
            if not src.is_file() or digest(src) != item["sha256"]:
                raise ValueError(f"Interrupted output is unavailable: {src}")
            dst.parent.mkdir(parents=True, exist_ok=True)
            # Staging is on the same volume. A hard link atomically creates an
            # exclusive destination, so an interruption cannot leave a partial
            # final artifact and a concurrent process cannot be overwritten.
            try:
                os.link(src, dst)
            except FileExistsError:
                if not dst.is_file() or digest(dst) != item["sha256"]:
                    raise ValueError(f"Concurrent analysis output differs; preserved: {dst}")
            if digest(dst) != item["sha256"]:
                raise ValueError(f"Published output checksum mismatch: {dst}")
    manifest = {k: v for k, v in pending.items() if k not in ("files", "manifest")}
    manifest["outputs"] = {item["target"]: item["sha256"] for item in pending["files"]}
    atomic_json(Path(pending["manifest"]), manifest)
    pending_path.unlink()


def execute(step, source, area, cfg, decision=None, continuation_command=None):
    directory, artifact_name, route = STEPS[step]
    output = area / directory
    output.mkdir(parents=True, exist_ok=True)
    state = area / ".workflow"
    state.mkdir(exist_ok=True)
    binding = {"step": step, "input": str(source), "input_sha256": digest(source),
               "configuration": stage_configuration(cfg, step), "decision": decision}
    manifest = state / f"{step}.json"
    pending_path = state / f"{step}.pending.json"
    final = (output if step=="regression-review" else area) / artifact_name
    log_name = {"pca": "pca.log", "harmony": "harmony.log", "regression-review": "regression_review.log",
                "cluster-review": "clustering.log", "finalize": "finalize.log"}[step]
    with (output / log_name).open("a", encoding="utf-8") as log:
        if pending_path.exists():
            pending = json.loads(pending_path.read_text(encoding="utf-8"))
            if pending["binding"] != binding:
                raise ValueError("An interrupted publication has different inputs/configuration; preserved")
            publish(pending_path)
        if manifest.exists():
            record = json.loads(manifest.read_text(encoding="utf-8"))
            if record["binding"] != binding:
                raise ValueError(f"Completed {step} has different inputs/configuration; use a separate project root")
            for path, checksum in record["outputs"].items():
                if not Path(path).is_file() or digest(path) != checksum:
                    raise ValueError(f"Completed output is missing or changed: {path}")
            for path in record["outputs"]:
                if Path(path).suffix == ".rds":
                    kind = "pca" if Path(path).name == "seurat_unintegrated.rds" else step
                    if "regression_sensitivity" not in path:
                        validate(Path(path), kind, log)
                    elif Path(path).name=="seurat_sensitivity.rds":
                        validate(Path(path), "sensitivity", log)
            checkpoint(state, **{step: record["status"]}, error=None)
            return record["status"]
        if final.exists():
            validate(final, step, log)
            raise ValueError(f"Existing {artifact_name} has no matching manifest; validated and preserved")
        # A completed computation can survive a validation/publication failure.
        # Recover it only with the launcher's successful exit record: an RDS or
        # PASS in stdout alone never upgrades a native crash to success.
        staging = None
        completion_name = f"{step.replace('-', '_')}_runtime_completion.json"
        for candidate in sorted((state / "staging").glob(f"{step}-*")):
            complete = candidate / completion_name
            request_path = candidate / "request.json"
            if not complete.is_file() or not request_path.is_file():
                continue
            previous = json.loads(request_path.read_text(encoding="utf-8"))
            runtime_record = json.loads(complete.read_text(encoding="utf-8"))
            if (previous["provenance"] == binding and runtime_record.get("status") == "PASS"
                    and runtime_record.get("exit_code") == 0):
                staging = candidate
                break
        recovered = staging is not None
        if not recovered:
            # An optional Harmony package is required only after its design gate.
            health_route = "analysis_pca" if step == "harmony" else route
            run_r45(stage="healthcheck", route=health_route, log_path=log, cwd=REPO)
            staging = state / "staging" / f"{step}-{uuid.uuid4().hex}"
            staging.mkdir(parents=True)
            atomic_json(staging / "request.json", {"configuration": cfg, "decision": decision,
                                                   "provenance": binding})
        request = staging / "request.json"
        stage_name = f"analysis_{step}"
        script = REPO / "analysis/run_analysis.R"
        before = binding["input_sha256"]
        try:
            if not recovered:
                computation = run_r45(script, [step, source, staging, request], stage=stage_name,
                                      log_path=log, cwd=REPO)
                record = getattr(computation, "record", {})
                atomic_json(staging / completion_name, {"status": "PASS", "exit_code": 0,
                    "runtime_run_id": record.get("run_id"), "input_sha256": before})
            result_path = staging / f"{step.replace('-', '_')}_result.json"
            result = json.loads(result_path.read_text(encoding="utf-8"))
            if result["status"] not in ("COMPLETE", "NEEDS_RESOLUTION_SELECTION", "SKIPPED_DESIGN", "USER_REVIEW_REQUIRED"):
                raise ValueError("Analysis stage did not complete")
            artifact = staging / artifact_name
            if step != "harmony" or result["status"] != "SKIPPED_DESIGN":
                if not artifact.is_file():
                    raise ValueError(f"Stage completed without {artifact_name}")
                validate(artifact, step, log)
            if step == "harmony":
                validate(staging / "seurat_unintegrated.rds", "pca", log)
            if digest(source) != before:
                raise ValueError("Analysis input changed during computation")
            def publication_target(file):
                relative = Path(file).relative_to(staging)
                if relative.as_posix()==artifact_name:
                    return final
                if relative.as_posix()=="seurat_unintegrated.rds":
                    return area / relative
                return output / relative
            def published_path(value):
                if isinstance(value, str) and value.startswith(str(staging)):
                    return str(publication_target(value))
                if isinstance(value, list): return [published_path(v) for v in value]
                if isinstance(value, dict): return {k: published_path(v) for k,v in value.items()}
                return value
            result = published_path(result)
            result["input"] = str(source)
            result["outputs"] = [str(publication_target(file)) for file in sorted(staging.rglob("*"))
                                 if file.is_file() and file!=request]
            result["output_object"] = str(area / "seurat_unintegrated.rds" if
                step=="harmony" and result["status"]=="SKIPPED_DESIGN" else final)
            if continuation_command:
                result["next_command"] = continuation_command
            atomic_json(result_path, result)
            items = []
            for file in sorted(staging.rglob("*")):
                if not file.is_file() or file == request:
                    continue
                target = publication_target(file)
                items.append({"source": str(file), "target": str(target), "sha256": digest(file)})
            atomic_json(pending_path, {"binding": binding, "status": result["status"],
                                     "files": items, "manifest": str(manifest)})
            publish(pending_path)
            checkpoint(state, **{step: result["status"]}, error=None)
            return result["status"]
        except BaseException as error:
            checkpoint(state, **{step: "FAILED"}, error=failure(error, stage_name, script))
            atomic_json(output / f"{step.replace('-', '_')}_failure.json", dict(stage=step,status="FAILED",
                input=str(source),output=str(final),parameters=stage_configuration(cfg,step),
                unresolved_issues=[str(error)],needs_user_decision=False,
                next_command="Resume this same fixed stage after resolving the concrete error"))
            raise
        finally:
            if digest(source) != before:
                raise RuntimeError("Analysis source RDS changed; results were not approved for use")


def run(gse, root=None, step="all", input_path=None, config=None, resolution_decision=None,
        output_root=None, run_id=None, regression_decision=None):
    dataset_name(gse)
    root = storage_root(root, output_root)
    runtime_output(root, gse)
    area = analysis_path(root, gse, run_id)
    area.mkdir(parents=True, exist_ok=True)
    cfg = configuration(config)
    if regression_decision:
        cfg["regression_decision"] = json.loads(Path(regression_decision).read_text(encoding="utf-8-sig"))
    if cfg["regression_decision"] is not None:
        d = cfg["regression_decision"]
        if (not run_id or d.get("confirmed_by")!="user" or not d.get("user_statement","").strip()
                or d.get("scheme") not in ("feature","count","cell_cycle")):
            raise ValueError("Adopting regression requires a new --run-id and a literal user scheme decision")
        evidence = json.loads(Path(d["comparison_config"]).read_text(encoding="utf-8-sig"))
        expected = {"feature":["pMT","nFeature_RNA"],"count":["pMT","nCount_RNA"],"cell_cycle":["pMT","S.Score","G2M.Score"]}[d["scheme"]]
        if (evidence.get("input_sha256") != digest(qc_input(root,gse,input_path)) or
                evidence.get("configuration",{}).get("assay") != cfg["assay"] or evidence.get("regress") != expected):
            raise ValueError("Regression decision comparison must use this same QC input and assay")
    steps = ["pca", "harmony", "cluster-review", "regression-review"] if step == "all" else [step]
    if step == "all" and resolution_decision:
        selection(resolution_decision)
        steps.append("finalize")
    status = None
    for current in steps:
        if current == "pca":
            source = qc_input(root, gse, input_path)
        elif input_path and step != "all":
            source = qc_input(root, gse, input_path)
        elif current == "harmony":
            source = area / "seurat_pca.rds"
        elif current == "cluster-review":
            source = (area / "seurat_harmony.rds" if (area / "seurat_harmony.rds").is_file()
                      else area / "seurat_unintegrated.rds")
            if not source.is_file():
                source = area / "seurat_pca.rds"
        else:
            source = area / "seurat_cluster_review.rds"
        if not source.is_file():
            raise ValueError(f"Required prior analysis RDS is missing: {source}")
        decision = selection(resolution_decision) if current == "finalize" else None
        comparison_manifest = area / ".workflow/regression-review.json"
        if current=="finalize" and comparison_manifest.is_file():
            comparison=json.loads(comparison_manifest.read_text(encoding="utf-8"))
            if comparison.get("status")=="USER_REVIEW_REQUIRED" and (
                    decision.get("regression_scheme")!="baseline" or not decision.get("regression_statement","").strip()):
                raise ValueError("Regression tradeoff needs the user's formal-scheme choice: record regression_scheme=baseline and literal regression_statement to keep this run, or adopt a candidate in a new --run-id")
        next_step={"pca":"harmony","harmony":"cluster-review","cluster-review":"regression-review",
                   "regression-review":"finalize","finalize":None}[current]
        command=[sys.executable,str(REPO / "analysis/run_workflow.py"),gse,
                 "--output-root" if flat_storage(root) else "--root",str(root)]
        if run_id: command += ["--run-id",run_id]
        if config: command += ["--config",str(Path(config).resolve())]
        if regression_decision: command += ["--regression-decision",str(Path(regression_decision).resolve())]
        if next_step: command += ["--step",next_step]
        if next_step=="finalize": command += ["--resolution-decision","<user-confirmed-decision.json>"]
        if current=="finalize": command=[sys.executable,str(REPO / "workflow_audit/run_audit.py"),gse,
            "--output-root" if flat_storage(root) else "--root",str(root)]+(["--run-id",run_id] if run_id else [])
        status = execute(current, source.resolve(), area, cfg, decision, command)
        if current == "regression-review" and step == "all":
            status = "NEEDS_RESOLUTION_SELECTION"
    # A later `all` invocation also validates a previously authorized final
    # selection; the stored literal user decision remains authoritative.
    final_manifest = area / ".workflow/finalize.json"
    if step == "all" and not resolution_decision and final_manifest.is_file():
        previous = json.loads(final_manifest.read_text(encoding="utf-8"))
        status = execute("finalize", (area / "seurat_cluster_review.rds").resolve(), area,
                         cfg, previous["binding"]["decision"])
    if step in ("all", "cluster-review", "finalize", "regression-review"):
        try:
            audit(gse, root, run_id) if run_id else audit(gse, root)
        except BaseException as error:
            checkpoint(area / ".workflow", workflow_audit="FAILED",
                       error=failure(error, "workflow_audit", REPO / "workflow_audit/run_audit.py"))
            raise
    return status


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gse")
    parser.add_argument("--root", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--run-id", help="Isolated analysis run; never reads another run's checkpoints")
    parser.add_argument("--regression-decision", type=Path)
    parser.add_argument("--step", choices=[*STEPS, "all"], default="all")
    parser.add_argument("--input", dest="input_path", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--resolution-decision", type=Path)
    options = parser.parse_args()
    print(run(**vars(options)))
