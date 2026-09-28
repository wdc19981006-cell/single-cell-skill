"""Run the confirmed Stage B pipeline and publish an audit in all terminal states."""
import argparse
import json
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from audit_run import finalize, restore_provenance, restore_workflow, run_id
from common import ROOT, dataset_paths, read_csv
from stage_a_policy import UnsupportedModality, assess_single_cell_modality

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "runtime/r45"))
from runtime import RRunError, run_r45
from workflow import checkpoint, failure


def error_record(error, stage):
    category = getattr(error, "category", None) or {
        "download": "CHECKSUM_ERROR" if "checksum" in str(error).lower() or "provenance" in str(error).lower() else "DOWNLOAD_ERROR",
        "prebuild_probe": "FORMAT_ERROR", "validation": "RAW_VALIDATION_ERROR",
        "build": "RAW_BUILD_ERROR",
    }.get(stage, "R_SCRIPT_ERROR")
    record = failure(error, stage, category=category)
    record["workflow_category"] = {"build": "RAW_BUILD_ERROR", "validation": "RAW_VALIDATION_ERROR"}.get(stage, category)
    return record


def raw_handoff(gse, paths):
    """Report the validated raw stage without initiating the independent QC stage."""
    info = paths["info"].read_text(encoding="utf-8").splitlines()
    if not info or info[0] != "STATUS: COMPLETE_RAW":
        raise RuntimeError("Raw handoff requires STATUS: COMPLETE_RAW")
    def value(label):
        index = info.index(label)
        return info[index + 1]
    rows = read_csv(paths["workflow"] / "sample_summary.csv")
    if not rows:
        raise RuntimeError("Raw handoff requires a sample summary")
    groups = {}
    for row in rows:
        groups[row["group"]] = groups.get(row["group"], 0) + int(row["cells"])
    return "\n".join((
        "STATUS: COMPLETE_RAW",
        f"GSE: {gse}",
        f"Samples: {value('Total samples:')}",
        f"Cells: {value('Total cells:')}",
        f"Genes: {value('Total genes:')}",
        "Groups (cells): " + "; ".join(f"{group}={cells}" for group, cells in groups.items()),
        f"Raw expression: data/{gse}/raw/",
        f"Raw Seurat: data/{gse}/seurat_raw.rds",
        "原始 Seurat 数据已经构建并验证完成。是否继续进行 QC（细胞质量过滤、DoubletFinder、细胞周期评分）？",
    ))


def execute(label, command, log, root):
    started = time.perf_counter()
    log.write(f"\n[{datetime.now(timezone.utc).isoformat()}] {label}: {command}\n")
    log.flush()
    import os
    result = subprocess.run(command, cwd=root, capture_output=True, text=True, encoding="utf-8", errors="replace",
                            env={**os.environ, "GEO_AUDIT_RUN_ACTIVE": "1", "PYTHONIOENCODING": "utf-8"})
    output = result.stdout + result.stderr
    log.write(output)
    log.write(f"[{datetime.now(timezone.utc).isoformat()}] {label} exit={result.returncode}\n")
    log.flush()
    if result.returncode:
        error = RuntimeError(f"{label} failed with exit code {result.returncode}: {output[-2000:]}")
        error.script, error.exit_code = command[1], result.returncode
        if label == "prebuild_probe" and any(word in output.lower() for word in ("mapping", "cell map", "cell_map")):
            error.category = "MAPPING_ERROR"
        raise error
    return time.perf_counter() - started


def execute_r(label, script, args, log, root):
    started = time.perf_counter()
    log.write(f"\n[{datetime.now(timezone.utc).isoformat()}] {label}: {script} {args}\n")
    log.flush()
    run_r45(script, args, stage=label, log_path=log, cwd=root, route="base" if script is None else None)
    log.write(f"[{datetime.now(timezone.utc).isoformat()}] {label} exit=0\n")
    log.flush()
    return time.perf_counter() - started


def run(gse, root=ROOT):
    paths = dataset_paths(root, gse)
    workflow = paths["workflow"]
    if not workflow.exists() and (paths["dataset"] / "group_confirmation.json").is_file():
        restore_workflow(gse, root)
    if not (workflow / "audit_run.json").exists():
        raise ValueError("Start this real GSE run with audit_run.py start before Stage A; existing workflow files require review")
    run_id(workflow)
    manifest = f"data/{gse}/.workflow/sample_manifest.csv"
    status, reason = "success", ""
    handoff = None
    independent_validation_seconds = None
    run_started = time.perf_counter()
    stage = "discovery"
    try:
        with (workflow / "execution.log").open("a", encoding="utf-8") as log:
            log.write(f"Run started UTC: {datetime.now(timezone.utc).isoformat()}\n")
            inspection_path = workflow / "inspection.json"
            if not inspection_path.exists():
                raise UnsupportedModality("Stage A inspection is missing; single-cell modality is not established")
            modality = assess_single_cell_modality(json.loads(inspection_path.read_text(encoding="utf-8")))
            (workflow / "modality_gate.json").write_text(json.dumps(modality, ensure_ascii=False, indent=2), encoding="utf-8")
            if modality["status"] != "single_cell":
                raise UnsupportedModality(modality["reason"] + "; STOP before processed expression download and Seurat construction")
            checkpoint(workflow, discovery="COMPLETE")
            if not paths["final"].is_file():
                stage = "download"
                restore_provenance(gse, root)
                execute("download", [sys.executable, str(SCRIPTS / "download_processed.py"), manifest, "--root", str(root)], log, root)
                checkpoint(workflow, download="DOWNLOAD_COMPLETE", group="CONFIRMED")
                stage = "prebuild_probe"
                execute("prebuild_probe", [sys.executable, str(SCRIPTS / "prebuild_probe.py"), manifest, "--root", str(root)], log, root)
                checkpoint(workflow, prebuild="COMPLETE")
            else:
                checkpoint(workflow, raw_build="ARTIFACT_PRESENT")
                log.write("Existing seurat_raw.rds: skip download, prebuild and build; independently validate.\n")
            stage = "healthcheck"
            execute_r("healthcheck", None, [], log, root)
            checkpoint(workflow, r_health="PASS")
            if not paths["final"].is_file():
                stage = "dependencies"
                execute_r("dependencies", SCRIPTS / "check_dependencies.R", [], log, root)
                stage = "build"
                execute_r("build", SCRIPTS / "build_seurat.R", [str(root), manifest], log, root)
                checkpoint(workflow, raw_build="COMPLETE")
            stage = "validation"
            independent_validation_seconds = execute_r("validation", SCRIPTS / "validate_seurat.R",
                                                       [str(root), manifest, f"data/{gse}/seurat_raw.rds", "--record-completion"], log, root)
            checkpoint(workflow, raw_validation="COMPLETE", raw_build="COMPLETE_RAW", group="CONFIRMED", error=None)
    except KeyboardInterrupt:
        status, reason = "interrupted", "User or process interruption"
        checkpoint(workflow, error=error_record(KeyboardInterrupt(reason), stage))
    except Exception as error:
        status, reason = ("interrupted" if getattr(error, "category", None) == "INTERRUPTED" else "failure"), str(error)
        if paths["final"].is_file():
            artifact_state = ("RAW_ARTIFACT_CREATED_BUT_RUNTIME_FAILED"
                              if getattr(error, "category", None) == "R_NATIVE_CRASH" or stage == "build"
                              else "RAW_ARTIFACT_PRESENT_BUT_VALIDATION_FAILED" if stage == "validation"
                              else "ARTIFACT_PRESENT")
        else:
            artifact_state = "FAILED" if stage == "build" else "NOT_COMPLETE"
        checkpoint(workflow, error=error_record(error, stage),
                   raw_build=artifact_state)
    if status != "success":
        (workflow / "validation.json").write_text(json.dumps({"status": status, "reason": reason,
            "error": json.loads((workflow / "state.json").read_text(encoding="utf-8")).get("error")}, indent=2), encoding="utf-8")
        info = paths["info"]
        if info.exists() and info.read_text(encoding="utf-8").startswith("STATUS: COMPLETE_RAW") and not paths["final"].is_file():
            old = info.read_text(encoding="utf-8")
            info.write_text(old.replace("STATUS: COMPLETE_RAW", "STATUS: BUILD_FAILED", 1) + f"\nFailure: {reason}\n", encoding="utf-8")
    if status == "success":
        (workflow / "validation.json").write_text(json.dumps({
            "status": "success", "validator": "validate_seurat.R",
            "validated_at_utc": datetime.now(timezone.utc).isoformat(),
            "rds_bytes": paths["final"].stat().st_size,
        }, indent=2), encoding="utf-8")
        try:
            build_profile = workflow / "build_profile.json"
            if build_profile.exists():
                profile = json.loads(build_profile.read_text(encoding="utf-8"))
                profile["independent_validation_seconds"] = independent_validation_seconds
                profile["pipeline_total_seconds"] = time.perf_counter() - run_started
                probe = workflow / "candidate_probe.json"
                if probe.exists():
                    profile["candidate_probe"] = json.loads(probe.read_text(encoding="utf-8"))
                build_profile.write_text(json.dumps(profile, indent=2), encoding="utf-8")
            handoff = raw_handoff(gse, paths)
        except Exception as error:
            handoff = f"HANDOFF_WARNING: {error}; validated raw object remains COMPLETE_RAW"
            with (workflow / "execution.log").open("a", encoding="utf-8") as log:
                log.write(handoff + "\n")
    try:
        url = finalize(gse, root, status, reason)
    except Exception as error:
        raise RuntimeError(f"Audit upload failed; local workflow retained for retry: {error}") from error
    print(f"Audit: {url}")
    if status != "success":
        raise RuntimeError(f"{status}: {reason}")
    print(handoff)
    return url


def main():
    def interrupted(_signum, _frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, interrupted)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gse")
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    run(args.gse, args.root.resolve())


if __name__ == "__main__":
    main()
