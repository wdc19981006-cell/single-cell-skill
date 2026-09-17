"""Run the confirmed Stage B pipeline and publish an audit in all terminal states."""
import argparse
import json
import os
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from audit_run import finalize, restore_provenance, run_id
from common import ROOT, dataset_paths

SCRIPTS = Path(__file__).resolve().parent


def execute(label, command, log, root):
    started = time.perf_counter()
    log.write(f"\n[{datetime.now(timezone.utc).isoformat()}] {label}: {command}\n")
    log.flush()
    result = subprocess.run(command, cwd=root, stdout=log, stderr=subprocess.STDOUT,
                            env={**os.environ, "GEO_AUDIT_RUN_ACTIVE": "1"})
    log.write(f"[{datetime.now(timezone.utc).isoformat()}] {label} exit={result.returncode}\n")
    log.flush()
    if result.returncode:
        raise RuntimeError(f"{label} failed with exit code {result.returncode}")
    return time.perf_counter() - started


def run(gse, root=ROOT, rscript="Rscript"):
    paths = dataset_paths(root, gse)
    workflow = paths["workflow"]
    if not (workflow / "audit_run.json").exists():
        raise ValueError("Start this real GSE run with audit_run.py start before Stage A; existing workflow files require review")
    run_id(workflow)
    manifest = f"data/{gse}/.workflow/sample_manifest.csv"
    status, reason = "success", ""
    run_started = time.perf_counter()
    try:
        with (workflow / "execution.log").open("a", encoding="utf-8") as log:
            log.write(f"Run started UTC: {datetime.now(timezone.utc).isoformat()}\n")
            restore_provenance(gse, root)
            execute("download", [sys.executable, str(SCRIPTS / "download_processed.py"), manifest, "--root", str(root)], log, root)
            execute("prebuild_probe", [sys.executable, str(SCRIPTS / "prebuild_probe.py"), manifest, "--root", str(root)], log, root)
            execute("dependencies", [rscript, str(SCRIPTS / "check_dependencies.R")], log, root)
            execute("build", [rscript, str(SCRIPTS / "build_seurat.R"), str(root), manifest], log, root)
            independent_validation_seconds = execute("validation", [rscript, str(SCRIPTS / "validate_seurat.R"), str(root), manifest, f"data/{gse}/seurat_raw.rds"], log, root)
        build_profile = workflow / "build_profile.json"
        if build_profile.exists():
            profile = json.loads(build_profile.read_text(encoding="utf-8"))
            profile["independent_validation_seconds"] = independent_validation_seconds
            profile["pipeline_total_seconds"] = time.perf_counter() - run_started
            probe = workflow / "candidate_probe.json"
            if probe.exists():
                profile["candidate_probe"] = json.loads(probe.read_text(encoding="utf-8"))
            build_profile.write_text(json.dumps(profile, indent=2), encoding="utf-8")
        (workflow / "validation.json").write_text(json.dumps({
            "status": "success", "validator": "validate_seurat.R",
            "validated_at_utc": datetime.now(timezone.utc).isoformat(),
            "rds_bytes": paths["final"].stat().st_size,
        }, indent=2), encoding="utf-8")
    except KeyboardInterrupt:
        status, reason = "interrupted", "User or process interruption"
    except Exception as error:
        status, reason = "failure", str(error)
    if status != "success":
        (workflow / "validation.json").write_text(json.dumps({"status": status, "reason": reason}, indent=2), encoding="utf-8")
        info = paths["info"]
        if info.exists() and info.read_text(encoding="utf-8").startswith("STATUS: COMPLETE"):
            old = info.read_text(encoding="utf-8")
            info.write_text(old.replace("STATUS: COMPLETE", "STATUS: BUILD_FAILED", 1) + f"\nFailure: {reason}\n", encoding="utf-8")
    try:
        url = finalize(gse, root, status, reason)
    except Exception as error:
        raise RuntimeError(f"Audit upload failed; local workflow retained for retry: {error}") from error
    print(f"Audit: {url}")
    if status != "success":
        raise RuntimeError(f"{status}: {reason}")
    return url


def main():
    def interrupted(_signum, _frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, interrupted)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gse")
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--rscript", default=os.environ.get("GEO_RSCRIPT", "Rscript"))
    args = parser.parse_args()
    run(args.gse, args.root.resolve(), args.rscript)


if __name__ == "__main__":
    main()
