"""Real existing-cohort regression in an isolated root; no GEO discovery/download.

The original raw RDS and literal prior QC decision are copied without edits.
The real-cohort cluster review deliberately awaits a user's final selection.
"""
import argparse
import json
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "runtime/r45"))
from runtime import run_r45
from workflow import digest


def main(gse, raw_source, decision_source, root, start="qc1"):
    if not re.fullmatch(r"GSE\d+", gse):
        raise ValueError("Expected GSE accession")
    root = Path(root).resolve()
    if not root.is_relative_to(ROOT / ".runtime"):
        raise ValueError("Acceptance outputs must stay inside this repository's ignored .runtime")
    raw_source, decision_source = Path(raw_source).resolve(), Path(decision_source).resolve()
    protected = {str(path): digest(path) for path in (raw_source, decision_source)}
    project = root / "results" / gse
    raw = project / "01_download/seurat_raw.rds"
    decision = project / "02_QC/qc_decision.json"
    for source, target in ((raw_source, raw), (decision_source, decision)):
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            if digest(target) != digest(source):
                raise ValueError(f"Existing fixture differs; preserved: {target}")
        else:
            shutil.copy2(source, target)
    record = dict(gse=gse, scope="complete existing raw cohort; no discovery or download",
        root=str(root), raw_source=str(raw_source), raw_sha256=digest(raw),
        original_decision=str(decision_source), protected_before=protected,
        started_at=datetime.now(timezone.utc).isoformat(), stages=[], status="RUNNING")
    path = root / "test_results.json"
    try:
        steps = ["qc1", "qc2", "qc3", "analysis", "validation", "resume"]
        with (root / "acceptance.log").open("a", encoding="utf-8") as log:
            for step in steps[steps.index(start):]:
                if step == "validation":
                    result = run_r45(ROOT / "tests/validate_staged_real.R", [ROOT, root, gse],
                        stage="staged_real_validation", cwd=ROOT, log_path=log)
                    record["stages"].append(dict(step=step, exit_code=result.returncode, status="PASS"))
                    continue
                if step == "analysis":
                    command = [sys.executable, str(ROOT / "analysis/run_workflow.py"), gse, "--root", str(root)]
                elif step == "resume":
                    command = [sys.executable, str(ROOT / "qc/run_workflow.py"), gse, "--root", str(root), "--step", "qc3"]
                else:
                    command = [sys.executable, str(ROOT / "qc/run_workflow.py"), gse, "--root", str(root), "--step", step]
                log.write("\nCOMMAND " + json.dumps(command) + "\n")
                log.flush()
                result = subprocess.run(command, cwd=ROOT, stdout=log, stderr=log)
                record["stages"].append(dict(step=step, exit_code=result.returncode,
                    status="PASS" if result.returncode == 0 else "FAILED"))
                path.write_text(json.dumps(record, indent=2), encoding="utf-8")
                if result.returncode:
                    raise RuntimeError(f"{step} failed with exit {result.returncode}; inspect acceptance.log")
                print("PASS " + step, flush=True)
            record["protected_after"] = {p: digest(Path(p)) for p in protected}
            if protected != record["protected_after"] or digest(raw) != record["raw_sha256"]:
                raise ValueError("A protected source changed during acceptance")
            record["status"] = "PASS_TO_CLUSTER_REVIEW"
            record["resolution_status"] = "NEEDS_USER_SELECTION; no real-cohort seurat_analysis.rds created"
    except BaseException as error:
        record["status"] = "INTERRUPTED" if isinstance(error, KeyboardInterrupt) else "FAILED"
        record["error"] = str(error)
        raise
    finally:
        record["ended_at"] = datetime.now(timezone.utc).isoformat()
        record["artifacts"] = {str(p.relative_to(root)): dict(bytes=p.stat().st_size, sha256=digest(p))
            for p in project.rglob("*") if p.is_file() and ".workflow" not in p.parts}
        path.write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8")
    print(path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gse")
    parser.add_argument("--raw-source", required=True, type=Path)
    parser.add_argument("--decision-source", required=True, type=Path)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--start", choices=("qc1", "qc2", "qc3", "analysis", "validation", "resume"), default="qc1")
    options = parser.parse_args()
    main(**vars(options))
