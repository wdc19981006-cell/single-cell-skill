"""Independent/resumable QC-1, QC-2 and QC-3 through the fixed R45 runtime."""
import argparse
import csv
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "runtime/r45"))
from runtime import RRunError, run_r45
from workflow import checkpoint, digest, failure

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


def resolve_paths(root, gse):
    dataset = root / "results" / gse
    canonical = dataset / "01_download/seurat_raw.rds"
    legacy = root / "data" / gse
    return dataset, canonical if canonical.is_file() else legacy / "seurat_raw.rds", legacy


def run(gse, root=REPO, precheck_only=False, step="all", input_path=None):
    if not re.fullmatch(r"GSE\d+", gse):
        raise ValueError("Expected GSE accession")
    if step not in (*STEPS, "all"):
        raise ValueError("Expected qc1, qc2, qc3 or all")
    if input_path is not None and step == "all":
        raise ValueError("--input requires one explicit --step")
    root = Path(root).resolve()
    dataset, raw, legacy = resolve_paths(root, gse)
    explicit = Path(input_path).resolve() if input_path else None
    if explicit and not explicit.is_file():
        raise ValueError(f"Input RDS missing: {explicit}")
    if not raw.is_file() and not (explicit and step in ("qc2", "qc3")):
        raise ValueError("未找到 seurat_raw.rds，请先完成数据下载和原始对象构建。")
    area = dataset / "02_QC"
    area.mkdir(parents=True, exist_ok=True)
    if area.resolve() != area:
        raise ValueError("QC directory redirects outside its location")
    state_dir = area / ".workflow"
    protected = {raw: digest(raw)} if raw.is_file() else {}
    if explicit:
        protected[explicit] = digest(explicit)
    stage, script = "healthcheck", None
    current_final = None
    try:
        # Check dependencies at the applicable step so a missing optional DecontX
        # dependency cannot prevent QC-1/QC-2 from being saved.
        with (area / "execution.log").open("a", encoding="utf-8") as log:
            run_r45(stage="healthcheck", route="qc1", log_path=log, cwd=REPO)
            checkpoint(state_dir, r_health="PASS", raw_sha256=protected.get(raw))
            if precheck_only:
                stage, script = "qc_precheck", REPO / "qc/qc_precheck.R"
                args = [root, gse, "qc1" if step == "qc1" else "qc2"]
                run_r45(script, args, stage=stage, log_path=log, cwd=REPO)
                status = report_status(area / "qc_report.csv")
                if status != "PASS":
                    raise ValueError(f"QC precheck returned {status}")
                checkpoint(state_dir, qc_precheck="COMPLETE")
                return status
            selected = list(STEPS) if step == "all" else [step]
            if step == "qc2" and explicit is None:
                prior = area / STEPS["qc1"][0] / STEPS["qc1"][1]
                legacy_final = legacy / "qc/seurat_qc.rds"
                if not prior.is_file() and not legacy_final.is_file():
                    selected.insert(0, "qc1")
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
                current_final = destination / filename
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
            return status
    except Exception as error:
        category = "QC_ERROR"
        needs_decision = report_status(area / "qc_report.csv") == "NEEDS_USER_DECISION"
        if needs_decision and isinstance(error, RRunError) and error.category == "R_SCRIPT_ERROR":
            error.category = category = "QC_NEEDS_USER_DECISION"
        checkpoint(state_dir, **{stage: category}, error=failure(error, stage, script, category),
                   qc="QC_ARTIFACT_CREATED_BUT_RUNTIME_FAILED" if current_final and current_final.is_file() else category)
        raise
    finally:
        for path, before in protected.items():
            after = digest(path)
            if path == raw:
                checkpoint(state_dir, raw_sha256_after=after)
            if after != before:
                raise RuntimeError(f"QC_ERROR: protected input changed during QC: {path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gse")
    parser.add_argument("--root", type=Path, default=REPO)
    parser.add_argument("--precheck", action="store_true")
    parser.add_argument("--step", choices=[*STEPS, "all"], default="all")
    parser.add_argument("--input", dest="input_path", type=Path)
    args = parser.parse_args()
    print(run(args.gse, args.root, args.precheck, args.step, args.input_path))
