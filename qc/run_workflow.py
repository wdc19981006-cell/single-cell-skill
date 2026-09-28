"""Run independent QC with runtime health, checkpoints and existing-artifact recovery."""
import argparse
import csv
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "runtime/r45"))
from runtime import RRunError, run_r45
from workflow import checkpoint, digest, failure


def report_status(path):
    if not path.is_file():
        return None
    with path.open(encoding="utf-8-sig", newline="") as source:
        rows = list(csv.DictReader(source))
    statuses = [r["value"] for r in rows if r["section"] == "RUN" and r["metric"] == "status"]
    return statuses[-1] if statuses else None


def run(gse, root=REPO, precheck_only=False):
    if not re.fullmatch(r"GSE\d+", gse):
        raise ValueError("Expected GSE accession")
    root = Path(root).resolve()
    dataset = root / "data" / gse
    if dataset.resolve() != dataset:
        raise ValueError("Dataset directory redirects outside its location")
    raw = dataset / "seurat_raw.rds"
    if not raw.is_file():
        raise ValueError("未找到 seurat_raw.rds，请先完成数据下载和原始对象构建。")
    area = dataset / "qc"
    area.mkdir(exist_ok=True)
    state_dir = dataset / ".workflow"
    final = area / "seurat_qc.rds"
    before = digest(raw)
    stage = "healthcheck"
    script = None
    try:
        with (area / "execution.log").open("a", encoding="utf-8") as log:
            run_r45(stage=stage, route="qc", log_path=log, cwd=REPO)
            checkpoint(state_dir, r_health="PASS", raw_sha256=before)
            if final.is_file():
                stage = "qc_validation"
                script = REPO / "qc/validate_qc.R"
                run_r45(script, [root, gse], stage=stage, log_path=log, cwd=REPO)
            else:
                stage = "qc_precheck"
                script = REPO / "qc/qc_precheck.R"
                run_r45(script, [root, gse], stage=stage, log_path=log, cwd=REPO)
                checkpoint(state_dir, qc_precheck="COMPLETE")
                if precheck_only:
                    return "PASS"
                stage = "qc"
                script = REPO / "qc/run_qc.R"
                run_r45(script, [root, gse], stage=stage, log_path=log, cwd=REPO)
            status = report_status(area / "qc_report.csv")
            if status not in ("COMPLETE_QC", "QC_COMPLETE_WITH_UNEVALUATED_DOUBLETS"):
                raise ValueError("QC completion report missing or incomplete")
            checkpoint(state_dir, qc=status, error=None, qc_sha256=digest(final))
            return status
    except Exception as error:
        category = "QC_ERROR"
        if (stage == "qc_precheck" and report_status(area / "qc_report.csv") == "NEEDS_USER_DECISION"
                and isinstance(error, RRunError) and error.category == "R_SCRIPT_ERROR"):
            error.category = category = "QC_NEEDS_USER_DECISION"
        checkpoint(state_dir, **{stage: category}, error=failure(error, stage, script, category),
                   qc="QC_ARTIFACT_CREATED_BUT_RUNTIME_FAILED" if final.is_file() else category)
        raise
    finally:
        after = digest(raw)
        checkpoint(state_dir, raw_sha256_after=after)
        if after != before:
            raise RuntimeError("QC_ERROR: seurat_raw.rds changed during QC")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gse")
    parser.add_argument("--root", type=Path, default=REPO)
    parser.add_argument("--precheck", action="store_true")
    args = parser.parse_args()
    print(run(args.gse, args.root, args.precheck))
