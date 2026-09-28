"""Run existing R loader/QC regressions in an isolated synthetic workspace."""
import json
import os
import shutil
import sys
import uuid
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "runtime/r45"))
from config import RUNTIME_LOG
from runtime import run_r45
from workflow import digest
import make_fixtures
import run_end_to_end


def main():
    fixture = RUNTIME_LOG / ("regression-" + uuid.uuid4().hex[:8])
    scripts = ".agents/skills/geo-single-cell-loader/scripts"
    shutil.copytree(ROOT / scripts, fixture / scripts, ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copytree(ROOT / "qc", fixture / "qc", ignore=shutil.ignore_patterns("__pycache__"))
    make_fixtures.ROOT = fixture
    make_fixtures.main()
    results = []
    env = dict(os.environ, GEO_SINGLE_CELL_PYTHON=sys.executable, PYTHONIOENCODING="utf-8",
               QC_TEST_RETAIN_FIXTURE="1")
    run_r45(route="qc", stage="regression_health", cwd=ROOT, env=env)
    run_end_to_end.main(fixture)
    results.append(dict(test="synthetic_E2E", status="PASS"))
    dataset = fixture / "data/GSE999999999"
    final = dataset / "seurat_raw.rds"
    before = digest(final)
    info = dataset / "sample_info.txt"
    original = info.read_text(encoding="utf-8")
    info.write_text(original.replace("STATUS: COMPLETE_RAW", "STATUS: BUILD_FAILED", 1), encoding="utf-8")
    run_r45(ROOT / scripts / "validate_seurat.R",
            [fixture, "data/GSE999999999/.workflow/sample_manifest.csv",
             "data/GSE999999999/seurat_raw.rds", "--record-completion"],
            stage="synthetic_raw_artifact_recovery", cwd=fixture, env=env,
            log_path=RUNTIME_LOG / "regression-output.log")
    assert info.read_text(encoding="utf-8").startswith("STATUS: COMPLETE_RAW")
    assert digest(final) == before
    results.append(dict(test="existing_raw_artifact_recovery", status="PASS", sha256=before))
    print("PASS existing raw artifact recovered by validation only", flush=True)
    for name in ("test_10x_literal_na.R", "test_final_object.R", "test_global_min_cells.R",
                 "test_pooled_inputs.R", "test_seurat.R", "test_stream_text.R",
                 "test_utf8_manifest.R", "test_qc.R", "test_qc_real_df.R"):
        args = [fixture, sys.executable] if name == "test_stream_text.R" else [fixture]
        result = run_r45(ROOT / "tests" / name, args, stage=name, cwd=fixture, env=env,
                         log_path=RUNTIME_LOG / "regression-output.log")
        results.append(dict(test=name, status="PASS", exit_code=result.returncode,
                            duration_seconds=result.record.get("duration_seconds")))
        (RUNTIME_LOG / "regressions.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
        print("PASS " + name, flush=True)
    qc_fixture = Path((fixture / "data/qc_fixture_path.txt").read_text().strip())
    spec = importlib.util.spec_from_file_location("qc_workflow", ROOT / "qc/run_workflow.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.run("GSE999999997", qc_fixture) == "COMPLETE_QC"
    results.append(dict(test="QC_existing_artifact_recovery", status="PASS"))
    (RUNTIME_LOG / "regressions.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print("PASS QC existing-artifact recovery", flush=True)
    print("Synthetic fixture retained for review: " + str(fixture), flush=True)


if __name__ == "__main__":
    main()
