"""Actual Windows R runtime acceptance tests; no package changes or downloads."""
import argparse
import json
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "runtime/r45"))
from config import RUNTIME_LOG
from runtime import RRunError, run_r45
from workflow import digest


def main(real_fixture=False):
    results = []
    health = run_r45(stage="runtime_test_health", route="qc", cwd=ROOT)
    print(health.stdout, flush=True)
    expressions = {"arithmetic": 'stopifnot(1+1==2); cat("PASS\\n")'}
    for package in ("cli", "ggplot2", "SeuratObject", "Seurat", "DoubletFinder"):
        expressions[package] = f'library({package}); stopifnot(Sys.getenv("PROCESSOR_ARCHITECTURE") != "ARM64"); cat("PASS\\n")'
    for name, expression in expressions.items():
        result = run_r45(expression=expression, stage="runtime_test_" + name, cwd=ROOT)
        assert result.returncode == 0
        assert "C.UTF-8 failed" not in result.stderr
        results.append(dict(test=name, status="PASS", exit_code=result.returncode))
        print(f"PASS {name}: exit 0", flush=True)
    try:
        run_r45(expression='library(cli); stop("RUNTIME_TEST_ERROR")', stage="runtime_test_error", cwd=ROOT)
    except RRunError as error:
        assert error.exit_code == 1 and "RUNTIME_TEST_ERROR" in error.output
        results.append(dict(test="intentional_error", status="PASS", exit_code=error.exit_code))
        print("PASS intentional error: exit 1", flush=True)
    else:
        raise AssertionError("R error was swallowed")
    with tempfile.TemporaryDirectory(prefix="r45-roundtrip-") as temporary:
        result = run_r45(Path(__file__).with_name("roundtrip.R"), [temporary], stage="runtime_test_roundtrip", cwd=ROOT)
        results.append(dict(test="RDS_UTF8_roundtrip", status="PASS", exit_code=result.returncode))
        print(result.stdout, flush=True)
    if real_fixture:
        artifact = ROOT / "data/GSE225857/seurat_raw.rds"
        before = digest(artifact)
        result = run_r45(ROOT / ".agents/skills/geo-single-cell-loader/scripts/validate_seurat.R",
                         ["--existing", ROOT, "GSE225857"], stage="read_only_regression", cwd=ROOT)
        assert before == digest(artifact)
        results.append(dict(test="GSE225857", status="PASS", exit_code=result.returncode, sha256=before))
        print(result.stdout, flush=True)
    else:
        results.append(dict(test="GSE225857", status="NOT RUN"))
    (RUNTIME_LOG / "acceptance.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--real-fixture", action="store_true")
    main(parser.parse_args().real_fixture)
