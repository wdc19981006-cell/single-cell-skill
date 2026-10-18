"""Actual synthetic R 4.5 test of analysis stages, finalization, audit and no-op resume.

The selected resolution is a test fixture, never a choice for a real dataset.
Use a new output directory so no prior artifact or decision can be replaced.
"""
import argparse
import importlib.util
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "runtime/r45"))
from runtime import run_r45
from workflow import digest

spec = importlib.util.spec_from_file_location("synthetic_analysis_workflow", REPO / "analysis/run_workflow.py")
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


def run(output):
    output = Path(output).resolve()
    if output.exists():
        raise ValueError("Use a new output directory; existing synthetic test artifacts are preserved")
    output.mkdir(parents=True)
    with (output / "smoke.log").open("w", encoding="utf-8") as log:
        run_r45(stage="healthcheck", route="analysis_pca", log_path=log, cwd=REPO)
        run_r45(REPO / "tests/test_analysis.R", [REPO, output],
                stage="analysis_synthetic_fixture", log_path=log, cwd=REPO)
    source = output / "synthetic_qc.rds"
    config = output / "analysis_config.json"
    original = digest(source)
    gse = "GSE999999"
    area = output / "results" / gse / "03_analysis"
    assert analysis.run(gse, output, input_path=source, config=config) == "NEEDS_RESOLUTION_SELECTION"
    assert not (area / "seurat_analysis.rds").exists(), "Formal labels appeared before selection"
    assert (area / "seurat_harmony.rds").is_file(), "Crossed technical batch did not exercise Harmony"
    decision = output / "synthetic_resolution_decision.json"
    decision.write_text(json.dumps({"confirmed_by": "user", "resolution": .6,
        "user_statement": "TEST FIXTURE ONLY: simulate a user selecting resolution 0.6 on synthetic cells"}),encoding="utf-8")
    assert analysis.run(gse, output, step="finalize", config=config, resolution_decision=decision) == "COMPLETE"
    final = area / "seurat_analysis.rds"
    final_hash = digest(final)
    assert analysis.run(gse, output, input_path=source, config=config) == "COMPLETE"
    assert digest(final) == final_hash
    assert digest(source) == original
    report = {"dataset": "SYNTHETIC_GSE999999_LABEL_ONLY", "cells": 160, "genes": 400,
              "r_version": "4.5.0", "status": "PASS", "formal_selection": "TEST FIXTURE ONLY",
              "input_sha256": original, "final_sha256": final_hash,
              "verified": ["fresh PCA", "crossed-batch Harmony", "same-coordinate resolution scan",
                           "clustree rendering", "user decision gate", "final RDS roundtrip",
                           "original counts preserved", "audit", "resume without recomputation"]}
    (output / "test_results.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args().output)
