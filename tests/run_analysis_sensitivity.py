"""Actual targeted-regression and explicit corrected-counts PCA branch tests."""
import argparse
import csv
import importlib.util
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("sensitivity_analysis_workflow", REPO / "analysis/run_workflow.py")
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


def rows(path):
    with path.open(encoding="utf-8-sig",newline="") as source:
        return list(csv.DictReader(source))


def run(fixture, output):
    output = Path(output).resolve()
    if output.exists():
        raise ValueError("Use a new output directory")
    output.mkdir(parents=True)
    fixture = Path(fixture).resolve()
    before = analysis.digest(fixture)
    config = output / "sensitivity.json"
    config.write_text(json.dumps({"nfeatures":200,"max_pcs":12,
        "assay":"RNA","sensitivity":["rna_total","cell_cycle"],
        "markers":["g1","g50"],"regression_review":{"auto_numeric":False}}),encoding="utf-8")
    root = output / "rna"
    analysis.run("GSE999999",root,input_path=fixture,config=config)
    area = root / "results/GSE999999/03_analysis/01_pca_regression"
    measured = rows(area / "regression_diagnostics.csv")
    max_rho = max(abs(float(row["rho"])) for row in measured if row["variable"]=="nCount_RNA" and row["rho"])
    assert max_rho >= .5, "Fixture must exhibit a measured RNA-total/PC risk before sensitivity"
    sensitivity = {row["target"]:row for row in rows(area / "regression_sensitivity/sensitivity_summary.csv")}
    assert sensitivity["count"]["status"]=="COMPLETE"
    assert sensitivity["cell_cycle"]["status"]=="SKIPPED_NO_CLEAR_RISK"
    assert (area / "regression_sensitivity/count/seurat_sensitivity.rds").is_file()
    assert not (area / "regression_sensitivity/cell_cycle").exists()
    corrected_config = output / "corrected.json"
    corrected_config.write_text(json.dumps({"nfeatures":200,"max_pcs":12,"assay":"decontX"}),encoding="utf-8")
    corrected_root = output / "corrected"
    analysis.run("GSE999999",corrected_root,step="pca",input_path=fixture,config=corrected_config)
    corrected = corrected_root / "results/GSE999999/03_analysis/01_pca_regression"
    assert {row["metric"]:row["value"] for row in rows(corrected / "pca_summary.csv")}["assay"]=="decontX"
    assert analysis.digest(fixture)==before
    result = {"status":"PASS","data":"synthetic fractional corrected-counts assay",
        "nCount_RNA_max_abs_pc_rho":max_rho,"threshold":.5,
        "rna_total_alias_to_count":"COMPLETE","cell_cycle_sensitivity":"SKIPPED_NO_CLEAR_RISK",
        "corrected_assay":"decontX","original_RNA_counts":"PRESERVED_AND_ROUNDTRIP_VALIDATED"}
    (output / "test_results.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    run(args.fixture,args.output)
