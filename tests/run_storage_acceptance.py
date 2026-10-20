"""Real R45/Seurat/Harmony acceptance on synthetic data; never run on a real GSE."""
import argparse
import importlib.util
import json
import sys
from pathlib import Path

REPO=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO/"runtime/r45"))
from runtime import run_r45
sys.path.insert(0,str(REPO/"runtime"))
from storage import sha256, read, atomic, Storage

spec=importlib.util.spec_from_file_location("storage_acceptance_analysis",REPO/"analysis/run_workflow.py")
analysis=importlib.util.module_from_spec(spec); spec.loader.exec_module(analysis)


def run(output,resume=False):
    output=Path(output).resolve()
    if output.exists() and not resume: raise ValueError("Use a new acceptance output directory or --resume")
    output.mkdir(parents=True,exist_ok=True)
    if not resume:
        with (output/"acceptance.log").open("w",encoding="utf-8") as log:
            run_r45(stage="healthcheck",route="analysis_pca",log_path=log,cwd=REPO)
            run_r45(REPO/"tests/test_seurat_storage.R",[REPO,output/"layers"],stage="storage_layers",log_path=log,cwd=REPO)
            run_r45(REPO/"tests/test_analysis.R",[REPO,output],stage="storage_fixture",log_path=log,cwd=REPO)
    source=output/"synthetic_qc.rds"; before=sha256(source)
    config=output/"analysis_config.json"; cfg=read(config)
    cfg.update(risk_correlation=0,sensitivity=["count"])
    if not resume: atomic(config,cfg)
    dataset=output/"results/SYNTHETIC_STORAGE"; area=dataset/"03_analysis"
    assert analysis.run("SYNTHETIC_STORAGE",output,input_path=source,config=config)=="NEEDS_RESOLUTION_SELECTION"
    work=dataset/".work/checkpoints/analysis/baseline"
    assert (work/"seurat_harmony.rds").exists()
    assert (work/"count/seurat_sensitivity.rds").exists()
    review=work/"seurat_cluster_review.rds"; reviewed=sha256(review)
    assert analysis.run("SYNTHETIC_STORAGE",output,config=config)=="NEEDS_RESOLUTION_SELECTION"
    assert sha256(review)==reviewed
    choice=output/"synthetic_choice.json"
    atomic(choice,dict(confirmed_by="user",user_statement="SYNTHETIC ACCEPTANCE ONLY: choose .6",resolution=.6,
        regression_scheme="baseline",regression_statement="SYNTHETIC ACCEPTANCE ONLY: retain baseline"))
    assert analysis.run("SYNTHETIC_STORAGE",output,step="finalize",config=config,resolution_decision=choice)=="COMPLETE"
    final=area/"seurat_analysis.rds"; final_hash=sha256(final)
    assert not review.exists(); assert not (work/"count/seurat_sensitivity.rds").exists()
    assert analysis.run("SYNTHETIC_STORAGE",output,config=config)=="COMPLETE"
    assert sha256(final)==final_hash and sha256(source)==before
    audit_spec=importlib.util.spec_from_file_location("storage_acceptance_audit",REPO/"workflow_audit/run_audit.py")
    audit=importlib.util.module_from_spec(audit_spec); audit_spec.loader.exec_module(audit)
    audit_rows,warnings,_=audit.collect(output,"SYNTHETIC_STORAGE")
    assert not [r for r in warnings if r["category"]=="FILE_STATE_INCONSISTENCY"]
    assert any(r["metric"]=="resolution_selection" and r["status"]=="RESOLVED" for r in audit_rows)
    with (output/"acceptance.log").open("a",encoding="utf-8") as log:
        run_r45(REPO/"runtime/validate_formal.R",[final,"analysis"],stage="storage_acceptance_independent",log_path=log,cwd=REPO)
    result=dict(status="PASS",fixture="SYNTHETIC_ONLY",r_version="4.5.0",source_sha256=before,
        final_sha256=final_hash,storage=read(dataset/"storage_summary.json"),
        verified=["exact RNA and fractional decontX counts", "normalized layers", "PCA/Harmony/unique UMAP",
          "pause/resume at resolution review", "literal synthetic selection", "independent final RDS readback",
          "FeaturePlot/FindMarkers/annotation API", "regression candidate cleanup", "manifest cleanup tombstones",
          "restart after cleanup", "existing final unchanged", "audit accepts cleaned dependencies"])
    atomic(output/"test_results.json",result)
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--resume",action="store_true")
    run(**vars(parser.parse_args()))
