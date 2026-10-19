"""Actual small-fixture acceptance in the configured output root, never a real GSE analysis."""
import argparse
import importlib.util
import json
import sys
from pathlib import Path

REPO=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("upgrade_acceptance",REPO/"analysis/run_workflow.py")
analysis=importlib.util.module_from_spec(spec);spec.loader.exec_module(analysis)
sys.path.insert(0,str(REPO/"runtime"))
from project_paths import storage_root, project_path, runtime_output


def run(fixture,dataset,output_root=None):
    root=storage_root(output_root=output_root);runtime_output(root,dataset)
    project=project_path(root,dataset)
    evidence=project/"validation";evidence.mkdir(parents=True,exist_ok=True)
    source=Path(fixture).resolve();before=analysis.digest(source)
    cfg=evidence/"default_sample_config.json"
    cfg.write_text(json.dumps(dict(nfeatures=200,max_pcs=12,markers=["g1","g50","g100"],
        proliferation_markers=["g2"],resolutions=[.2,.3,.6],regression_review=dict(comparison_resolution=.3))),encoding="utf-8")
    assert analysis.run(dataset,input_path=source,config=cfg,output_root=root,run_id="default-sample-actual") == "NEEDS_RESOLUTION_SELECTION"
    area=analysis.analysis_path(root,dataset,"default-sample-actual")
    assert (area/"seurat_harmony.rds").is_file()
    decision=evidence/"fixture_resolution_03.json"
    decision.write_text(json.dumps(dict(confirmed_by="user",resolution=.3,
        user_statement="TEST FIXTURE ONLY: simulate selecting 0.3; not a real dataset decision",
        regression_scheme="baseline",regression_statement="TEST FIXTURE ONLY: retain pMT baseline after reviewing synthetic count comparison")),encoding="utf-8")
    assert analysis.run(dataset,step="finalize",config=cfg,output_root=root,run_id="default-sample-actual",resolution_decision=decision)=="COMPLETE"
    final=area/"seurat_analysis.rds";final_sha=analysis.digest(final)
    assert analysis.run(dataset,input_path=source,config=cfg,output_root=root,run_id="default-sample-actual")=="COMPLETE"
    assert analysis.digest(final)==final_sha and analysis.digest(source)==before
    selected=area/"03_clustering/selected_resolution_0.3"
    assert all((selected/name).is_file() for name in ["marker_dotplot.pdf","marker_featureplots.pdf","qc_cycle_distributions.pdf","sample_group_composition.pdf"])
    # A targeted cycle comparison requires an explicit fixture-only marker interpretation.
    cycle=evidence/"cycle_fixture_config.json"
    cycle.write_text(json.dumps(dict(sensitivity=["cell_cycle"],risk_correlation=0,
        markers=["g1","g50"],proliferation_markers=["g2"],regression_review=dict(auto_numeric=False,comparison_resolution=.3,
            cell_cycle_review=dict(confirmed_by="user",user_statement="TEST FIXTURE ONLY: compare cycle regression",
                marker_interpretation="Synthetic test genes only; no biological interpretation")))),encoding="utf-8")
    analysis.run(dataset,step="regression-review",input_path=area/"seurat_cluster_review.rds",config=cycle,
        output_root=root,run_id="cycle-fixture-comparison")
    candidate=analysis.analysis_path(root,dataset,"cycle-fixture-comparison")/"01_pca_regression/regression_sensitivity/cell_cycle"
    assert (candidate/"seurat_sensitivity.rds").is_file()
    # Adoption must start a new formal run and recompute PC choice, with real evidence binding.
    adoption=evidence/"fixture_regression_choice.json"
    adoption.write_text(json.dumps(dict(scheme="cell_cycle",confirmed_by="user",user_statement="TEST FIXTURE ONLY: adopt comparison in isolated PCA",
        comparison_config=str(candidate/"comparison_config.json"))),encoding="utf-8")
    assert analysis.run(dataset,step="pca",input_path=source,config=cfg,output_root=root,run_id="adoption-fixture-pca",
        regression_decision=adoption)=="COMPLETE"
    # Original RNA is an explicit independent alternative, never a silent fallback.
    rna_cfg=evidence/"rna_config.json";rna_cfg.write_text(json.dumps(dict(assay="RNA",nfeatures=200,max_pcs=12)),encoding="utf-8")
    assert analysis.run(dataset,step="pca",input_path=source,config=rna_cfg,output_root=root,run_id="rna-fixture-pca")=="COMPLETE"
    assert analysis.digest(final)==final_sha and analysis.digest(source)==before
    result=dict(status="PASS",data="SYNTHETIC fixture with actual installed celda correction",dataset=dataset,
        verified=["flat default output","default sample Harmony despite confounding","fractional decontX input",
            "resolution 0.3 marker plots","resolved historical gate","no-op resume","full cycle comparison with fixture review",
            "isolated approved regression adoption","explicit RNA PCA alternative","immutable source/final"],
        input_sha256=before,final_sha256=final_sha)
    (evidence/"acceptance_results.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--fixture",type=Path,required=True)
    p.add_argument("--dataset",default="TEST_upgrade_20261019");p.add_argument("--output-root",type=Path)
    run(**vars(p.parse_args()))
