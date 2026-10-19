"""Read-only QC-3 diagnosis in a new versioned directory, without rerunning DecontX."""
import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(REPO/"runtime"),str(REPO/"runtime/r45")]
from project_paths import storage_root, project_path, runtime_output, dataset_name
from runtime import run_r45
from workflow import digest
from run_r45 import atomic_json


def run(dataset, input_path, output_root=None, review_id=None, config=None):
    root=storage_root(output_root=output_root); runtime_output(root,dataset)
    source=Path(input_path).resolve()
    before=digest(source)
    review_id=dataset_name(review_id or datetime.now(timezone.utc).strftime("review-%Y%m%dT%H%M%S%fZ"))
    area=project_path(root,dataset)/"02_QC/03_ambient_RNA/reviews"/review_id
    if area.resolve()!=area: raise ValueError("Review directory redirects outside output location")
    area.mkdir(parents=True,exist_ok=False)
    try:
        with (area/"review.log").open("w",encoding="utf-8") as log:
            args=[source,area]+([Path(config).resolve()] if config else [])
            run_r45(REPO/"qc/review_ambient.R",args,stage="ambient_review",cwd=REPO,log_path=log)
        if digest(source)!=before: raise ValueError("QC-3 source changed during read-only review")
        atomic_json(area/"review_provenance.json",dict(input=str(source),input_sha256=before,
            outputs={p.name:digest(p) for p in area.iterdir() if p.is_file()},status="COMPLETE"))
    except BaseException as error:
        atomic_json(area/"failure.json",dict(status="FAILED",message=str(error)))
        raise
    return area


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset"); parser.add_argument("--input",dest="input_path",type=Path,required=True)
    parser.add_argument("--output-root",type=Path); parser.add_argument("--review-id"); parser.add_argument("--config",type=Path)
    print(run(**vars(parser.parse_args())))
