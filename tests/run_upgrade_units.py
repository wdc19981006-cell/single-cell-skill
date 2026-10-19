"""Run only regression suites affected by the workflow upgrade."""
import argparse
import json
import os
import unittest
from pathlib import Path
from unittest.mock import patch

REPO=Path(__file__).resolve().parents[1]
PATTERNS=("test_loader.py","test_download_layout.py","test_qc_workflow.py","test_analysis_workflow.py",
          "test_workflow_audit.py","test_output_paths.py","test_runtime.py")


def run(output):
    suites=[unittest.defaultTestLoader.discover(str(REPO/"tests"),pattern=pattern) for pattern in PATTERNS]
    with patch.dict(os.environ):
        result=unittest.TextTestRunner(verbosity=1).run(unittest.TestSuite(suites))
    record=dict(status="PASS" if result.wasSuccessful() else "FAILED",tests=result.testsRun,
        failures=len(result.failures),errors=len(result.errors),skipped=len(result.skipped),suites=list(PATTERNS))
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(record,indent=2),encoding="utf-8")
    print(json.dumps(record));return 0 if result.wasSuccessful() else 1


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--output",type=Path,required=True)
    raise SystemExit(run(parser.parse_args().output))
