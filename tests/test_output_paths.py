"""Flat output roots, safe names, isolation and explicit legacy input compatibility."""
import importlib.util
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO / "runtime"))
from project_paths import storage_root, project_path, analysis_path, dataset_name, DEFAULT_OUTPUT_ROOT
sys.path.insert(0,str(REPO / ".agents/skills/geo-single-cell-loader/scripts"))
from common import dataset_paths, validate_manifest

spec=importlib.util.spec_from_file_location("path_analysis",REPO / "analysis/run_workflow.py")
analysis=importlib.util.module_from_spec(spec); spec.loader.exec_module(analysis)


class OutputPathsTests(unittest.TestCase):
    def test_default_is_central_desktop_config(self):
        with patch.dict(os.environ,{},clear=True):
            self.assertEqual(storage_root(),DEFAULT_OUTPUT_ROOT.resolve())
            self.assertEqual(project_path(storage_root(),"GSE149614"),DEFAULT_OUTPUT_ROOT.resolve()/"GSE149614")
            self.assertEqual(storage_root(REPO),REPO)

    def test_independent_stages_and_datasets(self):
        with tempfile.TemporaryDirectory() as tmp, patch.dict(os.environ):
            root=storage_root(REPO,tmp)
            download=dataset_paths(root,"GSE123")
            self.assertEqual(download["dataset"],root/"GSE123/01_download")
            self.assertEqual(analysis_path(root,"GSE123","decontx-v1"),root/"GSE123/03_analysis/runs/decontx-v1")
            self.assertNotEqual(project_path(root,"Study_A"),project_path(root,"Study_B"))
            row=dict(database="GSE123",sample="S1",tissue="x",disease="x",source_type="scRNA",group="test",
                local_path="GSE123/01_download/raw/S1",file_type="10x_mtx",count_source="raw",count_evidence="fixture",metadata_evidence="fixture")
            validate_manifest([row],root)
            with self.assertRaises(ValueError): validate_manifest([dict(row,local_path="GSE124/01_download/raw/S1")],root)

    def test_non_geo_and_windows_name_safety(self):
        for name in ("patient_A", "Study-2026.1", "GSE149614"):
            self.assertEqual(dataset_name(name),name)
        for name in ("../escape", "a/b", "a\\b", "C:/out", "a.", "CON", "NUL.txt", "COM1", "", ".hidden"):
            with self.subTest(name=name),self.assertRaises(ValueError): dataset_name(name)

    def test_old_explicit_workspace_does_not_flatten_or_move(self):
        with tempfile.TemporaryDirectory() as tmp, patch.dict(os.environ):
            root=Path(tmp).resolve()
            old=root/"results/GSE123/02_QC/03_ambient_RNA/seurat_QC-3.rds"
            old.parent.mkdir(parents=True); old.write_bytes(b"immutable historical fixture")
            target=storage_root(REPO,root/"new-output")
            self.assertEqual(analysis.qc_input(target,"GSE123",old),old)
            self.assertEqual(analysis.qc_input(root,"GSE123"),old)
            self.assertEqual(old.read_bytes(),b"immutable historical fixture")
            self.assertFalse(project_path(target,"GSE123").exists())

    def test_analysis_run_does_not_discover_prior_run_artifacts(self):
        with tempfile.TemporaryDirectory() as tmp, patch.dict(os.environ):
            root=storage_root(REPO,tmp)
            old=analysis_path(root,"Study_A")
            old.mkdir(parents=True); (old/"seurat_harmony.rds").write_bytes(b"old")
            with self.assertRaisesRegex(ValueError,"prior analysis RDS"):
                analysis.run("Study_A",output_root=root,run_id="new",step="cluster-review")
            self.assertEqual((old/"seurat_harmony.rds").read_bytes(),b"old")


if __name__=="__main__": unittest.main()
