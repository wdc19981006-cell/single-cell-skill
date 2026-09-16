import gzip
import importlib.util
import csv
import io
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "examples" / "gse149614" / "prepare.py"
SPEC = importlib.util.spec_from_file_location("gse149614_prepare", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class GSE149614MappingTests(unittest.TestCase):
    def test_metadata_requires_unique_cells_and_matching_patient(self):
        by_title = {"HCC01T": {"gsm": "GSM1", "patient": "HCC01", "site": "Tumor"}}
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "metadata.tsv.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("Cell\tsample\tsite\tpatient\n")
                handle.write("HCC01T_AAAC\tHCC01T\tTumor\tHCC01\n")
                handle.write("HCC01T_AAAC\tHCC01T\tTumor\tHCC01\n")
            with self.assertRaisesRegex(ValueError, "duplicate"):
                MODULE.read_metadata(path, by_title, expected_cells=2)
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("Cell\tsample\tsite\tpatient\n")
                handle.write("HCC01T_AAAC\tHCC01T\tTumor\tHCC02\n")
            with self.assertRaisesRegex(ValueError, "Patient/site mismatch"):
                MODULE.read_metadata(path, by_title, expected_cells=1)

    def test_no_group_is_created_from_sample_site(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "group_confirmation.json"
            path.write_text(json.dumps({"confirmed_by": "user", "user_statement": "HCC01T 为 A", "groups": {"GSM1": "A"}}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "every GSM"):
                MODULE.check_confirmation(path, {"GSM1", "GSM2"})
            path.unlink()
            with self.assertRaisesRegex(ValueError, "group 尚未创建"):
                MODULE.check_confirmation(path, {"GSM1", "GSM2"})

    def test_two_samples_share_one_explicit_pooled_input(self):
        fields = ["sample", "metadata_evidence"]
        rows = [{"sample": "GSM1", "metadata_evidence": "GEO"},
                {"sample": "GSM2", "metadata_evidence": "GEO"}]
        text = MODULE.prepare_report(rows, fields)
        result = list(csv.DictReader(io.StringIO(text)))
        self.assertEqual(result[0]["local_path"], result[1]["local_path"])
        self.assertEqual(result[0]["cell_map_path"], result[1]["cell_map_path"])
        self.assertEqual(result[0]["files_json"], result[1]["files_json"])
        self.assertEqual(result[0]["feature_column"], "__row_names__")
        self.assertEqual(result[0]["count_source"], "counts")


if __name__ == "__main__":
    unittest.main()
