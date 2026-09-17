"""Regression coverage for split gzip raw counts with a large CSV matrix."""
import gzip
import io
import json
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

SCRIPTS = Path(__file__).resolve().parents[1] / ".agents/skills/geo-single-cell-loader/scripts"
sys.path.insert(0, str(SCRIPTS))
from common import validate_manifest, write_csv
from build_manifest import confirm
from download_processed import assemble_parts, download, sha256
from stream_text_to_sparse import convert


class SplitRawCountsTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.raw = self.root / "data/GSE123/raw"
        self.raw.mkdir(parents=True)

    def tearDown(self):
        self.temporary.cleanup()

    def test_ordered_parts_reconstruct_gzip_and_stream_comma_counts(self):
        content = b"Gene,cell-a,cell-b\ngene-1,1,0\ngene-2,0,3\n"
        compressed = gzip.compress(content)
        cut = len(compressed) // 2
        names = ["data/GSE123/raw/counts.csv.gz00", "data/GSE123/raw/counts.csv.gz01"]
        records = {}
        for name, payload in zip(names, (compressed[:cut], compressed[cut:])):
            path = self.root / name
            path.write_bytes(payload)
            records[name] = {"sha256": sha256(path)}
        destination = "data/GSE123/raw/counts.csv.gz"
        assembly = assemble_parts(self.root, "GSE123", destination, names, records)
        self.assertEqual(assembly["status"], "ASSEMBLED")
        self.assertEqual(gzip.decompress((self.root / destination).read_bytes()), content)
        records[destination] = assembly
        self.assertEqual(assemble_parts(self.root, "GSE123", destination, names, records)["status"], "REUSED")
        output = self.root / "triplets"
        result = convert(self.root / destination, output, "comma")
        self.assertEqual((result["genes"], result["cells"], result["nonzero"]), (2, 2, 2))
        self.assertEqual((output / "cells.tsv").read_text(), "cell-a\ncell-b\n")
        self.assertEqual(struct.unpack("<2d", (output / "x.bin").read_bytes()), (1.0, 3.0))

    def test_rejects_changed_part_and_incomplete_manifest(self):
        names = ["data/GSE123/raw/part0", "data/GSE123/raw/part1"]
        records = {}
        for name in names:
            path = self.root / name
            path.write_bytes(b"source")
            records[name] = {"sha256": sha256(path)}
        (self.root / names[1]).write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "verified download provenance"):
            assemble_parts(self.root, "GSE123", "data/GSE123/raw/counts.csv.gz", names, records)
        row = dict(database="GSE123", sample="GSM1", tissue="Stomach", disease="GC",
                   source_type="tumor", group="Tumor",
                   local_path="data/GSE123/raw/counts.csv.gz", file_type="text",
                   count_source="counts", count_evidence="author raw UMI",
                   metadata_evidence="author barcode suffix", orientation="genes_by_cells",
                   assembly_parts_json=json.dumps(names),
                   files_json=json.dumps([{"url": "https://example.org/part0", "local_path": names[0]}]))
        with self.assertRaisesRegex(ValueError, "declared download destinations"):
            validate_manifest([row], self.root)

    def test_tab_reader_still_accepts_existing_layout(self):
        source = self.raw / "counts.tsv.gz"
        source.write_bytes(gzip.compress(b"gene\tcell-a\ngene-1\t2\n"))
        result = convert(source, self.root / "tab-triplets")
        self.assertEqual((result["genes"], result["cells"], result["nonzero"]), (1, 1, 1))

    def test_truncated_http_fragment_is_never_finalized(self):
        class Response(io.BytesIO):
            headers = {"Content-Length": "10"}
        target = self.raw / "counts.csv.gz01"
        with patch("download_processed.urllib.request.urlopen", side_effect=lambda *_args, **_kwargs: Response(b"short")) as get, \
                patch("download_processed.time.sleep") as sleep:
            with self.assertRaisesRegex(OSError, "Incomplete HTTP response"):
                download("https://example.org/part01", target, expected_bytes=10)
        self.assertEqual(get.call_count, 4)
        self.assertEqual([call.args[0] for call in sleep.call_args_list], [2, 5, 10])
        self.assertFalse(target.exists())
        self.assertFalse((self.raw / "counts.csv.gz01.part").exists())

    def test_confirmation_accepts_persistent_cell_map(self):
        workflow = self.root / "data/GSE123/.workflow"
        workflow.mkdir(parents=True)
        cell_map = self.root / "data/GSE123/cell_map/cells.csv"
        cell_map.parent.mkdir(parents=True)
        cell_map.write_text("cell,sample\ncell-a,GSM1\n", encoding="utf-8")
        report = workflow / "sample_report.csv"
        write_csv(report, [dict(database="GSE123", sample="GSM1", tissue="Stomach",
                                disease="GC", source_type="tumor",
                                local_path="data/GSE123/raw/counts.csv.gz", file_type="text",
                                count_source="counts", count_evidence="author raw UMI",
                                metadata_evidence="barcode suffix",
                                cell_map_path="data/GSE123/cell_map/cells.csv",
                                orientation="genes_by_cells", delimiter="comma",
                                feature_column="__row_names__", files_json="[]")])
        (workflow / "inspection.json").write_text(
            json.dumps({"complete": True, "samples": [{"accession": "GSM1"}]}),
            encoding="utf-8")
        approval = workflow / "group_confirmation.json"
        approval.write_text(json.dumps({"confirmed_by": "user",
                                        "user_statement": "Tumor",
                                        "groups": {"GSM1": "Tumor"}}), encoding="utf-8")
        rows = confirm(report, approval, workflow / "sample_manifest.csv", self.root)
        self.assertTrue(rows[0]["cell_map_md5"])

    def test_large_comma_counts_select_streaming_automatically(self):
        rscript = shutil.which("Rscript")
        if not rscript:
            windows_r = Path("D:/R/R-4.5.3/bin/Rscript.exe")
            rscript = str(windows_r) if windows_r.is_file() else None
        if not rscript:
            self.skipTest("Rscript unavailable")
        source = (SCRIPTS / "seurat_common.R").as_posix()
        code = (
            f'source("{source}"); Sys.unsetenv("GEO_SINGLE_CELL_STREAM_TEXT"); '
            'row <- list(count_source="counts",delimiter="comma",'
            'orientation="genes_by_cells",feature_column="__row_names__",drop_columns=""); '
            'stopifnot(stream_text_candidate(row,174415491)); '
            'row$delimiter <- "tab"; stopifnot(!stream_text_candidate(row,174415491))'
        )
        subprocess.run([rscript, "-e", code], check=True, capture_output=True)


if __name__ == "__main__":
    unittest.main()
