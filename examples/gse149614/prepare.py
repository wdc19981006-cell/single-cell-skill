"""Prepare and verify GSE149614's pooled cell-to-GSM mapping.

Run after MCP-first discovery has created the canonical inspection and report.
``prepare`` requires the user's exact group confirmation; it never chooses groups.
"""

import argparse
import csv
import gzip
import io
import json
from collections import Counter
from pathlib import Path


GSE = "GSE149614"
EXPECTED_SAMPLES = 21
EXPECTED_CELLS = 71915
COUNT_NAME = "GSE149614_HCC.scRNAseq.S71915.count.txt.gz"
COUNT_URL = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/{GSE}/suppl/{COUNT_NAME}"
METADATA_NAME = "GSE149614_HCC.metadata.updated.txt.gz"
METADATA_URL = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/{GSE}/suppl/{METADATA_NAME}"
SITE_FROM_GEO_TISSUE = {
    "primary tumor": "Tumor",
    "adjacent non-tumor liver": "Normal",
    "portal vein tumor thrombus (PVTT)": "PVTT",
    "metastatic lymph node": "Lymph",
}
EXTRA_FIELDS = (
    "local_path", "file_type", "count_source", "count_evidence", "delimiter",
    "orientation", "feature_column", "drop_columns", "cell_map_path",
    "cell_map_cell_column", "cell_map_sample_column", "files_json",
)


def load_inspection(path):
    inspection = json.loads(path.read_text(encoding="utf-8"))
    samples = inspection.get("samples", [])
    if inspection.get("gse") != GSE or inspection.get("complete") is not True or len(samples) != EXPECTED_SAMPLES:
        raise ValueError("A complete GSE149614 inspection with 21 GSM records is required")
    by_title = {}
    gsms = set()
    for record in samples:
        gsm = record.get("accession") or record.get("gsm")
        title = record.get("title")
        traits = {item.get("key"): item.get("value") for item in record.get("characteristics", [])}
        if not gsm or not title or gsm in gsms or title in by_title:
            raise ValueError("Missing or duplicate GEO GSM/title mapping")
        site = SITE_FROM_GEO_TISSUE.get(traits.get("tissue"))
        if not site or not traits.get("patient") or traits.get("tumor") != "Hepatocellular carcinoma (HCC)":
            raise ValueError(f"Incomplete clinical source facts for {gsm}")
        by_title[title] = {"gsm": gsm, "patient": traits["patient"], "site": site}
        gsms.add(gsm)
    return by_title


def read_metadata(path, by_title, expected_cells=EXPECTED_CELLS):
    seen = set()
    mapping = []
    counts = Counter()
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Cell", "sample", "site", "patient"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Public metadata lacks Cell/sample/site/patient columns")
        for row in reader:
            title = row["sample"]
            if title not in by_title:
                raise ValueError(f"Unknown author sample in metadata: {title}")
            facts = by_title[title]
            cell = row["Cell"]
            if not cell or cell in seen or not cell.startswith(title + "_"):
                raise ValueError(f"Blank, duplicate, or mismatched cell ID: {cell}")
            if row["patient"] != facts["patient"] or row["site"] != facts["site"]:
                raise ValueError(f"Patient/site mismatch for {title}")
            seen.add(cell)
            mapping.append((cell, facts["gsm"]))
            counts[facts["gsm"]] += 1
    if len(mapping) != expected_cells or set(counts) != {item["gsm"] for item in by_title.values()}:
        raise ValueError("Metadata cell or sample inventory is incomplete")
    return mapping, counts


def read_report(path, by_title):
    original = path.read_text(encoding="utf-8-sig")
    reader = csv.DictReader(io.StringIO(original))
    fields = list(reader.fieldnames or [])
    rows = list(reader)
    if len(rows) != EXPECTED_SAMPLES or {row["sample"] for row in rows} != {item["gsm"] for item in by_title.values()}:
        raise ValueError("Stage A sample report does not cover all 21 GSM records")
    for row in rows:
        facts = by_title.get(row.get("author_sample"))
        if not facts or row["sample"] != facts["gsm"] or row.get("patient") != facts["patient"]:
            raise ValueError("Report author-sample/GSM/patient mapping differs from GEO")
        if row.get("disease") != "HCC":
            raise ValueError("All patients have HCC, including adjacent non-tumor samples")
    return rows, fields, original


def check_confirmation(path, samples):
    if not path.is_file():
        raise ValueError("group 尚未创建，请先由用户确认全部 21 个 GSM 的分组")
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    groups = value.get("groups")
    if (value.get("confirmed_by") != "user" or not str(value.get("user_statement", "")).strip()
            or not isinstance(groups, dict) or set(groups) != samples
            or any(not isinstance(group, str) or not group.strip() for group in groups.values())):
        raise ValueError("Explicit user confirmation must cover every GSM with a nonblank group")


def cell_map_text(mapping):
    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(("cell", "sample"))
    writer.writerows(mapping)
    return output.getvalue()


def prepare_report(rows, fields):
    if any(field in fields for field in EXTRA_FIELDS):
        raise ValueError("Report already has a download plan; do not silently rewrite it")
    local_path = f"data/{GSE}/raw/{COUNT_NAME}"
    cell_map_path = f"data/{GSE}/.workflow/cell_map.csv"
    plan = json.dumps([{"url": COUNT_URL, "local_path": local_path}], separators=(",", ":"))
    for row in rows:
        row.update({
            "local_path": local_path,
            "file_type": "text",
            "count_source": "counts",
            "count_evidence": (
                f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={row['sample']} "
                "Sample_data_processing identifies the pooled .count.txt.gz as Cell Ranger raw counts; "
                "the separate normalized/log-transformed file is excluded"
            ),
            "delimiter": "tab",
            "orientation": "genes_by_cells",
            "feature_column": "__row_names__",
            "drop_columns": "",
            "cell_map_path": cell_map_path,
            "cell_map_cell_column": "cell",
            "cell_map_sample_column": "sample",
            "files_json": plan,
        })
        row["metadata_evidence"] += f"; {METADATA_URL} Cell/sample/patient/site; exact author sample to GSM from GEO Sample_title"
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields + list(EXTRA_FIELDS), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def write_once(path, content):
    if path.exists():
        if path.read_text(encoding="utf-8") != content:
            raise ValueError(f"Existing file differs; review before changing: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".part")
    temporary.write_text(content, encoding="utf-8", newline="")
    temporary.replace(path)


def replace_report(path, original, content):
    if path.read_text(encoding="utf-8-sig") != original:
        raise ValueError("Stage A report changed during preparation; review before retrying")
    temporary = path.with_name(path.name + ".part")
    temporary.write_text(content, encoding="utf-8", newline="")
    temporary.replace(path)


def verify_download(path, mapping):
    expected = {cell for cell, _ in mapping}
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if len(header) != len(expected) or len(set(header)) != len(header) or set(header) != expected:
            raise ValueError("Count matrix header does not equal the metadata cell set")
        genes = 0
        for line in handle:
            if line.count("\t") != len(header):
                raise ValueError(f"Count matrix row {genes + 1} has the wrong number of columns")
            genes += 1
    if genes == 0:
        raise ValueError("Count matrix contains no genes")
    return genes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("check", "prepare", "verify-download"))
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    root = args.root.resolve()
    workflow = root / "data" / GSE / ".workflow"
    by_title = load_inspection(workflow / "inspection.json")
    mapping, counts = read_metadata(workflow / "public_metadata" / METADATA_NAME, by_title)
    rows, fields, original_report = read_report(workflow / "sample_report.csv", by_title)
    if args.action == "check":
        print(json.dumps({"gse": GSE, "samples": len(counts), "cells": len(mapping), "per_sample": counts}, indent=2))
        return
    check_confirmation(workflow / "group_confirmation.json", {row["sample"] for row in rows})
    if args.action == "prepare":
        if (workflow / "sample_manifest.csv").exists():
            raise ValueError("Manifest already exists; do not silently rebuild a confirmed workflow")
        report_text = prepare_report(rows, fields)
        write_once(workflow / "cell_map.csv", cell_map_text(mapping))
        replace_report(workflow / "sample_report.csv", original_report, report_text)
        print("Prepared exact cell map and pooled-count report for confirmed groups")
    else:
        genes = verify_download(root / "data" / GSE / "raw" / COUNT_NAME, mapping)
        print(json.dumps({"gse": GSE, "genes": genes, "cells": len(mapping), "cell_ids_match": True}))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        raise SystemExit(str(error)) from None
