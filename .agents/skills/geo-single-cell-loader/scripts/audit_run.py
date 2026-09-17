"""Publish one real GEO run audit, then remove temporary workflow files."""
import argparse
import csv
import json
import re
import shutil
import subprocess
import tempfile
import uuid
from datetime import datetime, timezone
from pathlib import Path

from common import ROOT, dataset_paths, read_csv, write_csv

AUDIT_REPO = "https://github.com/wdc19981006-cell/single-cell-skill-audit.git"
AUDIT_FILES = (
    "decision_trace.md", "stage_a.json", "file_selection.csv",
    "download_profile.csv", "input_routing.csv", "build_profile.json",
    "validation.json", "sample_summary.csv", "run_summary.txt", "execution.log",
)
DOWNLOAD_COLUMNS = (
    "local_path", "url", "member", "bytes", "sha256", "downloaded_at",
    "extracted_sha256", "extracted_bytes", "status", "download_seconds",
    "mb_per_second", "retry", "sha_seconds", "extract_seconds",
)
ROUTING_COLUMNS = (
    "local_path", "file_type", "reader", "unique_inputs", "expression_reads",
    "cell_map_reads", "features", "cells", "read_seconds", "mapping_seconds",
    "dense_conversion",
)


def git(*args, cwd=None):
    result = subprocess.run(["git", *args], cwd=cwd, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(f"git {' '.join(args)} failed: {result.stderr.strip()}")
    return result.stdout.strip()


def clone_audit():
    temporary = tempfile.TemporaryDirectory(prefix="geo-audit-")
    checkout = Path(temporary.name) / "repo"
    git("clone", "--quiet", AUDIT_REPO, str(checkout))
    return temporary, checkout


def run_id(workflow, require_clean=False):
    state = workflow / "audit_run.json"
    if state.exists():
        value = json.loads(state.read_text(encoding="utf-8"))
        if re.fullmatch(r"\d{8}T\d{6}Z-[0-9a-f]{8}", value["run_id"]):
            return value["run_id"]
        raise ValueError("Invalid audit run id")
    if require_clean and workflow.exists() and any(workflow.iterdir()):
        raise ValueError("Preexisting workflow lacks an audit run id; review and migrate it before starting a new real run")
    value = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid.uuid4().hex[:8]
    workflow.mkdir(parents=True, exist_ok=True)
    state.write_text(json.dumps({"run_id": value}, indent=2), encoding="utf-8")
    with (workflow / "execution.log").open("a", encoding="utf-8") as log:
        log.write(f"[{datetime.now(timezone.utc).isoformat()}] audit run started: {value}\n")
    return value


def restore_provenance(gse, root=ROOT):
    """Restore checksums for safe REUSED checks after previous local cleanup."""
    workflow = dataset_paths(root, gse)["workflow"]
    target = workflow / "download.json"
    if target.exists():
        return target
    temporary, checkout = clone_audit()
    try:
        records = {}
        parent = checkout / gse
        if parent.exists():
            for profile in sorted(parent.glob("*/download_profile.csv")):
                for row in read_csv(profile):
                    if row.get("local_path") and row.get("sha256") and row.get("url"):
                        records[row["local_path"]] = {
                            key: int(row[key]) if key in ("bytes", "extracted_bytes") and row.get(key) else row.get(key, "")
                            for key in ("local_path", "url", "member", "bytes", "sha256",
                                        "downloaded_at", "extracted_sha256", "extracted_bytes")
                        }
        workflow.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(list(records.values()), indent=2), encoding="utf-8")
        return target
    finally:
        temporary.cleanup()


def inventory(inspection):
    result = []
    for item in inspection.get("files", []):
        result.append(dict(url=item.get("url") or item.get("download_url") or item.get("source_url") or "",
                           owner_accession=item.get("owner_accession", item.get("scope", "")),
                           candidate_format=item.get("candidate_format", item.get("format", ""))))
    return result


def selected_files(workflow):
    manifest = workflow / "sample_manifest.csv"
    selected = {}
    if manifest.exists():
        for row in read_csv(manifest):
            for item in json.loads(row.get("files_json") or "[]"):
                selected[item["url"]] = item.get("local_path", "")
    return selected


def write_decisions(workflow, output, inspection, selection, status, reason):
    lines = ["# Reproducible decision trace", "", f"Status: {status}", f"Stop or exception reason: {reason or 'none'}", "",
             "## Discovery and MCP", f"Discovery method: {inspection.get('discovery_method', 'unavailable')}",
             f"MCP tools used: {', '.join(inspection.get('mcp_tools_used', [])) or 'none recorded'}",
             f"Fallback reason: {inspection.get('fallback_reason', 'none')}", "",
             "## Discovered and selected files"]
    for row in selection:
        lines.append(f"- {row['url']} | {row['selection']} | {row['reason']}")
    lines += ["", "## Raw count, structure, sample/GSM/cell mapping and reader evidence"]
    manifest = workflow / "sample_manifest.csv"
    if manifest.exists():
        for row in read_csv(manifest):
            lines += [f"- {row.get('sample', '')}: {row.get('local_path', '')}; format={row.get('file_type', '')}; "
                      f"count_source={row.get('count_source', '')}; raw_count_evidence={row.get('count_evidence', '')}; "
                      f"cell_map={row.get('cell_map_path', '')}; source_GSM={row.get('source_gsm', '')}; "
                      f"metadata_evidence={row.get('metadata_evidence', '')}"]
    else:
        lines.append("- Manifest unavailable; no expression reader or mapping decision was executed.")
    routing = workflow / "input_routing.csv"
    if routing.exists():
        lines += ["", "## Observed reader routing"]
        for row in read_csv(routing):
            lines.append(f"- {row.get('local_path', '')}: format={row.get('file_type', '')}; reader={row.get('reader', '')}; "
                         f"reason=verified input structure and manifest reader fields; "
                         f"dimensions={row.get('features', '')}x{row.get('cells', '')}; "
                         f"dense_conversion={row.get('dense_conversion', '')}")
    for name, title in (("text_schema_probe.json", "Text schema probe"),
                        ("pooled_mapping_probe.json", "Pooled mapping evidence check")):
        source = workflow / name
        if source.exists():
            lines += ["", f"## {title}", source.read_text(encoding="utf-8")]
    notes = workflow / "decision_notes.md"
    if notes.exists():
        lines += ["", "## Additional observed evidence", notes.read_text(encoding="utf-8")]
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def prepare_bundle(gse, root, status, reason, destination):
    paths = dataset_paths(root, gse)
    workflow = paths["workflow"]
    destination.mkdir(parents=True, exist_ok=True)
    inspection_file = workflow / "inspection.json"
    inspection = json.loads(inspection_file.read_text(encoding="utf-8")) if inspection_file.exists() else {}
    (destination / "stage_a.json").write_text(json.dumps(inspection or {"status": "unavailable"}, ensure_ascii=False, indent=2), encoding="utf-8")
    selected = selected_files(workflow)
    stated_reasons = {}
    reasons_file = workflow / "file_selection_reasons.csv"
    if reasons_file.exists():
        stated_reasons = {row["url"]: row.get("reason", "") for row in read_csv(reasons_file)}
    discovered = inventory(inspection)
    seen = {row["url"] for row in discovered}
    discovered.extend(dict(url=url, owner_accession="", candidate_format="") for url in selected if url not in seen)
    selection = [dict(**row, selection="SELECTED" if row["url"] in selected else "EXCLUDED",
                      local_path=selected.get(row["url"], ""),
                      reason=stated_reasons.get(row["url"]) or ("selected expression input" if row["url"] in selected else "not chosen; explicit reason unavailable"))
                 for row in discovered]
    write_csv(destination / "file_selection.csv", selection,
              ["url", "owner_accession", "candidate_format", "selection", "local_path", "reason"])
    write_decisions(workflow, destination / "decision_trace.md", inspection, selection, status, reason)
    for name, fields in (("download_profile.csv", DOWNLOAD_COLUMNS), ("input_routing.csv", ROUTING_COLUMNS),
                         ("sample_summary.csv", ("sample", "cells", "group"))):
        source = workflow / name
        if source.exists():
            shutil.copy2(source, destination / name)
        else:
            write_csv(destination / name, [], fields)
    for name, fallback in (
        ("build_profile.json", {"status": "unavailable", "reason": reason}),
        ("validation.json", {"status": status, "reason": reason}),
    ):
        source = workflow / name
        if source.exists():
            shutil.copy2(source, destination / name)
        else:
            (destination / name).write_text(json.dumps(fallback, ensure_ascii=False, indent=2), encoding="utf-8")
    summary = workflow / "run_summary.txt"
    (destination / "run_summary.txt").write_text(summary.read_text(encoding="utf-8") if summary.exists()
                                                   else f"GSE: {gse}\nStatus: {status}\nReason: {reason}\n", encoding="utf-8")
    log = workflow / "execution.log"
    (destination / "execution.log").write_text(log.read_text(encoding="utf-8") if log.exists()
                                                else f"Status: {status}\nReason: {reason}\n", encoding="utf-8")
    if set(p.name for p in destination.iterdir()) != set(AUDIT_FILES):
        raise ValueError("Audit bundle file set is incomplete")
    return selection


def preserve_essentials(gse, root, strict=True):
    paths = dataset_paths(root, gse)
    workflow, dataset = paths["workflow"], paths["dataset"]
    confirmation = workflow / "group_confirmation.json"
    if confirmation.exists():
        target = dataset / "group_confirmation.json"
        if target.exists() and target.read_bytes() != confirmation.read_bytes():
            if strict:
                raise ValueError("Existing group confirmation differs; refusing overwrite")
        else:
            shutil.copy2(confirmation, target)
    manifest = workflow / "sample_manifest.csv"
    if manifest.exists():
        for row in read_csv(manifest):
            relative = row.get("cell_map_path", "")
            if not relative:
                continue
            source = root / relative
            if source.is_file() and source.resolve().is_relative_to((dataset / "cell_map").resolve()):
                continue
            if not source.is_file() or not source.resolve().is_relative_to(workflow.resolve()):
                if strict:
                    raise ValueError("Required cell map is absent or outside workflow")
                continue
            target = dataset / "cell_map" / source.name
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.exists() and target.read_bytes() != source.read_bytes():
                if strict:
                    raise ValueError("Existing cell map differs; refusing overwrite")
            else:
                shutil.copy2(source, target)


def finalize(gse, root=ROOT, status="success", reason=""):
    if status not in ("success", "failure", "interrupted"):
        raise ValueError("Invalid run status")
    paths = dataset_paths(root, gse)
    workflow = paths["workflow"]
    if not workflow.exists():
        raise ValueError("No active workflow to audit")
    if status == "success":
        validation = workflow / "validation.json"
        if not paths["final"].is_file() or not validation.is_file() or json.loads(validation.read_text(encoding="utf-8")).get("status") != "success":
            raise ValueError("Successful audit requires final RDS and independent validation result")
    identifier = run_id(workflow)
    preserve_essentials(gse, root, strict=status == "success")
    temporary, checkout = clone_audit()
    try:
        destination = checkout / gse / identifier
        selection = prepare_bundle(gse, root, status, reason, destination)
        if status == "success" and any(not row["url"] for row in selection):
            raise ValueError("Discovered file lacks a source URL; complete Stage A inventory before successful finalization")
        if status == "success" and any(row["reason"] == "not chosen; explicit reason unavailable" for row in selection):
            raise ValueError("Record an explicit exclusion reason in file_selection_reasons.csv before successful finalization")
        git("add", "--", f"{gse}/{identifier}", cwd=checkout)
        git("-c", "user.name=single-cell-skill-audit", "-c", "user.email=audit@users.noreply.github.com",
            "commit", "-m", f"{gse} {identifier} {status}", cwd=checkout)
        git("push", "origin", "HEAD", cwd=checkout)
        # A successful push is the deletion boundary. Keep all files if upload fails.
        resolved = workflow.resolve()
        expected = (root.resolve() / "data" / gse / ".workflow").resolve()
        if resolved != expected or not resolved.is_relative_to(root.resolve()):
            raise ValueError("Refusing to remove workflow outside the repository")
        shutil.rmtree(resolved)
        return f"{AUDIT_REPO[:-4]}/tree/main/{gse}/{identifier}"
    finally:
        temporary.cleanup()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    for action in ("start", "restore-provenance", "finalize"):
        command = sub.add_parser(action)
        command.add_argument("gse")
        command.add_argument("--root", type=Path, default=ROOT)
        if action == "finalize":
            command.add_argument("--status", choices=("success", "failure", "interrupted"), required=True)
            command.add_argument("--reason", default="")
    args = parser.parse_args()
    root = args.root.resolve()
    if args.action == "start":
        print(run_id(dataset_paths(root, args.gse)["workflow"], require_clean=True))
    elif args.action == "restore-provenance":
        print(restore_provenance(args.gse, root))
    else:
        print(finalize(args.gse, root, args.status, args.reason))


if __name__ == "__main__":
    main()
