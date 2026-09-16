"""Testable Stage A orchestration policy; the Codex agent supplies real MCP callables."""
from __future__ import annotations

import re
import time
from datetime import datetime, timezone
from typing import Any, Callable


class IncompleteMCPResult(ValueError):
    """Raised when an MCP response is successful at transport level but incomplete."""


def _urls(value: Any) -> list[str]:
    if isinstance(value, dict):
        return [url for child in value.values() for url in _urls(child)]
    if isinstance(value, (list, tuple)):
        return [url for child in value for url in _urls(child)]
    if isinstance(value, str) and value.startswith(("https://", "http://", "ftp://")):
        return [value]
    return []


def validate_series(gse: str, result: dict[str, Any]) -> None:
    accession = result.get("accession") or result.get("gse")
    samples = result.get("samples")
    if accession != gse:
        raise IncompleteMCPResult("geo/get_geo_info returned a different accession")
    if result.get("complete") is not True:
        raise IncompleteMCPResult("geo/get_geo_info returned complete=false")
    if not isinstance(samples, list) or result.get("sample_count") != len(samples):
        raise IncompleteMCPResult("geo/get_geo_info sample_count does not match the GSM inventory")
    accessions = [sample.get("accession") for sample in samples if isinstance(sample, dict)]
    if len(accessions) != len(samples) or any(not re.fullmatch(r"GSM\d+", accession or "") for accession in accessions) or len(set(accessions)) != len(accessions):
        raise IncompleteMCPResult("geo/get_geo_info returned an invalid GSM inventory")
    if not _urls(result):
        raise IncompleteMCPResult("geo/get_geo_info returned no source URLs")


def validate_files(gse: str, result: dict[str, Any]) -> None:
    files = result.get("files")
    if not isinstance(files, list):
        raise IncompleteMCPResult("geo/list_geo_files returned no files list")
    for entry in files:
        if not isinstance(entry, dict) or not entry.get("owner_accession") or not entry.get("source_level") or not entry.get("candidate_format"):
            raise IncompleteMCPResult("geo/list_geo_files returned a file without ownership, source level, or candidate format")
        if not str(entry["owner_accession"]).startswith((gse, "GSM")):
            raise IncompleteMCPResult("geo/list_geo_files returned an unrelated owner accession")


def discover_known_gse(
    gse: str,
    get_geo_info: Callable[[str], dict[str, Any]],
    list_geo_files: Callable[[str], dict[str, Any]],
    fallback: Callable[[str, str], Any],
) -> dict[str, Any]:
    """Call get_geo_info then list_geo_files; fallback exactly once only on real failure."""
    if not re.fullmatch(r"GSE\d+", gse):
        raise ValueError("Expected GSE accession")
    started_at = datetime.now(timezone.utc)
    timer_started = time.perf_counter()
    tools_used: list[str] = []
    try:
        series = get_geo_info(gse)
        tools_used.append("get_geo_info")
        validate_series(gse, series)
        files = list_geo_files(gse)
        tools_used.append("list_geo_files")
        validate_files(gse, files)
    except Exception as error:
        reason = f"geo MCP failure: {type(error).__name__}: {error}"
        finished_at = datetime.now(timezone.utc)
        return {
            "discovery_method": "inspect_geo_fallback",
            "mcp_tools_used": tools_used,
            "mcp_complete": False,
            "fallback_reason": reason,
            "fallback_result": fallback(gse, reason),
            "discovery_started_at": started_at.isoformat(),
            "discovery_finished_at": finished_at.isoformat(),
            "discovery_seconds": time.perf_counter() - timer_started,
        }
    finished_at = datetime.now(timezone.utc)
    return {
        "gse": gse,
        "series": series,
        "samples": series["samples"],
        "files": files["files"],
        "input_plans": series.get("input_plans", {}),
        "papers": series.get("papers", series.get("publications", [])),
        "total_samples": series["sample_count"],
        "inspected_samples": len(series["samples"]),
        "complete": True,
        "source": "global geo MCP",
        "discovery_method": "geo_mcp",
        "mcp_tools_used": tools_used,
        "mcp_complete": True,
        "discovery_started_at": started_at.isoformat(),
        "discovery_finished_at": finished_at.isoformat(),
        "discovery_seconds": time.perf_counter() - timer_started,
    }


def discover_keyword(
    query: str,
    search_geo: Callable[[str], Any],
    select_gse: Callable[[Any], str],
    get_geo_info: Callable[[str], dict[str, Any]],
    list_geo_files: Callable[[str], dict[str, Any]],
    fallback: Callable[[str, str], Any],
) -> dict[str, Any]:
    """Keyword discovery always searches MCP before applying the known-GSE flow."""
    started_at = datetime.now(timezone.utc)
    timer_started = time.perf_counter()
    try:
        candidates = search_geo(query)
        gse = select_gse(candidates)
    except Exception as error:
        reason = f"geo MCP failure: {type(error).__name__}: {error}"
        finished_at = datetime.now(timezone.utc)
        return {
            "discovery_method": "inspect_geo_fallback",
            "mcp_tools_used": ["search_geo"],
            "mcp_complete": False,
            "fallback_reason": reason,
            "fallback_result": fallback(query, reason),
            "discovery_started_at": started_at.isoformat(),
            "discovery_finished_at": finished_at.isoformat(),
            "discovery_seconds": time.perf_counter() - timer_started,
        }
    result = discover_known_gse(gse, get_geo_info, list_geo_files, fallback)
    result["mcp_tools_used"] = ["search_geo", *result.get("mcp_tools_used", [])]
    result["discovery_started_at"] = started_at.isoformat()
    result["discovery_finished_at"] = datetime.now(timezone.utc).isoformat()
    result["discovery_seconds"] = time.perf_counter() - timer_started
    return result
