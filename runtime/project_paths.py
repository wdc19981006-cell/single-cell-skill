"""Shared storage contract; --root retains explicit legacy workspace compatibility."""
import json
import os
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_ROOT = Path(json.loads(Path(__file__).with_name("output_config.json").read_text())["output_root"])


def dataset_name(value):
    if (not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,99}", value)
            or value.endswith((".", " ")) or value.split(".")[0].upper() in
            {"CON", "PRN", "AUX", "NUL", *(f"COM{i}" for i in range(1, 10)), *(f"LPT{i}" for i in range(1, 10))}):
        raise ValueError("Unsafe dataset/run name")
    return value


def storage_root(root=None, output_root=None):
    if output_root is not None or root is None:
        target = Path(output_root or os.environ.get("SINGLE_CELL_OUTPUT_ROOT") or DEFAULT_OUTPUT_ROOT).resolve()
        os.environ["SINGLE_CELL_OUTPUT_ROOT"] = str(target)
        return target
    return Path(root).resolve()


def flat_storage(root):
    root = Path(root).resolve()
    declared = os.environ.get("SINGLE_CELL_OUTPUT_ROOT")
    return root == DEFAULT_OUTPUT_ROOT.resolve() or bool(declared and root == Path(declared).resolve())


def project_path(root, name):
    name = dataset_name(name)
    root = Path(root).resolve()
    path = root / name if flat_storage(root) else root / "results" / name
    if path.resolve() != path or not path.is_relative_to(root):
        raise ValueError("Dataset path redirects outside declared output location")
    return path


def analysis_path(root, name, run_id=None):
    path = project_path(root, name) / "03_analysis"
    if run_id is not None:
        path = path / "runs" / dataset_name(run_id)
    if path.resolve() != path:
        raise ValueError("Analysis path redirects outside declared output location")
    return path


def runtime_output(root, name):
    os.environ["SINGLE_CELL_RUNTIME_ROOT"] = str(project_path(root, name) / ".runtime/r45")
