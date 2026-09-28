"""Stage checkpoints and error records shared by raw and QC workflows."""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from run_r45 import atomic_json


def checkpoint(workflow, **changes):
    path = Path(workflow) / "state.json"
    state = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    state.update(changes)
    state["updated_at"] = datetime.now(timezone.utc).isoformat()
    atomic_json(path, state)
    return state


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def failure(error, stage, script="", category="R_SCRIPT_ERROR"):
    return dict(category=getattr(error, "category", category), stage=stage,
                script=getattr(error, "script", str(script)),
                exit_code=getattr(error, "exit_code", None), message=str(error),
                timestamp=datetime.now(timezone.utc).isoformat())
