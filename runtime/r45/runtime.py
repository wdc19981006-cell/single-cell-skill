"""Shared Python caller for the repository's sole R launcher."""
import subprocess
import sys
import json
import uuid
from pathlib import Path
from config import RUNTIME_LOG

LAUNCHER = Path(__file__).with_name("run_r45.py")


class RRunError(RuntimeError):
    def __init__(self, stage, script, exit_code, category, output=""):
        self.stage, self.script, self.exit_code, self.category = stage, str(script), exit_code, category
        self.output = output
        detail = "\n" + output[-1500:] if output else ""
        super().__init__(f"{category}: {stage} failed with exit code {exit_code}: {script}{detail}")


def run_r45(script=None, args=(), stage="adhoc", log_path=None, cwd=None, route=None, env=None, expression=None):
    record_path = RUNTIME_LOG / "requests" / (uuid.uuid4().hex + ".json")
    command = [sys.executable, str(LAUNCHER), "--stage", stage, "--record", str(record_path)]
    if route:
        command += ["--route", route]
    if expression is not None:
        if script is not None:
            raise ValueError("Choose a script or an expression")
        command += ["--expr", expression]
    else:
        command += ["--healthcheck"] if script is None else [str(script), *map(str, args)]
    result = subprocess.run(command, cwd=cwd, env=env, capture_output=True, text=True, encoding="utf-8", errors="replace")
    record = json.loads(record_path.read_text(encoding="utf-8")) if record_path.exists() else {}
    result.record = record
    output = (result.stdout or "") + (result.stderr or "")
    if log_path is not None:
        if hasattr(log_path, "write"):
            log_path.write(output)
            log_path.flush()
        else:
            with Path(log_path).open("a", encoding="utf-8") as log:
                log.write(output)
    if result.returncode:
        category = record.get("status") or ("R_NATIVE_CRASH" if result.returncode in
                    (139, -1073741819, 3221225477) else "R_SCRIPT_ERROR")
        raise RRunError(stage, script or "healthcheck.R", result.returncode, category,
                        output)
    return result
