"""Launch every project R process with the fixed R 4.5.0 profile and audit it."""
import argparse
import json
import os
import signal
import subprocess
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

from config import PROFILE, REPO, R_HOME, R_LIBRARY, RSCRIPT, RUNTIME_LOG, R_VERSION

NATIVE_CRASH = {139, -1073741819, 3221225477}


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    pending = path.with_name(path.name + "." + uuid.uuid4().hex + ".pending")
    pending.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")
    pending.replace(path)


def run(script=None, args=(), *, expression=None, stage="adhoc", route=None, cwd=None, record_path=None):
    if script is None and expression is None:
        raise ValueError("Provide an R script or --expr")
    if script is not None and expression is not None:
        raise ValueError("Choose an R script or --expr")
    if any(value in ("--vanilla", "--no-init-file", "--no-environ") for value in args):
        raise ValueError("R options that bypass the project profile are forbidden")
    cwd = Path(cwd or Path.cwd()).resolve()
    script_path = (cwd / script).resolve() if script else None
    command = [str(RSCRIPT), "--no-environ", "--no-site-file", "--no-restore", "--no-save"]
    command += ["-e", expression] if expression is not None else [str(script_path), *map(str, args)]
    env = os.environ.copy()
    for key in ("LC_ALL", "LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_NUMERIC", "LC_TIME",
                "R_ENVIRON_USER", "R_PROFILE", "R_LIBS_SITE"):
        env.pop(key, None)
    env.update(R_HOME=str(R_HOME), R_PROFILE_USER=str(PROFILE), CLI_NO_THREAD="1",
               R_LIBS=str(R_LIBRARY), R_LIBS_USER=str(R_LIBRARY), R_LIBS_SITE=str(R_LIBRARY), LC_ALL="C",
               R45_EXPECTED_HOME=R_HOME.as_posix(), R45_EXPECTED_LIBRARY=R_LIBRARY.as_posix(),
               R45_EXPECTED_VERSION=R_VERSION)
    record = dict(r_executable=str(RSCRIPT), r_home=str(R_HOME), r_libs=env["R_LIBS"],
                  r_libs_user=env["R_LIBS_USER"], path=env.get("PATH", ""),
                  script=str(script_path) if script_path else "--expr", args=list(map(str, args)) if script_path else [expression],
                  stage=stage, route=route, cwd=str(cwd), started_at=timestamp(), command=command)
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S") + "-" + uuid.uuid4().hex[:8]
    run_path = RUNTIME_LOG / "runs" / f"{run_id}.json"
    record.update(status="RUNNING", run_id=run_id, cli_exit_workaround=True,
                  locale=env["LC_ALL"], profile=str(PROFILE))
    atomic_json(run_path, record)
    atomic_json(RUNTIME_LOG / "latest.json", record)
    started = time.monotonic()
    if not RSCRIPT.is_file():
        record.update(status="R_RUNTIME_MISSING", exit_code=2,
                      stdout="", stderr=f"R_RUNTIME_MISSING: {RSCRIPT}")
    elif not PROFILE.is_file():
        record.update(status="R_SCRIPT_ERROR", exit_code=2,
                      stdout="", stderr=f"Required runtime profile missing: {PROFILE}")
    elif script_path and not script_path.is_file():
        record.update(status="R_SCRIPT_ERROR", exit_code=2,
                      stdout="", stderr=f"R script missing: {script_path}")
    else:
        try:
            stdout_path = run_path.with_suffix(".stdout.log")
            stderr_path = run_path.with_suffix(".stderr.log")
            with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
                child = subprocess.Popen(command, cwd=cwd, env=env, stdout=stdout, stderr=stderr)
                try:
                    code = child.wait()
                except KeyboardInterrupt:
                    if os.name == "nt":
                        subprocess.run(["taskkill", "/PID", str(child.pid), "/T", "/F"], capture_output=True)
                    else:
                        child.terminate()
                    child.wait()
                    code = 130
            output = stdout_path.read_text(encoding="utf-8", errors="replace")
            errors = stderr_path.read_text(encoding="utf-8", errors="replace")
            status = "PASS" if code == 0 else (
                "R_NATIVE_CRASH" if code in NATIVE_CRASH else
                "INTERRUPTED" if code == 130 else
                "R_PACKAGE_MISSING" if "R_PACKAGE_MISSING" in errors else "R_SCRIPT_ERROR")
            record.update(status=status, exit_code=code, stdout=output, stderr=errors)
        except OSError as error:
            record.update(status="R_RUNTIME_MISSING", exit_code=2, stdout="", stderr=str(error))
    record.update(ended_at=timestamp(), duration_seconds=round(time.monotonic() - started, 3))
    atomic_json(run_path, record)
    atomic_json(RUNTIME_LOG / "latest.json", record)
    if record_path:
        atomic_json(Path(record_path), record)
    return record


def health_record(record):
    packages = {}
    details = {}
    for line in record["stdout"].splitlines():
        fields = line.split("\t")
        if fields[0] == "R45_PACKAGE" and len(fields) == 3:
            packages[fields[1]] = fields[2]
        elif fields[0].startswith("R45_") and len(fields) == 2:
            details[fields[0][4:].lower()] = fields[1]
    health = dict(status=record["status"], r_version=details.get("short_version"),
                  r_version_string=details.get("version"),
                  r_home=details.get("home"), library=R_LIBRARY.as_posix(),
                  lib_paths=details.get("libraries", "").split(";"), packages=packages,
                  route=details.get("route"), cli_exit_workaround=record["status"] == "PASS",
                  checked_at=record["ended_at"], exit_code=record["exit_code"])
    atomic_json(RUNTIME_LOG / "health.json", health)


def main():
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")
    signal.signal(signal.SIGTERM, lambda *_: (_ for _ in ()).throw(KeyboardInterrupt()))
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--expr", "-e")
    mode.add_argument("--healthcheck", action="store_true")
    parser.add_argument("--route", choices=("base", "qc", "text", "10x_h5", "h5ad_native"), default="base")
    parser.add_argument("--stage", default="adhoc")
    parser.add_argument("--record", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("script", nargs="?")
    parser.add_argument("args", nargs=argparse.REMAINDER)
    options = parser.parse_args()
    if options.healthcheck:
        if options.script or options.args:
            parser.error("--healthcheck does not accept a script")
        record = run(Path(__file__).with_name("healthcheck.R"), [options.route], stage="healthcheck", route=options.route, record_path=options.record)
        health_record(record)
    elif options.expr is not None:
        if options.script or options.args:
            parser.error("--expr does not accept a script")
        record = run(expression=options.expr, stage=options.stage, record_path=options.record)
    else:
        if not options.script:
            parser.error("An R script is required")
        record = run(options.script, options.args, stage=options.stage, record_path=options.record)
    if record["stdout"]:
        print(record["stdout"], end="", flush=True)
    if record["stderr"]:
        print(record["stderr"], end="", file=sys.stderr, flush=True)
    if record["status"] != "PASS":
        print(f"{record['status']}: stage={record['stage']} script={record['script']} exit={record['exit_code']}",
              file=sys.stderr, flush=True)
    return record["exit_code"]


if __name__ == "__main__":
    sys.exit(main())
