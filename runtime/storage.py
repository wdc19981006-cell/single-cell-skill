"""File-only lifecycle decisions. Never open an RDS to decide retention.

Only enrolled, newly created outputs are owned. Existing datasets are inventory-only.
Cleanup is a write-ahead transaction: manifests record a verified successor before
any unlink. An interrupted unlink is safe to retry; absent owned files are never
silently treated as corrupt or as usable inputs.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
from datetime import datetime, timezone
from uuid import uuid4
from contextlib import contextmanager

CONFIG = Path(__file__).with_name("storage_config.json")
RETAINED, CLEANED, INCOMPLETE = "COMPLETE_RETAINED", "COMPLETE_CLEANED", "INCOMPLETE"


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def atomic(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + "." + uuid4().hex + ".pending")
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")
    os.replace(tmp, path)


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


class Storage:
    def __init__(self, dataset, mode=None, retain=None):
        self.dataset = Path(dataset).resolve()
        self.work = self.dataset / ".work/checkpoints"
        self.policy_path = self.dataset / ".work/storage_policy.json"
        self.registry_path = self.work / "registry.json"
        self.active = self.policy_path.is_file()
        cfg = read(CONFIG)
        if mode is not None and mode not in ("compact", "full"):
            raise ValueError("Storage mode must be compact or full")
        if not self.active:
            # Raw-only datasets can enroll. Any pre-existing downstream object
            # (including an alternate run) makes this an old, protected dataset.
            old = any(p.is_file() for stage in ("02_QC", "03_analysis")
                      for p in (self.dataset / stage).rglob("*.rds"))
            if self.dataset.parent.name=="results":
                old = old or any((self.dataset.parents[1]/"data"/self.dataset.name/"qc").glob("*.rds"))
            if not old:
                self.dataset.mkdir(parents=True, exist_ok=True)
                atomic(self.policy_path, dict(cfg, mode=mode or cfg["mode"],
                    retain_checkpoints=list(retain or cfg["retain_checkpoints"]),
                    created_at=datetime.now(timezone.utc).isoformat(), ownership="new_outputs_only"))
                self.active = True
        self.policy = read(self.policy_path) if self.active else dict(cfg, mode="protected_legacy")
        if self.policy["mode"] not in ("compact","full","protected_legacy"):
            raise ValueError("Invalid saved storage mode")
        if retain is not None and (not isinstance(retain,(list,tuple)) or any(not isinstance(v,str) or not v for v in retain)):
            raise ValueError("retain_checkpoints must be a list of stage names or basenames")
        if self.active and (mode is not None or retain is not None):
            if mode is not None: self.policy["mode"] = mode
            if retain is not None: self.policy["retain_checkpoints"] = list(retain)
            atomic(self.policy_path, self.policy)
        self.registry = read(self.registry_path) if self.registry_path.is_file() else {"files": {}, "references": {}}

    def save(self):
        if self.active: atomic(self.registry_path, self.registry)

    def checkpoint_root(self, scope):
        path = self.work / scope
        if path.resolve() != path or not path.is_relative_to(self.dataset):
            raise ValueError("Checkpoint path redirected")
        return path

    def reference(self, role, path):
        path = Path(path).resolve()
        value = {"path": str(path), "sha256": sha256(path), "size_bytes": path.stat().st_size,
                 "external": not path.is_relative_to(self.dataset), "never_delete": True}
        previous = self.registry["references"].get(role)
        if previous and previous != value:
            raise ValueError(f"Referenced input changed: {path}; use a new run")
        self.registry["references"][role] = value
        self.save()

    def referenced(self, role):
        value = self.registry["references"].get(role)
        if not value: return None
        path = Path(value["path"])
        if not path.is_file() or sha256(path) != value["sha256"]:
            raise ValueError(f"Referenced input unavailable or changed: {path}")
        return path

    def own(self, path, stage, scope, formal=False, manifest=None):
        if not self.active: return
        path = Path(path).absolute()
        if path.resolve() != path or not path.is_relative_to(self.dataset):
            raise ValueError("Refusing ownership of redirected/external output")
        # Only temporary RDS inside the checkpoint tree can ever be deleted.
        temporary = not formal and path.is_relative_to(self.work)
        self.registry["files"][str(path)] = dict(path=str(path), stage=stage, scope=scope,
            size_bytes=path.stat().st_size, sha256=sha256(path), formal=formal,
            temporary=temporary, lifecycle=RETAINED, manifest=str(manifest) if manifest else None)
        self.save()

    def incomplete(self, directory, stage, scope):
        if not self.active: return
        for path in Path(directory).rglob("*.rds"):
            if str(path) in self.registry["files"]: continue
            self.own(path,stage,scope)
            self.registry["files"][str(path)]["lifecycle"] = INCOMPLETE
        self.save()

    def release_staging(self):
        # Publication transferred these bytes to an owned, retained checkpoint.
        self.registry["files"] = {k:v for k,v in self.registry["files"].items()
            if not (v["lifecycle"]==INCOMPLETE and not Path(k).exists())}
        self.save()

    def receipt(self, final, stage, scope, validator="independent_R_readback"):
        final = Path(final).resolve()
        proof = dict(path=str(final), sha256=sha256(final), stage=stage, scope=scope,
                     validator=validator, validated=True)
        atomic(self.checkpoint_root(scope) / "formal_validation.json", proof)
        return proof

    @staticmethod
    def valid_proof(proof):
        path = Path(proof.get("path", ""))
        return (proof.get("validated") is True and proof.get("validator")=="independent_R_readback"
                and path.is_file() and sha256(path)==proof.get("sha256"))

    def proof(self, scope):
        path = self.checkpoint_root(scope) / "formal_validation.json"
        return read(path) if path.is_file() else None

    def clean(self, scope, proof):
        if not self.active or self.policy["mode"] != "compact": return
        if not self.valid_proof(proof): raise ValueError("Cleanup requires a validated, unchanged downstream final")
        formal = self.registry["files"].get(proof["path"],{})
        if proof.get("scope")!=scope or not formal.get("formal") or formal.get("scope")!=scope:
            raise ValueError("Cleanup successor belongs to a different run or is not an owned formal object")
        keep = set(self.policy["retain_checkpoints"])
        referenced = {v["path"] for v in self.registry["references"].values() if v.get("never_delete")}
        for item in self.registry["files"].values():
            if item["scope"]==scope and item["temporary"] and item["lifecycle"]==RETAINED and not Path(item["path"]).is_file():
                raise ValueError(f"Unexplained missing checkpoint: {item['path']}")
        candidates = [v for v in self.registry["files"].values() if v["scope"]==scope
            and v["temporary"] and not v["formal"] and v["lifecycle"]!=INCOMPLETE
            and v["path"] not in referenced
            and not (v["stage"] in keep or f"{scope}:{v['stage']}" in keep or Path(v["path"]).name in keep)
            and (Path(v["path"]).exists() or v.get("cleanup_pending"))]
        if not candidates: return
        event = dict(scope=scope,started_at=datetime.now(timezone.utc).isoformat(),
                     total_before_bytes=self.total_bytes(),removed_checkpoint_bytes=0)
        for key, item in self.registry["files"].items():
            if item["scope"] != scope or not item["temporary"] or item["formal"]: continue
            if key in referenced: continue
            if item["lifecycle"]==INCOMPLETE: continue
            if item["stage"] in keep or f"{scope}:{item['stage']}" in keep or Path(key).name in keep: continue
            path = Path(key)
            if path.resolve() != path or not path.is_relative_to(self.work):
                raise ValueError("Cleanup path redirected or outside owned checkpoints")
            if path.is_file() and sha256(path) != item["sha256"]:
                raise ValueError(f"Checkpoint changed; preserved: {path}")
            if not path.exists() and item["lifecycle"] != CLEANED:
                raise ValueError(f"Unexplained missing checkpoint: {path}")
            # Persist tombstone before unlink. COMPLETE_CLEANED means the old
            # object is no longer an input, even if interruption leaves bytes.
            if item.get("manifest"):
                manifest = Path(item["manifest"])
                record = read(manifest)
                checksum = record.get("outputs", {}).pop(key, None)
                if checksum:
                    record.setdefault("cleaned_outputs", {})[key] = dict(sha256=checksum, successor=proof)
                record["lifecycle"] = CLEANED
                atomic(manifest, record)
            item.update(lifecycle=CLEANED, successor=proof, cleanup_pending=path.exists())
            self.save()
            if path.exists():
                path.unlink()
                event["removed_checkpoint_bytes"] += item["size_bytes"]
            item["cleanup_pending"] = False
            self.save()
        event["total_after_bytes"] = self.total_bytes()
        event["net_reduction_bytes"] = event["total_before_bytes"]-event["total_after_bytes"]
        self.registry.setdefault("cleanup_events",[]).append(event)
        self.save()

    def total_bytes(self):
        return sum(p.stat().st_size for p in self.dataset.rglob("*") if p.is_file()
                   and p.name not in ("storage_report.csv","storage_summary.json"))

    def report(self, output_dir=None):
        destination = Path(output_dir) if output_dir else self.dataset
        destination.mkdir(parents=True, exist_ok=True)
        owned = self.registry["files"]
        rows = []
        total = 0
        for path in sorted(self.dataset.rglob("*")):
            if not path.is_file() or path.name in ("storage_report.csv", "storage_summary.json"): continue
            size = path.stat().st_size
            total += size
            item = owned.get(str(path), {})
            relative = path.relative_to(self.dataset)
            formal = item.get("formal", path.name in ("seurat_raw.rds", "seurat_QC-3.rds", "seurat_analysis.rds"))
            rows.append(dict(file_path=str(path), size_bytes=size, file_type=path.suffix,
                stage=item.get("stage", relative.parts[0]), formal_retained=formal,
                temporary_checkpoint=item.get("temporary", False), safely_cleaned=False,
                lifecycle=item.get("lifecycle", "UNMANAGED_PROTECTED"),
                cleanup_candidate=path.suffix.lower()==".rds" and not formal,
                requires_user_confirmation=not self.active,
                dependencies=json.dumps(item.get("successor", {"status":"UNKNOWN: no dependency manifest found"}))))
        saved = sum(v["size_bytes"] for v in owned.values() if v["lifecycle"]==CLEANED and not Path(v["path"]).exists())
        for item in owned.values():
            if item["lifecycle"]==CLEANED and not Path(item["path"]).exists():
                rows.append(dict(file_path=item["path"],size_bytes=item["size_bytes"],file_type=".rds",
                    stage=item["stage"],formal_retained=False,temporary_checkpoint=True,safely_cleaned=True,
                    lifecycle=CLEANED,cleanup_candidate=False,requires_user_confirmation=False,dependencies=json.dumps(item["successor"])))
        # Read manifests only; no RDS contents are loaded for inventory/dependencies.
        for manifest in self.dataset.rglob("*.json"):
            if manifest.name in ("registry.json", "storage_summary.json"): continue
            try: record = read(manifest)
            except (ValueError, OSError): continue
            binding = record.get("binding", {}) if isinstance(record, dict) else {}
            source = binding.get("input") or (record.get("input") if isinstance(record,dict) else None)
            if source:
                for row in rows:
                    if row["file_path"]==source:
                        row["dependencies"] = json.dumps({"dependent_manifest":str(manifest),"outputs":list(record.get("outputs", {}))})
        events = self.registry.get("cleanup_events",[])
        before = events[-1]["total_before_bytes"] if events else total
        after = events[-1]["total_after_bytes"] if events else total
        for row in rows:
            row.update(total_before_bytes=before,total_after_bytes=after,current_total_bytes=total,saved_bytes=saved)
        fields = list(rows[0]) if rows else ["file_path","size_bytes","file_type","stage","formal_retained",
            "temporary_checkpoint","safely_cleaned","lifecycle","cleanup_candidate","requires_user_confirmation","dependencies",
            "total_before_bytes","total_after_bytes","current_total_bytes","saved_bytes"]
        with (destination / "storage_report.csv").open("w",encoding="utf-8-sig",newline="") as out:
            writer = csv.DictWriter(out,fieldnames=fields); writer.writeheader(); writer.writerows(rows)
        summary = dict(schema_version=1,mode=self.policy["mode"],managed=self.active,
            total_before_bytes=before,total_after_bytes=after,current_total_bytes=total,saved_bytes=saved,
            estimated_reclaimable_bytes=sum(r["size_bytes"] for r in rows if r["cleanup_candidate"]),
            file_count=len(rows), states={s:sum(v["lifecycle"]==s for v in owned.values()) for s in (RETAINED,CLEANED,INCOMPLETE)},
            report=str(destination / "storage_report.csv"),references=self.registry["references"],
            cleanup_events=events,
            accounting="logical file bytes; before/after are observed at last completed cleanup (metadata updates included); current includes later reports/logs; saved is cumulative removed checkpoint bytes; hardlinks may share allocation")
        atomic(destination / "storage_summary.json", summary)
        return summary


def manifest_valid(record):
    for path, checksum in record.get("outputs", {}).items():
        if not Path(path).is_file() or sha256(path)!=checksum: return False
    for value in record.get("cleaned_outputs", {}).values():
        if not Storage.valid_proof(value["successor"]): return False
    return True


@contextmanager
def session(storage):
    """One writer per dataset; OS lock releases even after a native/process crash."""
    keys = ("SINGLE_CELL_STORAGE_MODE", "SINGLE_CELL_QC_CHECKPOINT_ROOT", "SINGLE_CELL_QC_RUN_ROOT")
    previous = {key:os.environ.get(key) for key in keys}
    lock = None
    try:
        if storage.active:
            storage.work.mkdir(parents=True,exist_ok=True)
            lock = (storage.work / "writer.lock").open("a+b")
            if lock.seek(0,2)==0: lock.write(b"0"); lock.flush()
            lock.seek(0)
            if os.name=="nt":
                import msvcrt
                msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
            else:
                import fcntl
                fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
            # Reload after acquiring the lock; never overwrite another writer.
            if storage.registry_path.exists(): storage.registry = read(storage.registry_path)
            os.environ[keys[0]] = storage.policy["mode"]
            os.environ[keys[1]] = storage.checkpoint_root("qc").as_posix()
            os.environ.pop(keys[2],None)
        else:
            for key in keys: os.environ.pop(key,None)
        yield storage
    finally:
        if lock: lock.close()
        for key,value in previous.items():
            if value is None: os.environ.pop(key,None)
            else: os.environ[key] = value


def formal_validate(storage, final, stage, scope, run_r45, repo, log):
    """Caller must first finish domain validation and successful runtime exit."""
    result_path = storage.checkpoint_root(scope) / "storage_validation_result.json"
    item = storage.registry["files"].get(str(Path(final).resolve()))
    try:
        run_r45(Path(repo) / "runtime/validate_formal.R",[final,stage],
                stage="storage_validation",log_path=log,cwd=repo)
    except BaseException as error:
        if item: item["lifecycle"] = INCOMPLETE
        storage.save()
        atomic(result_path,dict(status=INCOMPLETE,stage=stage,path=str(final),error=str(error),cleanup_allowed=False))
        raise
    if item: item["lifecycle"] = RETAINED
    storage.save()
    proof = storage.receipt(final,stage,scope)
    atomic(result_path,dict(status=RETAINED,stage=stage,path=str(final),cleanup_allowed=True,proof=proof))
    return proof


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Read-only inventory of old outputs; never enroll or delete")
    parser.add_argument("dataset",type=Path)
    parser.add_argument("--output-dir",type=Path,help="Write inventory here without changing old dataset")
    parser.add_argument("--collection",action="store_true",help="Inventory immediate dataset directories separately")
    args = parser.parse_args()
    if not args.dataset.is_dir(): parser.error("Dataset/collection directory does not exist")
    if args.collection and not args.output_dir: parser.error("Collection inventory needs --output-dir")
    datasets = sorted(p for p in args.dataset.resolve().iterdir() if p.is_dir() and not p.name.startswith(".")) if args.collection else [args.dataset]
    summaries = []
    for dataset in datasets:
        # Construct without enrollment. Reports are the only writes.
        storage = Storage.__new__(Storage)
        storage.dataset = dataset.resolve()
        storage.work = storage.dataset / ".work/checkpoints"
        storage.registry_path = storage.work / "registry.json"
        storage.registry = read(storage.registry_path) if storage.registry_path.is_file() else {"files":{},"references":{}}
        storage.active = False
        storage.policy = {"mode":"protected_legacy"}
        destination = args.output_dir / dataset.name if args.collection else args.output_dir
        summaries.append(dict(dataset=str(dataset),**storage.report(destination)))
    print(json.dumps(summaries if args.collection else summaries[0],ensure_ascii=False,indent=2))
