"""Download only manifest-authorized files after grouping confirmation."""
import argparse
import hashlib
import json
import shutil
import tarfile
import urllib.request
import zipfile
import csv
import time
import os
from datetime import datetime, timezone
from pathlib import Path
from common import ROOT, dataset_paths, local, read_csv, validate_manifest, verify_confirmation, workflow_path

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''): h.update(block)
    return h.hexdigest()

def profile(log, **fields):
    path = log.parent / 'download_profile.csv'
    columns = ['local_path','url','member','bytes','sha256','downloaded_at','extracted_sha256','extracted_bytes','status','download_seconds','mb_per_second','retry','sha_seconds','extract_seconds']
    exists = path.exists()
    with path.open('a', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        if not exists: writer.writeheader()
        writer.writerow({key: fields.get(key, '') for key in columns})

def download(url, target, expected=None, previous=None):
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        if not previous: raise ValueError('Unverified preexisting file: ' + str(target))
        old = previous
        sha_started = time.perf_counter()
        actual = sha256(target)
        sha_seconds = time.perf_counter() - sha_started
        if old['url'] != url or old['sha256'] != actual or int(old['bytes']) != target.stat().st_size: raise ValueError('Existing download does not match provenance')
        if expected and old['sha256'] != expected: raise ValueError('Checksum mismatch')
        return dict(old, status='REUSED', download_seconds=0, mb_per_second=0, retry=0, sha_seconds=sha_seconds, extract_seconds=0)
    part = target.with_name(target.name + '.part')
    request = urllib.request.Request(url, headers={'User-Agent': 'geo-single-cell-loader/1.0'})
    started = time.perf_counter()
    with urllib.request.urlopen(request, timeout=120) as response, part.open('wb') as out:
        shutil.copyfileobj(response, out, length=1024 * 1024)
    download_seconds = time.perf_counter() - started
    sha_started = time.perf_counter()
    actual = sha256(part)
    sha_seconds = time.perf_counter() - sha_started
    if expected and actual != expected: raise ValueError('Checksum mismatch; partial retained for diagnosis')
    part.replace(target)
    size = target.stat().st_size
    record = dict(url=url, sha256=actual, downloaded_at=datetime.now(timezone.utc).isoformat(), bytes=size,
                  status='DOWNLOADED',download_seconds=download_seconds,mb_per_second=size/1000000/download_seconds if download_seconds else 0,
                  retry=0,sha_seconds=sha_seconds,extract_seconds=0)
    return record

def extract_member(archive, member, destination):
    """Never extract archive paths directly; stream one explicitly mapped regular file."""
    if member.startswith(('/', '\\')) or '..' in Path(member.replace('\\', '/')).parts: raise ValueError('Unsafe archive member')
    if destination.exists(): raise ValueError('Refusing to overwrite existing archive member')
    destination.parent.mkdir(parents=True, exist_ok=True)
    part = destination.with_name(destination.name + '.part')
    if zipfile.is_zipfile(archive):
        with zipfile.ZipFile(archive) as z:
            matches = [i for i in z.infolist() if i.filename == member]
            if len(matches) != 1 or matches[0].is_dir() or (matches[0].external_attr >> 16) & 0o170000 == 0o120000: raise ValueError('Ambiguous/nonregular archive member')
            with z.open(matches[0]) as src, part.open('wb') as dst: shutil.copyfileobj(src, dst)
    else:
        with tarfile.open(archive) as t:
            matches = [i for i in t.getmembers() if i.name == member]
            if len(matches) != 1 or not matches[0].isfile(): raise ValueError('Ambiguous/nonregular archive member')
            with t.extractfile(matches[0]) as src, part.open('wb') as dst: shutil.copyfileobj(src, dst)
    part.replace(destination)

def run(manifest, root):
    rows = validate_manifest(read_csv(manifest), root); verify_confirmation(manifest)
    gse = rows[0]['database']
    workflow_path(root, manifest, gse)
    paths = dataset_paths(root, gse)
    log = paths['workflow'] / 'download.json'
    records = json.loads(log.read_text(encoding='utf-8')) if log.exists() else []
    by_path = {r['local_path']: r for r in records}
    if len(by_path) != len(records): raise ValueError('Duplicate download provenance records')
    def record(item):
        profile(log, **item)
        by_path[item['local_path']] = item
        log.parent.mkdir(parents=True, exist_ok=True)
        pending = log.with_suffix('.pending.json')
        pending.write_text(json.dumps(list(by_path.values()), indent=2), encoding='utf-8')
        pending.replace(log)
    plans = {}
    for row in rows:
        for item in json.loads(row.get('files_json') or '[]'): plans[item['local_path']] = item
    for item in plans.values():
        target = local(root, item['local_path'], f'data/{gse}/raw')
        old = by_path.get(item['local_path'])
        if item.get('member'):
            if target.exists():
                sha_started = time.perf_counter()
                existing_sha = sha256(target)
                sha_seconds = time.perf_counter() - sha_started
                if not old or old.get('url') != item['url'] or old.get('member') != item['member'] or old.get('extracted_sha256') != existing_sha or int(old.get('extracted_bytes') or -1) != target.stat().st_size:
                    raise ValueError('Existing archive member does not match provenance')
                if item.get('sha256') and old.get('sha256') != item['sha256']: raise ValueError('Checksum mismatch')
                profile(log, **dict(old, status='REUSED',download_seconds=0,mb_per_second=0,retry=0,sha_seconds=sha_seconds,extract_seconds=0))
                continue
            cache_rel = f'data/{gse}/raw/_archives/' + hashlib.sha256(item['url'].encode()).hexdigest() + '.archive'
            cache = local(root, cache_rel, f'data/{gse}/raw')
            archive_record = download(item['url'], cache, item.get('sha256'), by_path.get(cache_rel))
            record(dict(archive_record, local_path=cache_rel))
            extract_started = time.perf_counter()
            extract_member(cache, item['member'], target)
            extract_seconds = time.perf_counter() - extract_started
            sha_started = time.perf_counter()
            extracted_sha = sha256(target)
            record(dict(archive_record, local_path=item['local_path'], member=item['member'], extracted_sha256=extracted_sha, extracted_bytes=target.stat().st_size,
                        status='DOWNLOADED',extract_seconds=extract_seconds,sha_seconds=archive_record['sha_seconds']+time.perf_counter()-sha_started))
        else:
            if old and old.get('member'): raise ValueError('Download mapping differs from existing provenance')
            record(dict(download(item['url'], target, item.get('sha256'), old), local_path=item['local_path']))
    if not log.exists():
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text('[]\n', encoding='utf-8')
    print(str(log))

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('manifest', type=Path); p.add_argument('--root', type=Path, default=ROOT)
    a = p.parse_args()
    if os.environ.get('GEO_AUDIT_RUN_ACTIVE') != '1':
        p.error('Real downloads must use run_confirmed.py so every terminal run is audited')
    run(a.manifest, a.root.resolve())
