"""Download only manifest-authorized files after grouping confirmation."""
import argparse
import hashlib
import json
import shutil
import tarfile
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from common import ROOT, dataset_paths, local, read_csv, validate_manifest, verify_confirmation, workflow_path

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''): h.update(block)
    return h.hexdigest()

def download(url, target, expected=None, previous=None):
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        if not previous: raise ValueError('Unverified preexisting file: ' + str(target))
        old = previous
        if old['url'] != url or old['sha256'] != sha256(target) or old['bytes'] != target.stat().st_size: raise ValueError('Existing download does not match provenance')
        if expected and old['sha256'] != expected: raise ValueError('Checksum mismatch')
        return old
    part = target.with_name(target.name + '.part')
    request = urllib.request.Request(url, headers={'User-Agent': 'geo-single-cell-loader/1.0'})
    with urllib.request.urlopen(request, timeout=120) as response, part.open('wb') as out:
        shutil.copyfileobj(response, out, length=1024 * 1024)
    actual = sha256(part)
    if expected and actual != expected: raise ValueError('Checksum mismatch; partial retained for diagnosis')
    part.replace(target)
    record = dict(url=url, sha256=actual, downloaded_at=datetime.now(timezone.utc).isoformat(), bytes=target.stat().st_size)
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
                if not old or old.get('url') != item['url'] or old.get('member') != item['member'] or old.get('extracted_sha256') != sha256(target) or old.get('extracted_bytes') != target.stat().st_size:
                    raise ValueError('Existing archive member does not match provenance')
                if item.get('sha256') and old.get('sha256') != item['sha256']: raise ValueError('Checksum mismatch')
                continue
            cache_rel = f'data/{gse}/raw/_archives/' + hashlib.sha256(item['url'].encode()).hexdigest() + '.archive'
            cache = local(root, cache_rel, f'data/{gse}/raw')
            archive_record = download(item['url'], cache, item.get('sha256'), by_path.get(cache_rel))
            record(dict(archive_record, local_path=cache_rel))
            extract_member(cache, item['member'], target)
            record(dict(archive_record, local_path=item['local_path'], member=item['member'], extracted_sha256=sha256(target), extracted_bytes=target.stat().st_size))
        else:
            if old and old.get('member'): raise ValueError('Download mapping differs from existing provenance')
            record(dict(download(item['url'], target, item.get('sha256'), old), local_path=item['local_path']))
    if not log.exists():
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text('[]\n', encoding='utf-8')
    print(str(log))

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('manifest', type=Path); p.add_argument('--root', type=Path, default=ROOT)
    a = p.parse_args(); run(a.manifest, a.root.resolve())
