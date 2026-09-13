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
from common import ROOT, local, read_csv, validate_manifest, verify_confirmation

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''): h.update(block)
    return h.hexdigest()

def download(url, target, expected=None):
    target.parent.mkdir(parents=True, exist_ok=True)
    stamp = target.with_name(target.name + '.download.json')
    if target.exists():
        if not stamp.exists(): raise ValueError('Unverified preexisting file: ' + str(target))
        old = json.loads(stamp.read_text())
        if old['url'] != url or old['sha256'] != sha256(target): raise ValueError('Existing download does not match provenance')
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
    stamp.write_text(json.dumps(record, indent=2), encoding='utf-8')
    return record

def extract_member(archive, member, destination):
    """Never extract archive paths directly; stream one explicitly mapped regular file."""
    if member.startswith(('/', '\\')) or '..' in Path(member.replace('\\', '/')).parts: raise ValueError('Unsafe archive member')
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
    plans = {}
    for row in rows:
        for item in json.loads(row.get('files_json') or '[]'): plans[item['local_path']] = item
    records, archives = [], {}
    for item in plans.values():
        target = local(root, item['local_path'], 'datasets/' + rows[0]['database'])
        if item.get('member'):
            cache = root / 'datasets' / rows[0]['database'] / '_archives' / (hashlib.sha256(item['url'].encode()).hexdigest() + '.archive')
            if item['url'] not in archives: archives[item['url']] = download(item['url'], cache, item.get('sha256'))
            extract_member(cache, item['member'], target)
            records.append(dict(**archives[item['url']], local_path=item['local_path'], member=item['member'], extracted_sha256=sha256(target)))
        else: records.append(dict(**download(item['url'], target, item.get('sha256')), local_path=item['local_path']))
    log = root / 'logs' / (rows[0]['database'] + '_download.json')
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(json.dumps(records, indent=2), encoding='utf-8')
    print(str(log))

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('manifest', type=Path); p.add_argument('--root', type=Path, default=ROOT)
    a = p.parse_args(); run(a.manifest, a.root.resolve())
