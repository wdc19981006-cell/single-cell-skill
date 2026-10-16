"""Resolve historical RDS from a local mirror only after Git LFS SHA verification."""
import csv
import hashlib
import shutil
from pathlib import Path
import argparse


def sha256(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference', type=Path)
    parser.add_argument('local', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('cache', type=Path)
    args = parser.parse_args()
    rows = []
    manifest = list(csv.DictReader((args.reference / 'file_manifest.csv').open(encoding='utf-8-sig')))
    for entry in manifest:
        if not entry['path'].endswith('.rds'):
            continue
        pointer = (args.reference / entry['path']).read_text()
        assert 'version https://git-lfs.github.com/spec/v1' in pointer
        assert 'oid sha256:' + entry['sha256'] in pointer
        assert 'size ' + entry['bytes'] in pointer
        gse, relative = entry['path'].split('/', 1)
        candidates = [folder / relative for folder in args.local.glob(gse + '*')]
        matches = [p for p in candidates if p.is_file() and p.stat().st_size == int(entry['bytes'])
                   and sha256(p) == entry['sha256']]
        if len(matches) != 1:
            raise RuntimeError(f'Expected exactly one verified real RDS: {entry["path"]}')
        cached = args.cache / entry['path']
        cached.parent.mkdir(parents=True, exist_ok=True)
        if cached.exists():
            assert sha256(cached) == entry['sha256']
        else:
            shutil.copyfile(matches[0], cached)
            assert sha256(cached) == entry['sha256']
        rows.append(dict(dataset=gse, version='original' if '/qc/' not in entry['path'] else 'old_skill',
                         path=cached.as_posix(), bytes=entry['bytes'], sha256=entry['sha256'],
                         lfs_verified=True))
        print('VERIFIED', entry['path'], entry['sha256'], flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for entry in manifest:
        if entry['path'].endswith('/qc/qc_report.csv'):
            source = args.reference / entry['path']
            assert sha256(source) == entry['sha256']
            destination = args.cache / entry['path']
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, destination)
    with args.output.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()
