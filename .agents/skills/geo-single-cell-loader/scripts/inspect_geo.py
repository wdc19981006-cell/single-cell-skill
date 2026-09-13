"""Metadata-only GEO discovery. Prefer bio-server MCP when available."""
import argparse
import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from common import ROOT, choose_filtered, file_type, write_csv

def propose_input(gse, sample, urls):
    """Use GEO's sample ownership, never infer biological IDs from filenames."""
    base = 'datasets/' + gse + '/' + sample
    trio = {}
    for role in ('matrix.mtx', 'features.tsv', 'barcodes.tsv'):
        matches = [u for u in urls if u.endswith((role, role + '.gz'))]
        if len(matches) == 1: trio[role] = matches[0]
    if len(trio) == 3:
        files = [dict(url=u, local_path=base + '/' + role + ('.gz' if u.endswith('.gz') else '')) for role,u in trio.items()]
        return dict(local_path=base, file_type='10x_mtx', count_source='counts', files_json=json.dumps(files))
    h5 = [u for u in urls if file_type(u) == '10x_h5_candidate']
    candidates = h5 or [u for u in urls if file_type(u) in ('h5ad', 'rds')]
    if candidates:
        try: selected = choose_filtered(candidates)
        except ValueError: return dict(status='ambiguous', candidates=candidates)
        kind = file_type(selected).replace('_candidate','')
        target = base + '/' + urllib.parse.unquote(urllib.parse.urlparse(selected).path.rsplit('/',1)[1])
        return dict(local_path=target, file_type=kind, count_source='counts' if kind=='10x_h5' else '', files_json=json.dumps([dict(url=selected,local_path=target)]), candidates=candidates)
    return dict(status='requires archive/text/pooled mapping inspection', candidates=urls)

def fetch(url, limit=32 * 1024 * 1024):
    for attempt in range(3):
        try:
            request = urllib.request.Request(url, headers={'User-Agent': 'geo-single-cell-loader/1.0 metadata-inspection'})
            with urllib.request.urlopen(request, timeout=60) as response:
                data = response.read(limit + 1)
                if len(data) > limit: raise ValueError('Metadata response exceeds safety limit')
                return data.decode('utf-8')
        except (OSError, TimeoutError):
            if attempt == 2: raise
            time.sleep(2 ** attempt)

def soft(accession, target='self'):
    url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?' + urllib.parse.urlencode(dict(acc=accession, targ=target, form='text', view='full'))
    text = fetch(url)
    if '^' not in text or '!Series_' not in text and '!Sample_' not in text: raise ValueError('GEO returned no SOFT metadata (possibly blocked)')
    return text, url

def parse_soft(text):
    records, current = [], None
    for line in text.splitlines():
        if line.startswith('^') and ' = ' in line:
            kind, accession = line[1:].split(' = ', 1)
            current = {'kind': kind, 'accession': accession, 'fields': {}}
            records.append(current)
        elif current and line.startswith('!') and ' = ' in line:
            key, value = line[1:].split(' = ', 1)
            current['fields'].setdefault(key, []).append(value)
    return records

def supplementary(fields):
    return [u.replace('ftp://ftp.ncbi.nlm.nih.gov/', 'https://ftp.ncbi.nlm.nih.gov/') for k, vals in fields.items() if 'supplementary_file' in k for u in vals if u.upper() != 'NONE']

def inspect(gse, root, max_samples=None):
    if not re.fullmatch(r'GSE\d+', gse): raise ValueError('Expected GSE accession')
    text, url = soft(gse)
    series = next(r for r in parse_soft(text) if r['kind'] == 'SERIES')
    fields = series['fields']
    samples = fields.get('Series_sample_id', [])
    chosen = samples if max_samples is None else samples[:max_samples]
    records, report, files, input_plans = [], [], [], {}
    for u in supplementary(fields): files.append(dict(scope=gse, url=u, format=file_type(u)))
    for gsm in chosen:
        time.sleep(0.4)
        raw, source = soft(gsm)
        record = next(r for r in parse_soft(raw) if r['kind'] == 'SAMPLE')
        records.append(record)
        f = record['fields']
        characteristics = f.get('Sample_characteristics_ch1', [])
        # Preserve author facts verbatim; semantic standardization requires evidence review.
        report.append(dict(database=gse, sample=gsm, author_sample='; '.join(f.get('Sample_title', [])), tissue='', disease='', source_type='', sample_description='; '.join(f.get('Sample_source_name_ch1', []) + characteristics), metadata_evidence=source))
        for u in supplementary(f): files.append(dict(scope=gsm, url=u, format=file_type(u)))
        input_plans[gsm] = propose_input(gse,gsm,supplementary(f))
    papers = []
    for pmid in fields.get('Series_pubmed_id', []):
        paper_url = 'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pubmed.cgi/BioC_xml/' + pmid + '/unicode'
        try: papers.append(dict(pmid=pmid, url=paper_url, text=fetch(paper_url, 4 * 1024 * 1024)))
        except (OSError, ValueError) as e: papers.append(dict(pmid=pmid, url=paper_url, error=str(e)))
    folder = root / 'manifests'; folder.mkdir(parents=True, exist_ok=True)
    result = dict(gse=gse, series_url=url, series=series, samples=records, files=files, input_plans=input_plans, papers=papers, total_samples=len(samples), inspected_samples=len(records), complete=len(records) == len(samples), note='group 尚未创建，请确认分组方式。Blank biological fields require evidence review, never infer group.')
    (folder / (gse + '_inspection.json')).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    write_csv(folder / (gse + '_sample_report.csv'), report, ['database','sample','author_sample','tissue','disease','source_type','sample_description','metadata_evidence'])
    print(json.dumps(dict(gse=gse, title=fields.get('Series_title'), total_samples=len(samples), inspected_samples=len(records), files=files, note=result['note']), ensure_ascii=False))
    return result

def search(query):
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?' + urllib.parse.urlencode(dict(db='gds', term=query + ' AND GSE[ETYP]', retmode='json', retmax=20))
    return fetch(url)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('gse', nargs='?'); p.add_argument('--search'); p.add_argument('--root', type=Path, default=ROOT)
    p.add_argument('--max-samples', type=int, help='Development discovery only; never treat this as a full sample report')
    a = p.parse_args()
    if a.search: print(search(a.search))
    elif a.gse: inspect(a.gse.upper(), a.root.resolve(), a.max_samples)
    else: p.error('Provide GSE or --search')

if __name__ == '__main__': main()
