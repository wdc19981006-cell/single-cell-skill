"""Collect small diagnostic artifacts; keep RDS and complete cell maps local."""
import argparse
import csv
import json
import shutil
from collections import Counter
from pathlib import Path


def read_csv(path):
    return list(csv.DictReader(path.open(encoding='utf-8-sig')))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workspace',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    latest=args.workspace/'latest/data/GSE155698'
    control=args.workspace/'legacy_pc_control/data/GSE155698'
    raw_rows=read_csv(args.workspace/'history/raw_metadata.csv')
    raw_counts=Counter(r['sample'] for r in raw_rows)
    group={r['sample']:r['group'] for r in raw_rows}
    for mode,folder in [('current',latest),('legacy_pc_control',control)]:
        rows=read_csv(folder/'diagnostics/pca_metrics.csv')
        for row in rows:
            row['raw_cells']=sum(raw_counts.values()) if row['sample']=='ALL' else raw_counts[row['sample']]
            row['group']='ALL' if row['sample']=='ALL' else group[row['sample']]
        with (args.output/(mode+'_pca_diagnostics.csv')).open('w',encoding='utf-8',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
        shutil.copyfile(folder/'diagnostics/diagnostic_run.json',args.output/(mode+'_diagnostic_run.json'))
        provenance=json.loads((folder/'validation_provenance.json').read_text())
        assert provenance['historical_files_unchanged']
        shutil.copyfile(folder/'validation_provenance.json',args.output/(mode+'_provenance.json'))
    # All 21 current elbows and spectra are small standalone diagnostics.
    for pattern in ('elbow_*.png','*_pca_spectrum.csv'):
        for source in (latest/'diagnostics').glob(pattern):
            target=args.output/('elbows' if source.suffix=='.png' else 'spectra')/source.name
            target.parent.mkdir(exist_ok=True)
            shutil.copyfile(source,target)
    for source in (args.workspace/'comparison').iterdir():
        if source.suffix in ('.json','.csv','.png'):
            shutil.copyfile(source,args.output/source.name)
    shutil.copyfile(args.workspace/'history/history_inventory.json',args.output/'history_inventory.json')
    print('Packaged diagnostics; no expression matrices or RDS copied')


if __name__=='__main__':
    main()
