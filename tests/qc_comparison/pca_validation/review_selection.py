"""Independently recompute current PC choices from the exported spectrum."""
import argparse
import csv
import json
import math
from pathlib import Path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('diagnostics',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    rows=list(csv.DictReader((args.diagnostics/'pca_metrics.csv').open(encoding='utf-8-sig')))
    assert 'ALL' in [r['sample'] for r in rows]
    assert len({r['sample'] for r in rows})==len(rows)
    assert all(r['mode']=='latest' for r in rows)
    checks=[]
    for row in rows:
        spectrum=list(csv.DictReader((args.diagnostics/(row['sample']+'_pca_spectrum.csv')).open()))
        sd=[float(r['stdev']) for r in spectrum]
        pct=[v/sum(sd)*100 for v in sd]
        variance=[v*v/sum(x*x for x in sd)*100 for v in sd]
        cumulative=0;variance_pick=None
        for i,v in enumerate(pct,1):
            cumulative+=v
            if variance_pick is None and cumulative>90 and v<5:
                variance_pick=i
        drops=[i+2 for i in range(len(pct)-1) if pct[i]-pct[i+1]>0.1]
        drop_pick=max(drops) if drops else None
        candidates=[x for x in (variance_pick,drop_pick) if x is not None]
        max_pc=min(50,int(row['cells'])-1,int(row['genes'])-1,int(row['HVGs'])-1)
        selected=max(2,min(min(candidates) if candidates else 20,max_pc))
        assert int(row['max_pc'])==max_pc and int(row['pc_used'])==selected
        recorded=lambda name: None if row[name] in ('','NA') else int(row[name])
        assert recorded('variance_pick')==variance_pick and recorded('drop_pick')==drop_pick
        assert all(math.isclose(float(r['sd_share']),pct[i],abs_tol=1e-9) for i,r in enumerate(spectrum))
        assert all(math.isclose(float(r['variance_explained']),variance[i],abs_tol=1e-9) for i,r in enumerate(spectrum))
        assert math.isclose(float(row['cumulative_sd_share']),sum(pct[:selected]),abs_tol=1e-9)
        assert math.isclose(float(row['cumulative_variance_explained']),sum(variance[:selected]),abs_tol=1e-9)
        checks.append(dict(sample=row['sample'],pc_used=selected,max_pc=max_pc,
                           variance_pick=variance_pick,drop_pick=drop_pick,status='PASS'))
    args.output.write_text(json.dumps(dict(status='PASS',independent_python_recalculation=True,
        variance_denominator='sum squared SD across computed PCs, not total transcriptome variance',
        checks=checks),indent=2),encoding='utf-8')
    print('PASS independent spectrum review:',len(checks),'PCA selections')


if __name__=='__main__':
    main()
