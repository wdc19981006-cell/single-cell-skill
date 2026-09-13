"""Generate tiny format fixtures in ignored datasets; never store binary inputs in Git."""
import csv
import gzip
import json
import sys
from pathlib import Path
import numpy as np
import h5py
import anndata
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'.agents/skills/geo-single-cell-loader/scripts'))
from common import write_csv
from build_manifest import confirm

def encoding(obj, kind, version='0.2.0'):
    obj.attrs['encoding-type']=kind; obj.attrs['encoding-version']=version

def strings(group,name,values):
    d=group.create_dataset(name,data=np.array(values,dtype=object),dtype=h5py.string_dtype('utf-8')); encoding(d,'string-array'); return d

def main():
    gse='GSE999999999'; base=ROOT/'datasets'/gse; base.mkdir(parents=True,exist_ok=True)
    genes=['Gene'+str(i) for i in range(250)]; cells=['AAAC'+str(i)+'-1' for i in range(4)]
    counts=np.ones((250,4),dtype=np.int32); counts[0,:]=[2,3,4,5]
    trio=base/'trio'; trio.mkdir(exist_ok=True)
    with gzip.open(trio/'matrix.mtx.gz','wt') as f:
        f.write('%%MatrixMarket matrix coordinate integer general\n250 4 1000\n')
        for c in range(4):
            for r in range(250): f.write(f'{r+1} {c+1} {counts[r,c]}\n')
    with gzip.open(trio/'features.tsv.gz','wt') as f: f.writelines(f'ENSG{i}\t{g}\tGene Expression\n' for i,g in enumerate(genes))
    with gzip.open(trio/'barcodes.tsv.gz','wt') as f: f.write('\n'.join(cells)+'\n')
    with h5py.File(base/'filtered_feature_bc_matrix.h5','w') as f:
        m=f.create_group('matrix'); m['data']=counts.T.ravel(); m['indices']=np.tile(np.arange(250,dtype=np.int32),4); m['indptr']=np.arange(0,1001,250,dtype=np.int32); m['shape']=np.array([250,4],dtype=np.int32)
        m.create_dataset('barcodes',data=np.array(cells,dtype='S'))
        features=m.create_group('features')
        for name,values in [('id',['ENSG'+str(i) for i in range(250)]),('name',genes),('feature_type',['Gene Expression']*250),('genome',['GRCh38']*250)]: features.create_dataset(name,data=np.array(values,dtype='S'))
    obs=pd.DataFrame(index=cells); var=pd.DataFrame(index=genes)
    adata=anndata.AnnData(X=np.log1p(counts.T.astype(float)),obs=obs,var=var)
    adata.layers['counts']=counts.T
    adata.raw=anndata.AnnData(X=counts.T,obs=obs.copy(),var=var.copy())
    adata.write_h5ad(base/'counts.h5ad')
    with (base/'counts.h5ad').open('rb') as src,gzip.open(base/'counts.h5ad.gz','wb') as dst: dst.write(src.read())
    with (base/'counts.csv').open('w',newline='') as f:
        w=csv.writer(f); w.writerow(['gene','description']+cells)
        for gene,values in zip(genes,counts): w.writerow([gene,'synthetic']+list(values))
    with (base/'transposed.tsv').open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(['cell']+genes)
        for cell,values in zip(cells,counts.T): w.writerow([cell]+list(values))
    rows=[]
    for sample,kind,path,source in [('FixtureZ','10x_mtx','trio','counts'),('FixtureA','10x_h5','filtered_feature_bc_matrix.h5','counts'),('FixtureH','h5ad','counts.h5ad.gz','counts'),('FixtureT','text','counts.csv','counts')]:
        rows.append(dict(database=gse,sample=sample,tissue='Synthetic',disease='Synthetic',source_type='Tissue',local_path=f'datasets/{gse}/{path}',file_type=kind,count_source=source,count_evidence='Generated integer fixture counts; H5AD X intentionally normalized',metadata_evidence='Synthetic test only',files_json='[]',delimiter='comma',feature_column='gene',drop_columns='description',orientation='genes_by_cells'))
    folder=ROOT/'manifests'; folder.mkdir(exist_ok=True)
    report=folder/'synthetic_report.csv'; write_csv(report,rows)
    approval=folder/'synthetic_approval.json'; approval.write_text(json.dumps(dict(confirmed_by='user',user_statement='SIMULATED FIXTURE APPROVAL ONLY; no real biological data',groups={r['sample']:('Fixture_A' if r['sample']=='FixtureA' else 'Fixture_B') for r in reversed(rows)})))
    output=folder/(gse+'_sample_manifest.csv')
    if output.exists(): output.unlink() # Only this generated, reserved fixture manifest.
    confirm(report,approval,output,ROOT)
    print('Generated synthetic fixtures and mock manifest:',gse)

if __name__=='__main__': main()
