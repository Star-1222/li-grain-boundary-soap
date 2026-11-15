"""Supplemental diagnostics; preserve the agreed two main distance definitions."""
import run_sweep as core
import json
import itertools
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

ROOT=core.ROOT
def main():
    pairs=dict(np.load(ROOT/'fixed_pairs.npz'));candidates=json.loads((ROOT/'candidates.json').read_text(encoding='utf-8'))
    categories=[]
    for c in candidates:
        hi=c['config_id'];lo=core.config_id(c['r_cut'],c['sigma'],10,c['l_max'])
        a=np.load(ROOT/'metrics'/f'{lo}.npz');b=np.load(ROOT/'metrics'/f'{hi}.npz')
        for metric in ['raw','norm']:
            for ci,name in enumerate(pairs['names']):
                mask=pairs['category']==ci;x=a['pair_'+metric][mask];y=b['pair_'+metric][mask]
                categories.append(dict(lower=lo,upper=hi,metric=metric,category=name,n=int(mask.sum()),
                    spearman=core.rho(x,y),median_relative_change=float(np.median(np.abs(x-y))/max(np.median(y),1e-12))))
    pd.DataFrame(categories).to_csv(ROOT/'radial_category_diagnostics.csv',index=False)
    # DScribe single-species ordering verified in v2.1.2 soapGTO.cpp getPD:
    # l outermost, then n, then n' >= n. Weight sqrt(2) on off-diagonal entries.
    # Diagnostic only: no change to the agreed raw/norm results or selections.
    rows=[]
    for rc,sigma in [(3.25,.5),(4,.3),(5,.3),(6,.3)]:
        values={}
        for n in [10,12]:
            cid=core.config_id(rc,sigma,n,10);x=np.load(ROOT/'vectors'/f'{cid}.npz')['evaluation']
            w=np.tile([1. if i==j else np.sqrt(2.) for i in range(n) for j in range(i,n)],11)
            for weighting in ['unique_terms','symmetric_frobenius']:
                z=x if weighting=='unique_terms' else x*w
                for metric in ['raw','norm']:
                    xx=z if metric=='raw' else core.unit(z)
                    d=cdist(xx,xx);np.fill_diagonal(d,np.inf)
                    values[n,weighting,metric]=(d[pairs['i'],pairs['j']],np.argpartition(d,9,axis=1)[:,:10])
        for weighting,metric in itertools.product(['unique_terms','symmetric_frobenius'],['raw','norm']):
            a,nn1=values[10,weighting,metric];b,nn2=values[12,weighting,metric]
            corrs=[];changes=[]
            for ci in range(10):
                mask=pairs['category']==ci
                corrs.append(core.rho(a[mask],b[mask]));changes.append(np.median(np.abs(a[mask]-b[mask]))/np.median(b[mask]))
            row=dict(r_cut=rc,sigma=sigma,l_max=10,lower_n=10,upper_n=12,
                weighting=weighting,metric=metric,spearman_min=min(corrs),relative_change_max=max(changes),
                nn_overlap_mean=float(core.overlap(nn1,nn2).mean()))
            row['passed']=row['spearman_min']>=.98 and row['relative_change_max']<=.05 and row['nn_overlap_mean']>=.8
            rows.append(row)
    pd.DataFrame(rows).to_csv(ROOT/'symmetry_weight_diagnostic.csv',index=False)
    # Independent algebraic check of the symmetric feature weighting convention.
    rng=np.random.default_rng(core.SEED+500);a=rng.normal(size=(8,8));a=(a+a.T)/2
    q,_=np.linalg.qr(rng.normal(size=(8,8)));b=q@a@q.T
    iu=np.triu_indices(8);weights=np.where(iu[0]==iu[1],1.,2.)
    np.testing.assert_allclose(np.sum(weights*a[iu]**2),np.linalg.norm(a,'fro')**2,rtol=1e-13)
    np.testing.assert_allclose(np.sum(weights*a[iu]**2),np.sum(weights*b[iu]**2),rtol=1e-13)
    core.save_json(ROOT/'metric_diagnostics.json',dict(main_distances_unchanged=True,
        symmetry_source='https://github.com/SINGROUP/dscribe/blob/v2.1.2/dscribe/ext/soapGTO.cpp#L2240',
        convention='single species; l-major unique n<=n_prime; multiply n!=n_prime vector components by sqrt(2)',
        role='supplemental comparison only, not a replacement of raw/norm or a full alternate sweep',
        radial_category_rows=len(categories),symmetry_comparisons=len(rows),
        symmetry_edges_passing=sum(r['passed'] for r in rows),weighting_algebra_verified=True))
    core.log('metric_diagnostics_complete',radial_category_rows=len(categories),symmetry_edges_passing=sum(r['passed'] for r in rows))

if __name__=='__main__':main()
