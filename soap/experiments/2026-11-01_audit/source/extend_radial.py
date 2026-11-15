"""Pre-registered follow-up: check upper radial boundary without altering main-grid criteria."""
import run_sweep as core
import json
import warnings
import traceback
import time
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

ROOT=core.ROOT
def main():
    output=ROOT/'radial_extension';output.mkdir(exist_ok=True)
    candidates=json.loads((ROOT/'candidates.json').read_text(encoding='utf-8'))
    protocol=dict(registered_at=core.stamp(),reason='No radial edge in n=6/8/10/12 passes original criteria.',
        n_max=[14,16],configurations='same r_cut, sigma, l_max as six pilot-selected audit candidates',
        frozen_main_grid=True,criteria='same 0.98 / 5% / 80% edge criteria; no relaxation',
        limitations='Larger bases may fail orthogonalization; success/finite output alone does not establish reliability.')
    if not (output/'protocol.json').exists():core.save_json(output/'protocol.json',protocol)
    atoms=core.structures();frame=core.population()
    pilot=core.ordered_records(frame[frame.pilot]);refs=core.ordered_records(frame[frame.role=='reference'])
    pairs=dict(np.load(ROOT/'fixed_pairs.npz'));rows=[];edges=[]
    for source in candidates:
        previous=source['config_id'];prev_metrics=dict(np.load(ROOT/'metrics'/f'{previous}.npz'))
        for n in [14,16]:
            c={k:source[k] for k in ['r_cut','sigma','n_max','l_max']}
            c['n_max']=n;c['config_id']=core.config_id(c['r_cut'],c['sigma'],n,c['l_max'])
            record=output/f"{c['config_id']}.json"
            if record.exists():
                row=json.loads(record.read_text(encoding='utf-8'))
                if row['status']=='ok':
                    prev_metrics=dict(np.load(output/f"{c['config_id']}_metrics.npz"));previous=c['config_id']
                else:prev_metrics=None
                rows.append(row);continue
            begin=time.perf_counter();row=dict(**c,base_candidate=source['config_id'],registered_at=protocol['registered_at'])
            caught=[]
            try:
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always');soap=core.make_soap(c)
                    row.update(core.basis_health(soap))
                    x=core.calc_vectors(soap,atoms,pilot);ref=core.calc_vectors(soap,atoms,refs)
                m,_=core.metrics(x,ref,pilot,pairs)
                np.savez(output/f"{c['config_id']}_vectors.npz",evaluation=x,reference=ref)
                np.savez(output/f"{c['config_id']}_metrics.npz",**m)
                row.update(status='ok',finite=True)
                if prev_metrics is not None:
                    for metric in ['raw','norm']:
                        a=prev_metrics['pair_'+metric];b=m['pair_'+metric]
                        correlations=[];changes=[]
                        for ci in range(len(pairs['names'])):
                            mask=pairs['category']==ci
                            correlations.append(core.rho(a[mask],b[mask]))
                            changes.append(float(np.median(np.abs(a[mask]-b[mask]))/max(np.median(b[mask]),1e-12)))
                        near=core.overlap(prev_metrics['nn_'+metric],m['nn_'+metric]).mean()
                        e=dict(lower=previous,upper=c['config_id'],metric=metric,spearman_min=min(correlations),
                            relative_change_max=max(changes),nn_overlap_mean=near)
                        e['passed']=e['spearman_min']>=.98 and e['relative_change_max']<=.05 and near>=.8
                        edges.append(e)
                previous=c['config_id'];prev_metrics=m
            except Exception as exc:
                row.update(status='failed',error=repr(exc),traceback=traceback.format_exc());prev_metrics=None
            row['warnings']=sorted(set(str(w.message) for w in caught));row['seconds']=time.perf_counter()-begin
            core.save_json(record,row);rows.append(row)
            core.log('radial_extension_progress',config_id=c['config_id'],status=row['status'],warnings=len(row['warnings']))
    pd.DataFrame([{k:v for k,v in r.items() if not isinstance(v,(list,dict))}|{'warning_count':len(r.get('warnings',[]))} for r in rows]).to_csv(output/'summary.csv',index=False)
    # Recompute edges on a resumed run from durable metrics as well.
    edges=[]
    for source in candidates:
        ids=[source['config_id']]+[core.config_id(source['r_cut'],source['sigma'],n,source['l_max']) for n in [14,16]]
        for lo,hi in zip(ids[:-1],ids[1:]):
            low=(ROOT/'metrics'/f'{lo}.npz') if lo==source['config_id'] else (output/f'{lo}_metrics.npz')
            high=output/f'{hi}_metrics.npz'
            if not low.exists() or not high.exists():continue
            a=np.load(low);b=np.load(high)
            for metric in ['raw','norm']:
                correlations=[];changes=[]
                for ci in range(len(pairs['names'])):
                    mask=pairs['category']==ci;x=a['pair_'+metric][mask];y=b['pair_'+metric][mask]
                    correlations.append(core.rho(x,y));changes.append(float(np.median(np.abs(x-y))/max(np.median(y),1e-12)))
                near=core.overlap(a['nn_'+metric],b['nn_'+metric]).mean()
                e=dict(lower=lo,upper=hi,metric=metric,spearman_min=min(correlations),relative_change_max=max(changes),nn_overlap_mean=near)
                e['passed']=e['spearman_min']>=.98 and e['relative_change_max']<=.05 and near>=.8
                edges.append(e)
    pd.DataFrame(edges,columns=['lower','upper','metric','spearman_min','relative_change_max','nn_overlap_mean','passed']).to_csv(output/'convergence_edges.csv',index=False)
    core.log('radial_extension_complete',attempted=len(rows),success=sum(r['status']=='ok' for r in rows),
        failed=sum(r['status']!='ok' for r in rows),passing_edges=sum(e['passed'] for e in edges))

if __name__=='__main__':main()
