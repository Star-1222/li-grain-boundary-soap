"""Reproducible Ni SOAP sweep. Run stages sequentially; every stage persists evidence."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
os.environ['MPLBACKEND'] = 'Agg'

import argparse
import csv
import hashlib
import importlib.metadata
import itertools
import json
import platform
import shutil
import sys
import time
import traceback
import warnings
from datetime import datetime, timezone, timedelta
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.linalg import norm as matrix_norm
from scipy.special import gamma
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr
from dscribe.descriptors import SOAP

SOURCE = Path(__file__).resolve().parent
WORKSPACE = SOURCE.parent.parent
ROOT = WORKSPACE / 'output' / 'ni_soap_sweep_20261101'
sys.path.insert(0, str(WORKSPACE / 'tmp' / 'soap_timing'))
from cfg_reader import read_cfg, get_auxiliary

RC = [3.25, 4., 5., 6.]
SIGMA = [.3, .5, .7]
NMAX = [6, 8, 10, 12]
LMAX = [4, 6, 8, 10, 12]
SEED = 20261101
BASELINE = 'r3p25_s0p5_n12_l12'
GROUPS = ['stable_bulk', 'stable_interface', 'metastable_bulk', 'metastable_interface']

def stamp():
    return datetime.now(timezone(timedelta(hours=8))).isoformat(timespec='seconds')

def clean(obj):
    if isinstance(obj, dict): return {str(k):clean(v) for k,v in obj.items()}
    if isinstance(obj, (list,tuple)): return [clean(v) for v in obj]
    if isinstance(obj, np.ndarray): return clean(obj.tolist())
    if isinstance(obj, np.generic): return clean(obj.item())
    if isinstance(obj, float) and not np.isfinite(obj): return None
    if isinstance(obj, Path): return str(obj)
    return obj

def save_json(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(clean(value), ensure_ascii=False, indent=2), encoding='utf-8')
    temporary.replace(path)

def log(event, **kwargs):
    row = dict(time=stamp(), event=event, **kwargs)
    print(json.dumps(clean(row), ensure_ascii=False), flush=True)
    ROOT.mkdir(parents=True, exist_ok=True)
    with (ROOT/'progress.jsonl').open('a',encoding='utf-8') as f:
        f.write(json.dumps(clean(row),ensure_ascii=False)+'\n')
    save_json(ROOT/'status.json', row)

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def config_id(rc, s, n, l):
    def text(v): return f'{v:g}'.replace('.', 'p')
    return f'r{text(rc)}_s{text(s)}_n{n}_l{l}'

def configurations():
    return [dict(config_id=config_id(r,s,n,l), r_cut=r, sigma=s, n_max=n, l_max=l)
            for r,s,n,l in itertools.product(RC,SIGMA,NMAX,LMAX)]

def make_soap(c):
    return SOAP(species=['Ni'], periodic=True, r_cut=c['r_cut'], sigma=c['sigma'],
                n_max=int(c['n_max']), l_max=int(c['l_max']), rbf='gto',
                average='off', sparse=False, dtype='float64', weighting=None,
                compression={'mode':'off','species_weighting':None})

def structures():
    return {s:read_cfg(ROOT/'inputs'/f'ni_{s}.cfg', coordinate_convention='scaled')
            for s in ('stable','metastable')}

def summarize(x):
    x=np.asarray(x); x=x[np.isfinite(x)]
    return dict(n=len(x),median=float(np.median(x)),p05=float(np.quantile(x,.05)),
                p95=float(np.quantile(x,.95)),mean=float(np.mean(x))) if len(x) else {}

def unit(x):
    lengths=np.linalg.norm(x,axis=1)
    if not np.isfinite(x).all() or np.any(lengths<=0):
        raise ValueError('Invalid SOAP vector or zero norm')
    return x/lengths[:,None]

def rho(x,y):
    if np.ptp(x)==0 or np.ptp(y)==0: return float('nan')
    return float(spearmanr(x,y).statistic)

def profile(at, axis):
    z=at.positions[:,axis]
    edges=np.linspace(z.min(),z.max(),201)
    bins=np.digitize(z,edges[1:-1])
    bad=(get_auxiliary(at,'c_CNA')!=1)|(get_auxiliary(at,'c_PTM_1')!=1)
    count=np.bincount(bins,minlength=200)
    fraction=np.bincount(bins,weights=bad,minlength=200)/np.maximum(count,1)
    mask=(fraction>.1)&(count>=50)
    return dict(edges=edges,bins=bins,count=count,fraction=fraction,mask=mask)

def stratified(indices, keys, count, rng):
    """Round-robin across spatial/type strata; selection never uses PE or SOAP."""
    pools={}
    for i in indices: pools.setdefault(tuple(keys[int(i)]),[]).append(int(i))
    names=sorted(pools)
    for k in names: rng.shuffle(pools[k])
    selected=[]
    while len(selected)<count:
        progressed=False
        for k in names:
            if pools[k] and len(selected)<count:
                selected.append(pools[k].pop()); progressed=True
        if not progressed: raise ValueError('Insufficient sampling population')
    return np.asarray(selected,dtype=int)

def pair_indices(labels, seed, per_category=2000):
    rng=np.random.default_rng(seed); ii=[];jj=[];cats=[];names=[]
    for ci,(ga,gb) in enumerate(itertools.combinations_with_replacement(GROUPS,2)):
        a=np.flatnonzero(labels==ga); b=np.flatnonzero(labels==gb)
        if ga==gb:
            u,v=np.triu_indices(len(a),1)
            choice=rng.choice(len(u),min(per_category,len(u)),replace=False)
            i=a[u[choice]];j=a[v[choice]]
        else:
            choice=rng.choice(len(a)*len(b),min(per_category,len(a)*len(b)),replace=False)
            i=a[choice//len(b)];j=b[choice%len(b)]
        ii.extend(i);jj.extend(j);cats.extend([ci]*len(i));names.append(ga+'__'+gb)
    return np.asarray(ii),np.asarray(jj),np.asarray(cats),names

def prepare():
    if (ROOT/'protocol.json').exists():
        log('prepare_already_complete');return
    (ROOT/'inputs').mkdir(parents=True,exist_ok=True)
    evidence=[]
    old=json.loads((WORKSPACE/'tmp/soap_timing/inputs.json').read_text())
    for item in old:
        src=WORKSPACE/'tmp/soap_timing'/item['file']
        assert sha(src)==item['sha256'], f"Input changed: {src}"
        dst=ROOT/'inputs'/src.name; shutil.copy2(src,dst)
        evidence.append(dict(file=dst.name,sha256=sha(dst),bytes=dst.stat().st_size))
    atoms=structures(); rng=np.random.default_rng(SEED)
    eval_records=[]; ref_records=[]; alt_records=[]; profiles={}; geometry={}
    for state,at in atoms.items():
        ps=[profile(at,axis) for axis in range(3)]
        normal=max(range(3),key=lambda k:(ps[k]['mask'].sum(),ps[k]['fraction'].max()))
        assert normal==2, f'Unexpected normal {normal}'
        p=ps[normal]; interface=p['mask'][p['bins']]
        chosen=np.flatnonzero(p['mask'])
        runs=np.split(chosen,np.where(np.diff(chosen)>1)[0]+1)
        runs=[r for r in runs if len(r)]
        assert len(runs)==2, 'Expected two localized interface bin runs'
        edges=p['edges']; z=at.positions[:,2]; period=float(at.cell[2,2])
        run_id=np.full(len(at),-1,dtype=int)
        distance=np.full(len(at),np.inf)
        ranges=[]
        for k,r in enumerate(runs):
            lo=float(edges[r[0]]);hi=float(edges[r[-1]+1]);ranges.append([lo,hi])
            run_id[np.isin(p['bins'],r)]=k
            for shift in (-period,0,period):
                distance=np.minimum(distance,np.maximum.reduce([lo+shift-z,z-hi-shift,np.zeros(len(at))]))
        cna=get_auxiliary(at,'c_CNA').astype(int);ptm=get_auxiliary(at,'c_PTM_1').astype(int)
        ids=get_auxiliary(at,'id').astype(int)
        safe=(~interface)&(distance>10)&(cna==1)&(ptm==1)
        bulk_pool=np.flatnonzero(safe)
        assert len(bulk_pool)>=900
        bulk_keys=[(int(v//20),) for v in z]
        sampled=stratified(bulk_pool,bulk_keys,900,rng)
        bulk_eval=sampled[:500];ref=sampled[500:700];alt=sampled[700:900]
        selected_eval=np.concatenate([bulk_eval,np.flatnonzero(interface)])
        keys=[(int(run_id[i]),int(p['bins'][i]),int(cna[i]),int(ptm[i])) for i in range(len(at))]
        pilot_interface=stratified(np.flatnonzero(interface),keys,300,rng)
        pilot_bulk=rng.choice(bulk_eval,200,replace=False)
        pilot_set=set(np.concatenate([pilot_bulk,pilot_interface]).tolist())
        pert_set=set(np.concatenate([rng.choice(pilot_bulk,30,replace=False),
                                    stratified(pilot_interface,keys,30,rng)]).tolist())
        for role,indices,records in [('evaluation',selected_eval,eval_records),('reference',ref,ref_records),('alternate_reference',alt,alt_records)]:
            for index in indices:
                i=int(index)
                records.append(dict(structure=state,atom_index=i,atom_id=int(ids[i]),role=role,
                    region='interface' if interface[i] else 'bulk', group=state+('_interface' if interface[i] else '_bulk'),
                    interface_run=int(run_id[i]),distance_to_interface_A=float(distance[i]),
                    x=float(at.positions[i,0]),y=float(at.positions[i,1]),z=float(z[i]),
                    cna=int(cna[i]),ptm=int(ptm[i]),pilot=role=='evaluation' and i in pilot_set,
                    perturb=role=='evaluation' and i in pert_set))
        profiles[state]=dict(axis='z',edges=edges,counts=p['count'],fraction=p['fraction'],interface_bins=p['mask'])
        geometry[state]=dict(atoms=len(at),cell=at.cell.array.tolist(),interface_count=int(interface.sum()),
                             runs_A=ranges,safe_bulk_count=len(bulk_pool),
                             minimum_reference_margin_A=float(distance[ref].min()))
    frame=pd.DataFrame(eval_records+ref_records+alt_records)
    frame['record_row']=np.arange(len(frame));frame.to_csv(ROOT/'population.csv',index=False)
    assert frame[frame.role=='evaluation'].shape[0]==4659
    assert frame.pilot.sum()==1000 and frame.perturb.sum()==120
    assert not frame.duplicated(['structure','atom_index']).any()
    pilot=frame[frame.pilot].reset_index(drop=True)
    ii,jj,cat,names=pair_indices(pilot.group.to_numpy(),SEED+1)
    np.savez(ROOT/'fixed_pairs.npz',i=ii,j=jj,category=cat,names=names)
    save_json(ROOT/'geometry.json',geometry);save_json(ROOT/'profiles.json',profiles)
    protocol=dict(created_at=stamp(),seed=SEED,grid=dict(r_cut=RC,sigma=SIGMA,n_max=NMAX,l_max=LMAX),
        count=240,implementation=dict(rbf='gto',species=['Ni'],periodic=True,average='off',
        weighting=None,compression='off',dtype='float64',threads_per_process=1),
        distance_definition='raw Euclidean and Euclidean on unit DScribe unique-term output, no symmetry reweighting',
        reference='pooled 200 deep bulk per state; no PEATOM filter; alternate independent 200 per state',
        evaluation_count=4659,pilot_count=1000,reference_count=400,pair_count=len(ii),
        masks='200 bins over observed position range, non-dominant CNA/PTM >0.1 and bin count >=50; historical mask reproduced',
        bulk_margin_A=10.,max_neighbor_search_A=6+.7*np.sqrt(-2*np.log(.001)),
        convergence=dict(spearman_min=.98,median_relative_change_max=.05,knn10_overlap_min=.8,
                         require_two_successive_edges_per_axis=True),
        perturbation=dict(random_rms_A=[.01,.03,.05],seeds=[SEED+101,SEED+102,SEED+103],
            isotropic_strain=[-.01,-.005,.005,.01],invariance_relative_tolerance=1e-8,
            note='synthetic probes, not temperatures or MD trajectories'),
        selection='no PEATOM or stable/metastable separation optimization; separate raw/norm convergence evidence; four rc representatives + alternate sigma + historical baseline',
        source_repository='https://github.com/Star-1222/li-grain-boundary-soap',
        repository_main_at_start='b18a6c59d115f29ac33189275a151eca7512bf1b',
        input_evidence=evidence,versions={p:importlib.metadata.version(p) for p in ['numpy','scipy','ase','dscribe','pandas','matplotlib']},
        python=platform.python_version(),source_sha256=sha(__file__),
        limitations=['Two CFG configurations only; atoms and periodic interfaces are not independent GB replicas.',
                     'CNA/PTM and fixed historical masks provide sampling labels, not truth for all structural differences.',
                     'No universal physical tolerance or unique best length scale is assumed.'])
    save_json(ROOT/'protocol.json',protocol)
    shutil.copy2(__file__,ROOT/'run_sweep_source.py')
    pd.DataFrame(configurations()).to_csv(ROOT/'grid.csv',index=False)
    log('prepared',geometry=geometry,pilot=len(pilot),reference=400)

def population():
    return pd.read_csv(ROOT/'population.csv')

def calc_vectors(soap,atoms,records):
    blocks=[]
    for state,sub in records.groupby('structure',sort=False):
        blocks.append(soap.create(atoms[state],centers=sub.atom_index.to_numpy(dtype=int),n_jobs=1))
    # records are grouped by structure at each call; preserve their explicit order.
    assert len(records.structure.unique())==1 or list(records.structure)==sorted(records.structure,key={'stable':0,'metastable':1}.get)
    return np.vstack(blocks)

def ordered_records(frame):
    return pd.concat([frame[frame.structure==s] for s in ('stable','metastable')],ignore_index=True)

def basis_health(soap):
    conditions=[];residuals=[]
    for l,alpha in enumerate(soap._alphas):
        s=.5*gamma(l+1.5)*(alpha[:,None]+alpha[None,:])**(-l-1.5)
        b=soap._betas[l]
        conditions.append(float(np.linalg.cond(s)))
        residuals.append(float(np.linalg.norm(b@s@b.T-np.eye(len(alpha)),ord=np.inf)))
    return dict(basis_condition_max=max(conditions),basis_residual_max=max(residuals),
                basis_conditions=conditions,basis_residuals=residuals)

def metrics(x,ref,pilot,pairs):
    result={}; summaries=[]
    for metric,xx,rr in [('raw',x,ref),('norm',unit(x),unit(ref))]:
        d=cdist(xx,xx); pair=d[pairs['i'],pairs['j']]
        np.fill_diagonal(d,np.inf)
        nn=np.argpartition(d,9,axis=1)[:,:10]
        dr=cdist(xx,rr); nearest=np.partition(dr,19,axis=1)[:,:20];nearest.sort(axis=1)
        result['pair_'+metric]=pair
        result['nn_'+metric]=nn
        for k in [1,5,10,20]:result[f'D{ k }_'+metric]=nearest[:,:k].mean(axis=1)
        for group in GROUPS:
            mask=pilot.group.to_numpy()==group
            summaries.append(dict(metric=metric,group=group,**summarize(result['D10_'+metric][mask])))
        for ci,name in enumerate(pairs['names']):
            result[f'pair_median_{metric}_{ci}']=np.median(pair[pairs['category']==ci])
    result['norms']=np.linalg.norm(x,axis=1)
    a=np.linalg.norm(x[pairs['i']],axis=1);b=np.linalg.norm(x[pairs['j']],axis=1)
    residual=result['pair_raw']**2-((a-b)**2+a*b*result['pair_norm']**2)
    result['identity_max_absolute_residual']=np.max(np.abs(residual))
    result['raw_norm_D_spearman']=rho(result['D10_raw'],result['D10_norm'])
    return result,summaries

def sweep():
    atoms=structures(); frame=population()
    pilot=ordered_records(frame[frame.pilot]);refs=ordered_records(frame[frame.role=='reference'])
    pairs=dict(np.load(ROOT/'fixed_pairs.npz'))
    for name in ['vectors','metrics','records']:(ROOT/name).mkdir(exist_ok=True)
    begin=time.perf_counter()
    for count,c in enumerate(configurations(),1):
        cid=c['config_id']; record_path=ROOT/'records'/f'{cid}.json'
        if record_path.exists():continue
        t=time.perf_counter(); row=dict(**c,started_at=stamp())
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always');soap=make_soap(c)
                health=basis_health(soap)
                x=calc_vectors(soap,atoms,pilot);ref=calc_vectors(soap,atoms,refs)
            create_seconds=time.perf_counter()-t
            unit(x);unit(ref)
            np.savez(ROOT/'vectors'/f'{cid}.npz',evaluation=x,reference=ref)
            result,stats=metrics(x,ref,pilot,pairs)
            np.savez(ROOT/'metrics'/f'{cid}.npz',**result)
            row.update(status='ok',dimension=x.shape[1],create_s=create_seconds,
                total_s=time.perf_counter()-t,finite=True,
                warnings=sorted(set(str(w.message) for w in caught)),**health,
                identity_max_absolute_residual=float(result['identity_max_absolute_residual']),
                raw_norm_D_spearman=float(result['raw_norm_D_spearman']),D_summary=stats)
        except Exception as exc:
            row.update(status='failed',error=repr(exc),traceback=traceback.format_exc(),total_s=time.perf_counter()-t)
        save_json(record_path,row)
        if count%10==0 or count==1 or row['status']=='failed':
            log('sweep_progress',done=count,total=240,elapsed_s=time.perf_counter()-begin,last=cid,status=row['status'])
    rows=[json.loads(p.read_text(encoding='utf-8')) for p in sorted((ROOT/'records').glob('*.json'))]
    flat=[{k:v for k,v in r.items() if not isinstance(v,(list,dict))} | {'warning_count':len(r.get('warnings',[]))} for r in rows]
    pd.DataFrame(flat).to_csv(ROOT/'sweep_summary.csv',index=False)
    log('sweep_complete',ok=sum(r['status']=='ok' for r in rows),failed=sum(r['status']!='ok' for r in rows))

def overlap(a,b):
    return np.asarray([len(set(x)&set(y))/10 for x,y in zip(a,b)])

def convergence():
    valid={r['config_id']:r for r in pd.read_csv(ROOT/'sweep_summary.csv').to_dict('records') if r['status']=='ok'}
    pairs=dict(np.load(ROOT/'fixed_pairs.npz'));pilot=ordered_records(population().query('pilot'))
    cache={cid:dict(np.load(ROOT/'metrics'/f'{cid}.npz')) for cid in valid}
    rows=[]
    for rc,s in itertools.product(RC,SIGMA):
        edges=[]
        for l in LMAX:
            edges.extend(('n',config_id(rc,s,n1,l),config_id(rc,s,n2,l)) for n1,n2 in zip(NMAX[:-1],NMAX[1:]))
        for n in NMAX:
            edges.extend(('l',config_id(rc,s,n,l1),config_id(rc,s,n,l2)) for l1,l2 in zip(LMAX[:-1],LMAX[1:]))
        for axis,lo,hi in edges:
            if lo not in valid or hi not in valid:continue
            for metric in ['raw','norm']:
                a=cache[lo]['pair_'+metric];b=cache[hi]['pair_'+metric]
                correlations=[];changes=[]
                for ci,name in enumerate(pairs['names']):
                    select=pairs['category']==ci
                    corr=rho(a[select],b[select])
                    scale=max(float(np.median(b[select])),1e-12)
                    change=float(np.median(np.abs(a[select]-b[select]))/scale)
                    correlations.append(corr);changes.append(change)
                near=overlap(cache[lo]['nn_'+metric],cache[hi]['nn_'+metric])
                group_overlaps=[float(near[pilot.group.to_numpy()==g].mean()) for g in GROUPS]
                row=dict(axis=axis,lower=lo,upper=hi,metric=metric,r_cut=rc,sigma=s,
                    spearman_min=float(np.nanmin(correlations)),relative_change_max=max(changes),
                    nn_overlap_mean=float(near.mean()),nn_overlap_group_min=min(group_overlaps),
                    D_spearman=rho(cache[lo]['D10_'+metric],cache[hi]['D10_'+metric]))
                row['passed']=row['spearman_min']>=.98 and row['relative_change_max']<=.05 and row['nn_overlap_mean']>=.8
                rows.append(row)
    edges=pd.DataFrame(rows);edges.to_csv(ROOT/'convergence_edges.csv',index=False)
    summaries=[]
    for cid,c in valid.items():
        cfg=next(x for x in configurations() if x['config_id']==cid)
        for metric in ['raw','norm']:
            evidence=edges[(edges.metric==metric)&((edges.lower==cid)|(edges.upper==cid))]
            axes={}
            for axis,values in [('n',NMAX),('l',LMAX)]:
                position=values.index(int(cfg['n_max' if axis=='n' else 'l_max']))
                # At an interior point use incoming+outgoing; at the upper edge use last two edges.
                adjacent=evidence[evidence.axis==axis]
                if position==len(values)-1:
                    if axis=='n':previous=config_id(cfg['r_cut'],cfg['sigma'],values[-2],cfg['l_max'])
                    else:previous=config_id(cfg['r_cut'],cfg['sigma'],cfg['n_max'],values[-2])
                    older=edges[(edges.metric==metric)&(edges.axis==axis)&(edges.upper==previous)]
                    adjacent=pd.concat([adjacent,older])
                axes[axis]=len(adjacent)>=2 and bool(adjacent.passed.all())
            deficit=[]
            for e in evidence.to_dict('records'):
                deficit.append(max(0,.98-e['spearman_min'])/.02+max(0,e['relative_change_max']-.05)/.05+max(0,.8-e['nn_overlap_mean'])/.2)
            summaries.append(dict(**cfg,metric=metric,plateau_n=axes['n'],plateau_l=axes['l'],
                converged=axes['n'] and axes['l'],local_pass_fraction=float(evidence.passed.mean()),
                deficit=float(np.mean(deficit)) if deficit else 1e9,
                dimension=int(c['dimension']),warning_count=int(c['warning_count']),
                basis_residual_max=c['basis_residual_max']))
    table=pd.DataFrame(summaries);table.to_csv(ROOT/'convergence_summary.csv',index=False)
    joint=[]
    for cid,part in table.groupby('config_id'):
        row=part.iloc[0].to_dict();row.update(both_converged=bool(part.converged.all()),
            converged_metrics=int(part.converged.sum()),mean_deficit=float(part.deficit.mean()),
            metric_evidence=part[['metric','converged','plateau_n','plateau_l','deficit']].to_dict('records'))
        joint.append(row)
    def priority(r):
        return (-r['converged_metrics'],r['warning_count']>0,r['mean_deficit'],r['dimension'],r['n_max'],r['l_max'])
    selected=[]
    for rc in RC:
        options=sorted([r for r in joint if r['r_cut']==rc and r['config_id']!=BASELINE],key=priority)
        pick=options[0].copy();pick['selection_reason']='best convergence evidence at this r_cut; provisional audit candidate';selected.append(pick)
    used={r['config_id'] for r in selected};sigmas={r['sigma'] for r in selected}
    options=[r for r in joint if r['config_id'] not in used and r['config_id']!=BASELINE]
    options.sort(key=lambda r:(r['sigma'] in sigmas,*priority(r)))
    pick=options[0].copy();pick['selection_reason']='alternate smoothing scale';selected.append(pick)
    pick=next(r.copy() for r in joint if r['config_id']==BASELINE)
    pick['selection_reason']='historical baseline, retained for diagnosis even if not converged';selected.append(pick)
    save_json(ROOT/'candidates.json',selected)
    log('convergence_complete',edges=len(edges),converged_by_metric=table.groupby('metric').converged.sum().to_dict(),
        candidates=[dict(config_id=r['config_id'],both_converged=r['both_converged']) for r in selected])

def main():
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','sweep','convergence'])
    args=parser.parse_args();globals()[args.stage]()

if __name__=='__main__':main()
