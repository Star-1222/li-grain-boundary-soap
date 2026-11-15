"""Perturbation probes and held-out population checks, never used to retune pilot selection."""
import run_sweep as core
import argparse
import json
import time
import warnings
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.spatial import cKDTree

ROOT=core.ROOT

def candidates():
    return json.loads((ROOT/'candidates.json').read_text(encoding='utf-8'))

def make_transform(at,kind,amplitude,seed):
    result=at.copy(); inverse=None
    if kind=='random':
        displacement=np.random.default_rng(seed).normal(size=(len(at),3))
        displacement-=displacement.mean(axis=0)
        displacement*=amplitude/np.sqrt(np.mean(np.sum(displacement**2,axis=1)))
        result.positions+=displacement
    elif kind=='strain':
        result.set_cell(at.cell.array*(1+amplitude),scale_atoms=True)
    elif kind=='translation':
        result.positions+=np.array([2.137,-1.419,3.071]);result.wrap()
    elif kind=='rotation':
        q,_=np.linalg.qr(np.random.default_rng(seed).normal(size=(3,3)))
        if np.linalg.det(q)<0:q[:,0]*=-1
        result.set_cell(at.cell.array@q.T,scale_atoms=False)
        result.positions=at.positions@q.T
    elif kind=='permutation':
        perm=np.random.default_rng(seed).permutation(len(at));result=at[perm]
        inverse=np.argsort(perm)
    return result,inverse

def perturb():
    frame=core.population();pilot=core.ordered_records(frame[frame.pilot])
    selected=np.flatnonzero(pilot.perturb.to_numpy())
    anchors=pilot.iloc[selected].reset_index(drop=True)
    atoms=core.structures();rows=[];stats=[]
    (ROOT/'perturb').mkdir(exist_ok=True)
    transforms=[('random',amplitude,seed) for amplitude in [.01,.03,.05] for seed in [core.SEED+101,core.SEED+102,core.SEED+103]]
    transforms += [('strain',amplitude,core.SEED+201) for amplitude in [-.01,-.005,.005,.01]]
    transforms += [(kind,0.,core.SEED+301) for kind in ['translation','rotation','permutation']]
    for position,c in enumerate(candidates(),1):
        cid=c['config_id'];path=ROOT/'perturb'/f'{cid}.csv'
        if path.exists():
            rows.extend(pd.read_csv(path).to_dict('records'));continue
        t=__import__('time').perf_counter()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore');soap=core.make_soap(c)
        saved=np.load(ROOT/'vectors'/f'{cid}.npz');x=saved['evaluation'];ref=saved['reference']
        original=x[selected];u=core.unit(x);ou=u[selected]
        m=dict(np.load(ROOT/'metrics'/f'{cid}.npz'))
        fixed=dict(np.load(ROOT/'fixed_pairs.npz'));names=list(fixed['names'])
        scales={}
        for state in ['stable','metastable']:
            cat=names.index(state+'_interface__'+state+'_interface')
            for metric in ['raw','norm']:
                scales[state,metric]=float(np.median(m['pair_'+metric][fixed['category']==cat]))
        local_rows=[]
        for kind,amplitude,seed in transforms:
            blocks=[]
            for state,sub in anchors.groupby('structure',sort=False):
                changed,inverse=make_transform(atoms[state],kind,amplitude,seed)
                idx=sub.atom_index.to_numpy(dtype=int)
                if inverse is not None:idx=inverse[idx]
                blocks.append(soap.create(changed,centers=idx,n_jobs=1))
            modified=np.vstack(blocks);um=core.unit(modified)
            for metric,xx,base,all_original,rr in [('raw',modified,original,x,ref),('norm',um,ou,u,core.unit(ref))]:
                delta=np.linalg.norm(xx-base,axis=1)
                compare=cdist(xx,all_original);compare[np.arange(len(selected)),selected]=np.inf
                nearest=np.argpartition(compare,9,axis=1)[:,:10]
                retention=core.overlap(nearest,m['nn_'+metric][selected])
                dr=cdist(xx,rr);d10=np.partition(dr,9,axis=1)[:,:10].mean(axis=1)
                relative=delta/np.maximum(np.linalg.norm(base,axis=1),1e-300)
                for i,row in anchors.iterrows():
                    local_rows.append(dict(config_id=cid,kind=kind,amplitude=amplitude,seed=seed,metric=metric,
                        structure=row.structure,region=row.region,group=row.group,atom_id=int(row.atom_id),
                        atom_index=int(row.atom_index),distance=float(delta[i]),relative_vector_change=float(relative[i]),
                        nn10_retention=float(retention[i]),D10=float(d10[i]),
                        original_D10=float(m['D10_'+metric][selected[i]]),
                        distance_over_interface_pair_median=float(delta[i]/max(scales[row.structure,metric],1e-300))))
        pd.DataFrame(local_rows).to_csv(path,index=False);rows.extend(local_rows)
        core.log('perturb_progress',done=position,total=len(candidates()),config_id=cid,seconds=__import__('time').perf_counter()-t)
    allrows=pd.DataFrame(rows);allrows.to_csv(ROOT/'perturbation_atomwise.csv',index=False)
    for keys,g in allrows.groupby(['config_id','kind','amplitude','metric','group'],sort=False):
        row=dict(zip(['config_id','kind','amplitude','metric','group'],keys))
        row.update(core.summarize(g.distance.to_numpy()))
        row.update(relative_change_max=float(g.relative_vector_change.max()),
            nn10_retention_mean=float(g.nn10_retention.mean()),
            noise_ratio_median=float(g.distance_over_interface_pair_median.median()),
            noise_ratio_p95=float(g.distance_over_interface_pair_median.quantile(.95)))
        stats.append(row)
    pd.DataFrame(stats).to_csv(ROOT/'perturbation_summary.csv',index=False)
    invariance=allrows[allrows.kind.isin(['translation','rotation','permutation'])]
    core.save_json(ROOT/'invariance.json',dict(relative_tolerance=1e-8,
        max_relative_vector_change=float(invariance.relative_vector_change.max()),
        passed=bool((invariance.relative_vector_change<1e-8).all())))
    core.log('perturb_complete',rows=len(allrows),invariance_max=float(invariance.relative_vector_change.max()))

def correlation_rows(cid,metric,variants,baseline,records):
    rows=[]
    for name,values in variants.items():
        for group in ['all']+core.GROUPS:
            for subset in ['all','heldout']:
                mask=np.ones(len(records),dtype=bool) if group=='all' else records.group.to_numpy()==group
                if subset=='heldout':mask &= ~records.pilot.to_numpy()
                if mask.sum()<3:continue
                x=baseline[mask];y=values[mask]
                rows.append(dict(config_id=cid,metric=metric,variant=name,group=group,subset=subset,
                    n=int(mask.sum()),spearman=core.rho(x,y),
                    median_relative_change=float(np.median(np.abs(x-y))/max(np.median(x),1e-12))))
    return rows

def review():
    frame=core.population();evaluate=core.ordered_records(frame[frame.role=='evaluation'])
    refs=core.ordered_records(frame[frame.role=='reference']);alt=core.ordered_records(frame[frame.role=='alternate_reference'])
    atoms=core.structures();summaries=[];robustness=[];geometry_rows=[]
    (ROOT/'review').mkdir(exist_ok=True)
    pairs=core.pair_indices(evaluate.group.to_numpy(),core.SEED+401)
    np.savez(ROOT/'review_fixed_pairs.npz',i=pairs[0],j=pairs[1],category=pairs[2],names=pairs[3])
    for position,c in enumerate(candidates(),1):
        cid=c['config_id'];path=ROOT/'review'/f'{cid}.csv'
        summary_path=ROOT/'review'/f'{cid}_summary.json'
        if summary_path.exists():
            saved=json.loads(summary_path.read_text(encoding='utf-8'))
            summaries.extend(saved['summaries']);robustness.extend(saved['robustness']);geometry_rows.extend(saved['geometry']);continue
        t=time.perf_counter()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore');soap=core.make_soap(c)
            x=core.calc_vectors(soap,atoms,evaluate)
            altx=core.calc_vectors(soap,atoms,alt)
        saved=np.load(ROOT/'vectors'/f'{cid}.npz');ref=saved['reference']
        full=evaluate.copy();full['soap_norm']=np.linalg.norm(x,axis=1)
        pair_values={};local_robust=[];local_summary=[];local_geometry=[]
        for metric,xx,rr,aa in [('raw',x,ref,altx),('norm',core.unit(x),core.unit(ref),core.unit(altx))]:
            dr=cdist(xx,rr);ordered=np.partition(dr,19,axis=1)[:,:20];ordered.sort(axis=1)
            for k in [1,5,10,20]:full[f'D{k}_{metric}']=ordered[:,:k].mean(axis=1)
            da=cdist(xx,aa);alternative=np.partition(da,9,axis=1)[:,:10].mean(axis=1)
            perstate=np.empty(len(evaluate))
            for state in ['stable','metastable']:
                mask=evaluate.structure.to_numpy()==state
                d=dr[mask][:,refs.structure.to_numpy()==state]
                perstate[mask]=np.partition(d,9,axis=1)[:,:10].mean(axis=1)
            full['D10_alternate_'+metric]=alternative;full['D10_perstate_'+metric]=perstate
            variants={f'k{k}':full[f'D{k}_{metric}'].to_numpy() for k in [1,5,20]}
            variants.update(alternate_reference=alternative,perstate_reference=perstate)
            local_robust.extend(correlation_rows(cid,metric,variants,full['D10_'+metric].to_numpy(),evaluate))
            pair_values[metric]=np.asarray([np.linalg.norm(xx[i:i+1]-xx[j:j+1]) for i,j in zip(pairs[0],pairs[1])])
            for group in core.GROUPS:
                for subset in ['all','pilot','heldout']:
                    mask=evaluate.group.to_numpy()==group
                    if subset=='pilot':mask &= evaluate.pilot.to_numpy()
                    elif subset=='heldout':mask &= ~evaluate.pilot.to_numpy()
                    local_summary.append(dict(config_id=cid,metric=metric,group=group,subset=subset,
                        **core.summarize(full.loc[mask,'D10_'+metric].to_numpy())))
            for state in ['stable','metastable']:
                for run in [0,1]:
                    mask=(evaluate.structure.to_numpy()==state)&(evaluate.interface_run.to_numpy()==run)
                    local_summary.append(dict(config_id=cid,metric=metric,group=state+f'_interface_run{run}',subset='all',
                        **core.summarize(full.loc[mask,'D10_'+metric].to_numpy())))
        # Coordination counts use a declared real-space radius, separate from padded SOAP search.
        counts=[]
        for state,sub in evaluate.groupby('structure',sort=False):
            at=atoms[state];lengths=at.cell.lengths()
            tree=cKDTree(at.positions%lengths,boxsize=lengths)
            counts.extend(tree.query_ball_point(sub[['x','y','z']].to_numpy()%lengths,c['r_cut'],return_length=True)-1)
        full['coordination_nominal_rc']=counts
        for group in core.GROUPS:
            mask=full.group.to_numpy()==group
            part=full[mask]
            local_geometry.append(dict(config_id=cid,group=group,n=len(part),
                rho_norm_coordination=core.rho(part.soap_norm,part.coordination_nominal_rc),
                rho_Draw_Dnorm=core.rho(part.D10_raw,part.D10_norm),
                rho_norm_Draw=core.rho(part.soap_norm,part.D10_raw),
                rho_norm_Dnorm=core.rho(part.soap_norm,part.D10_norm)))
        full.to_csv(path,index=False)
        np.savez(ROOT/'review'/f'{cid}_vectors.npz',evaluation=x,alternate_reference=altx)
        np.savez(ROOT/'review'/f'{cid}_pairs.npz',**pair_values)
        # Explicit squared-distance decomposition against the common raw mean reference.
        mean=ref.mean(axis=0);norm_mean=np.linalg.norm(mean);lengths=np.linalg.norm(x,axis=1)
        raw2=np.sum((x-mean)**2,axis=1)
        magnitude=(lengths-norm_mean)**2
        fraction=np.divide(magnitude,raw2,out=np.full_like(raw2,np.nan),where=raw2>1e-20)
        full_decomp=evaluate[['structure','atom_id','group','pilot']].copy()
        full_decomp['magnitude_fraction_squared_distance']=fraction
        full_decomp.to_csv(ROOT/'review'/f'{cid}_decomposition.csv',index=False)
        for group in core.GROUPS:
            local_geometry.append(dict(config_id=cid,group=group,
                magnitude_fraction_median=float(np.nanmedian(fraction[evaluate.group.to_numpy()==group]))))
        evidence=dict(summaries=local_summary,robustness=local_robust,geometry=local_geometry,
                      review_seconds=time.perf_counter()-t)
        core.save_json(summary_path,evidence)
        summaries.extend(local_summary);robustness.extend(local_robust);geometry_rows.extend(local_geometry)
        core.log('review_progress',done=position,total=len(candidates()),config_id=cid,seconds=time.perf_counter()-t)
    pd.DataFrame(summaries).to_csv(ROOT/'review_summary.csv',index=False)
    pd.DataFrame(robustness).to_csv(ROOT/'reference_k_robustness.csv',index=False)
    pd.DataFrame(geometry_rows).to_csv(ROOT/'physical_audit.csv',index=False)
    compare=[]
    all_candidates=candidates()
    for a,b in __import__('itertools').combinations(all_candidates,2):
        fa=pd.read_csv(ROOT/'review'/f"{a['config_id']}.csv")
        fb=pd.read_csv(ROOT/'review'/f"{b['config_id']}.csv")
        for metric in ['raw','norm']:
            for group in core.GROUPS:
                for subset in ['pilot','heldout']:
                    mask=fa.group.to_numpy()==group
                    mask &= fa.pilot.to_numpy() if subset=='pilot' else ~fa.pilot.to_numpy()
                    compare.append(dict(first=a['config_id'],second=b['config_id'],metric=metric,group=group,subset=subset,
                        spearman=core.rho(fa.loc[mask,'D10_'+metric],fb.loc[mask,'D10_'+metric])))
    pd.DataFrame(compare).to_csv(ROOT/'candidate_agreement.csv',index=False)
    core.log('review_complete',evaluation=len(evaluate),heldout=int((~evaluate.pilot).sum()),candidates=len(all_candidates))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['perturb','review'])
    args=parser.parse_args();globals()[args.stage]()
