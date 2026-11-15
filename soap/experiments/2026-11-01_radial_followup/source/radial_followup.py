"""Six polynomial Ni tests and interpretable FCC probes; main-grid artifacts remain immutable."""
import run_sweep as core
import argparse
import inspect
import itertools
import json
import shutil
import time
import traceback
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from ase.build import bulk
from ase.io import write
from dscribe.descriptors import SOAP

ROOT=core.WORKSPACE/'output'/'ni_soap_radial_followup_20261101'
MAIN=core.ROOT

def log(event,**kwargs):
    row=dict(time=core.stamp(),event=event,**kwargs)
    print(json.dumps(core.clean(row),ensure_ascii=False),flush=True)
    ROOT.mkdir(parents=True,exist_ok=True)
    with (ROOT/'progress.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps(core.clean(row),ensure_ascii=False)+'\n')
    core.save_json(ROOT/'status.json',row)

def configurations():
    return [dict(config_id=core.config_id(rc,.3,n,10),r_cut=rc,sigma=.3,n_max=n,l_max=10)
            for rc,n in itertools.product([4.,5.],[6,8,10])]

def make_soap(c,rbf):
    return SOAP(species=['Ni'],periodic=True,r_cut=c['r_cut'],sigma=c['sigma'],
        n_max=c['n_max'],l_max=c['l_max'],rbf=rbf,average='off',sparse=False,
        dtype='float64',weighting=None,compression={'mode':'off','species_weighting':None})

def polynomial_health(soap,c):
    n=c['n_max'];rc=c['r_cut']
    i=np.arange(1,n+1)[:,None];j=np.arange(1,n+1)[None,:]
    s=2*rc**(7+i+j)/((5+i+j)*(6+i+j)*(7+i+j))
    # Inspect the actual DScribe tabulated radial functions on its 100-point GL grid.
    rx,gss=soap.get_basis_poly(rc,n)
    nodes,weights=np.polynomial.legendre.leggauss(100)
    assert np.max(np.abs(rx-rc*.5*(nodes+1)))<1e-12
    gram=(gss*(weights*rc*.5*rx**2)[None,:])@gss.T
    diagonal=np.sqrt(np.diag(s));scaled=s/diagonal[:,None]/diagonal[None,:]
    return dict(basis_condition_max=float(np.linalg.cond(s)),
        scaled_overlap_condition=float(np.linalg.cond(scaled)),
        basis_residual_max=float(np.linalg.norm(gram-np.eye(n),ord=np.inf)),
        gram_eigenvalue_min=float(np.linalg.eigvalsh(gram).min()),
        health_method='actual DScribe radial functions integrated with matching 100-point Gauss-Legendre quadrature')

def prepare():
    if (ROOT/'protocol.json').exists():log('already_prepared');return
    ROOT.mkdir(parents=True,exist_ok=True)
    files=['protocol.json','population.csv','fixed_pairs.npz']
    provenance=[]
    for name in files:
        path=MAIN/name
        provenance.append(dict(path=str(path),sha256=core.sha(path)))
    for c in configurations():
        for folder in ['vectors','metrics']:
            p=MAIN/folder/(c['config_id']+'.npz');assert p.exists()
            provenance.append(dict(path=str(p),sha256=core.sha(p)))
    protocol=dict(registered_at=core.stamp(),purpose='diagnose radial representation and metric interpretation; no retrospective relaxation of main-grid criteria',
        new_polynomial_configurations=configurations(),reused_gto_configurations=6,
        main_distances=['raw','norm'],main_weights='unweighted DScribe unique-term vector, as in original sweep',
        supplemental_weights='sqrt(2) on off-diagonal n,n_prime terms; supplementary Frobenius metric only',
        basis_comparison='compare induced distances, neighborhoods and D, never coordinates across different bases',
        convergence_criteria=dict(spearman_min=.98,median_relative_change_max=.05,knn10_overlap_min=.8),
        frozen_population='reuse original 1000 pilot centers and 400 disjoint references; no new sampling or energy filter',
        fcc_probes=dict(lattice_constant='estimate from pooled deep-bulk nearest-12 median distance times sqrt(2)',
            conventional_cell_repeat=[8,8,8],isotropic_strain=[-.02,-.01,-.005,.005,.01,.02],
            angular_degrees=[-10,-5,-2,-1,-.5,.5,1,2,5,10],first_shell_vacancies=[1,2],
            invariance=['translation','rotation','same-element permutation']),
        interpretation='controlled geometry responses are not thermodynamic labels, temperatures, or independent validation',
        main_input_evidence=provenance,source_sha256=core.sha(__file__))
    core.save_json(ROOT/'protocol.json',protocol)
    (ROOT/'implementation').mkdir(exist_ok=True)
    (ROOT/'implementation'/'installed_get_basis_poly.py').write_text(inspect.getsource(SOAP.get_basis_poly),encoding='utf-8')
    shutil.copy2(MAIN/'population.csv',ROOT/'population.csv')
    shutil.copy2(MAIN/'fixed_pairs.npz',ROOT/'fixed_pairs.npz')
    log('prepared',polynomial=6,reused_gto=6)

def scan():
    atoms=core.structures();frame=core.population()
    pilot=core.ordered_records(frame[frame.pilot]);refs=core.ordered_records(frame[frame.role=='reference'])
    pairs=dict(np.load(MAIN/'fixed_pairs.npz'));records=[]
    for folder in ['vectors','metrics','records']:(ROOT/folder).mkdir(exist_ok=True)
    for k,c in enumerate(configurations(),1):
        cid=c['config_id'];record=ROOT/'records'/f'polynomial_{cid}.json'
        if record.exists():records.append(json.loads(record.read_text(encoding='utf-8')));continue
        begin=time.perf_counter();row=dict(**c,rbf='polynomial');caught=[]
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always');soap=make_soap(c,'polynomial')
                row.update(polynomial_health(soap,c))
                x=core.calc_vectors(soap,atoms,pilot);ref=core.calc_vectors(soap,atoms,refs)
            core.unit(x);core.unit(ref)
            np.savez(ROOT/'vectors'/f'polynomial_{cid}.npz',evaluation=x,reference=ref)
            m,stats=core.metrics(x,ref,pilot,pairs)
            np.savez(ROOT/'metrics'/f'polynomial_{cid}.npz',**m)
            row.update(status='ok',dimension=x.shape[1],finite=True,D_summary=stats)
        except Exception as exc:row.update(status='failed',error=repr(exc),traceback=traceback.format_exc())
        row.update(warnings=sorted(set(str(w.message) for w in caught)),seconds=time.perf_counter()-begin)
        core.save_json(record,row);records.append(row)
        log('scan_progress',done=k,total=6,config_id=cid,status=row['status'],seconds=row['seconds'],basis_residual=row.get('basis_residual_max'))
    pd.DataFrame([{k:v for k,v in r.items() if not isinstance(v,(list,dict))}|{'warning_count':len(r['warnings'])} for r in records]).to_csv(ROOT/'polynomial_summary.csv',index=False)
    log('scan_complete',success=sum(r['status']=='ok' for r in records),failed=sum(r['status']=='failed' for r in records))

def vectors(rbf,c):
    folder=MAIN if rbf=='gto' else ROOT
    prefix='' if rbf=='gto' else 'polynomial_'
    return dict(np.load(folder/'vectors'/f"{prefix}{c['config_id']}.npz"))

def diagonal_weights(c):
    n=c['n_max']
    return np.tile([1. if a==b else np.sqrt(2.) for a in range(n) for b in range(a,n)],c['l_max']+1)

def metric_set(x,ref,c,weighting):
    if weighting=='symmetric_frobenius':
        w=diagonal_weights(c);x=x*w;ref=ref*w
    return {metric:(xx,rr) for metric,xx,rr in [('raw',x,ref),('norm',core.unit(x),core.unit(ref))]}

def relative_positions(at,center):
    relative=(at.positions-at.positions[center])@np.linalg.inv(at.cell.array)
    relative-=np.round(relative)
    return relative@at.cell.array

def make_fcc_cases():
    atoms=core.structures();frame=core.population();estimates=[];by_state={}
    for state in ['stable','metastable']:
        at=atoms[state];lengths=at.cell.lengths();tree=cKDTree(at.positions%lengths,boxsize=lengths)
        idx=frame[(frame.structure==state)&(frame.role=='reference')].atom_index.to_numpy(dtype=int)
        distances,_=tree.query(at.positions[idx]%lengths,k=13)
        estimates.extend(distances[:,1:].ravel()*np.sqrt(2))
        by_state[state]=float(np.median(distances[:,1:])*np.sqrt(2))
    a=float(np.median(estimates));original=bulk('Ni','fcc',a=a,cubic=True).repeat((8,8,8))
    original.set_array('original_index',np.arange(len(original)))
    center=int(np.argmin(np.linalg.norm(original.positions-original.cell.lengths()/2,axis=1)))
    rel=relative_positions(original,center);radii=np.linalg.norm(rel,axis=1)
    order=np.argsort(radii);first=order[1:13];first_radius=float(np.median(radii[first]))
    cutoff=float((radii[order[12]]+radii[order[13]])/2)
    assert np.ptp(radii[first])<1e-10 and len(original)==2048
    v=rel[first[0]];axis=np.cross(v,[0.,0.,1.])
    if np.linalg.norm(axis)<1e-10:axis=np.cross(v,[0.,1.,0.])
    axis/=np.linalg.norm(axis)
    opposite=int(first[np.argmin(rel[first]@v)])
    cases=[dict(case_id='fcc_baseline',kind='baseline',amplitude=0.,unit='none',atoms=original)]
    for strain in [-.02,-.01,-.005,.005,.01,.02]:
        at=original.copy();at.set_cell(original.cell.array*(1+strain),scale_atoms=True)
        cases.append(dict(case_id=f'strain_{strain:+g}',kind='strain',amplitude=strain,unit='fraction',atoms=at))
    for angle in [-10,-5,-2,-1,-.5,.5,1,2,5,10]:
        theta=np.deg2rad(angle)
        rotated=v*np.cos(theta)+np.cross(axis,v)*np.sin(theta)+axis*np.dot(axis,v)*(1-np.cos(theta))
        at=original.copy();at.positions[first[0]]=original.positions[center]+rotated;at.wrap()
        cases.append(dict(case_id=f'angle_{angle:+g}',kind='angle',amplitude=angle,unit='degree',atoms=at))
    for number in [1,2]:
        remove=[int(first[0])] if number==1 else [int(first[0]),opposite]
        mask=np.ones(len(original),dtype=bool);mask[remove]=False
        cases.append(dict(case_id=f'vacancy_{number}',kind='vacancy',amplitude=number,unit='neighbors_removed',atoms=original[mask]))
    at=original.copy();at.positions += [1.137,-2.223,3.319];at.wrap()
    cases.append(dict(case_id='translation',kind='invariance',amplitude=0.,unit='translation',atoms=at))
    q,_=np.linalg.qr(np.random.default_rng(core.SEED+801).normal(size=(3,3)))
    if np.linalg.det(q)<0:q[:,0]*=-1
    at=original.copy();at.set_cell(original.cell.array@q.T,scale_atoms=False);at.positions=original.positions@q.T
    cases.append(dict(case_id='rotation',kind='invariance',amplitude=0.,unit='rotation',atoms=at))
    perm=np.random.default_rng(core.SEED+802).permutation(len(original))
    cases.append(dict(case_id='permutation',kind='invariance',amplitude=0.,unit='permutation',atoms=original[perm]))
    geometry=[];folder=ROOT/'fcc_structures';folder.mkdir(exist_ok=True)
    for case in cases:
        at=case['atoms'];source_ids=at.get_array('original_index')
        ci=int(np.flatnonzero(source_ids==center)[0]);case['center_index']=ci
        new_rel=relative_positions(at,ci);new_radius=np.linalg.norm(new_rel,axis=1)
        common=[int(i) for i in first if i in source_ids]
        indices=[int(np.flatnonzero(source_ids==i)[0]) for i in common]
        old_vectors=rel[common]/radii[common,None]
        new_vectors=new_rel[indices]/new_radius[indices,None]
        cosine_change=float(np.max(np.abs(old_vectors@old_vectors.T-new_vectors@new_vectors.T)))
        radial_change=float(np.max(np.abs(new_radius-radii[source_ids])))
        count=int(np.sum((new_radius>1e-8)&(new_radius<cutoff)))
        record=dict(case_id=case['case_id'],kind=case['kind'],amplitude=case['amplitude'],unit=case['unit'],
            atoms=len(at),center_original_index=center,center_index=ci,volume_A3=float(at.get_volume()),
            first_shell_count=count,first_shell_cutoff_A=cutoff,max_radial_change_A=radial_change,
            common_first_shell_cosine_change=cosine_change)
        if case['kind']=='angle':assert radial_change<1e-10 and count==12
        if case['kind']=='strain':assert cosine_change<1e-10 and count==12
        if case['kind']=='vacancy':assert count==12-case['amplitude'] and radial_change<1e-10
        at.info.update(case_id=case['case_id'],kind=case['kind'],amplitude=case['amplitude'],center_original_index=center)
        write(folder/(case['case_id'].replace('+','p').replace('-','m')+'.xyz'),at,format='extxyz')
        geometry.append(record)
    pd.DataFrame(geometry).to_csv(ROOT/'fcc_geometry.csv',index=False)
    core.save_json(ROOT/'fcc_reference.json',dict(lattice_constant_A=a,per_state_estimates_A=by_state,
        method='sqrt(2) times median of first 12 neighbor distances of original deep-bulk reference atoms',
        perfect_fcc_atoms=2048,repeat=[8,8,8],center_original_index=center,
        nearest_neighbor_A=first_radius,first_shell_cutoff_A=cutoff,cases=len(cases)))
    return cases

def probes():
    cases=make_fcc_cases();rows=[];failures=[]
    folder=ROOT/'probe_vectors';folder.mkdir(exist_ok=True)
    for rbf,c in itertools.product(['gto','polynomial'],configurations()):
        cid=c['config_id'];caught=[]
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always');soap=make_soap(c,rbf)
                x=np.vstack([soap.create(case['atoms'],centers=[case['center_index']],n_jobs=1) for case in cases])
            core.unit(x)
            np.savez(folder/f'{rbf}_{cid}.npz',vectors=x,case_ids=[case['case_id'] for case in cases])
            for weighting in ['unique_terms','symmetric_frobenius']:
                xx=x if weighting=='unique_terms' else x*diagonal_weights(c)
                lengths=np.linalg.norm(xx,axis=1);uu=core.unit(xx)
                for metric,values in [('raw',xx),('norm',uu)]:
                    distance=np.linalg.norm(values-values[0],axis=1)
                    relative=distance/max(np.linalg.norm(values[0]),1e-300)
                    vacancy_index=next(i for i,case in enumerate(cases) if case['case_id']=='vacancy_1')
                    for i,case in enumerate(cases):
                        rows.append(dict(rbf=rbf,config_id=cid,r_cut=c['r_cut'],n_max=c['n_max'],
                            weighting=weighting,metric=metric,case_id=case['case_id'],kind=case['kind'],
                            amplitude=case['amplitude'],unit=case['unit'],distance=float(distance[i]),
                            distance_over_baseline_vector_norm=float(relative[i]),
                            distance_over_single_vacancy=float(distance[i]/max(distance[vacancy_index],1e-300)),
                            soap_norm=float(lengths[i]),norm_ratio=float(lengths[i]/lengths[0])))
            log('fcc_probes_progress',rbf=rbf,config_id=cid,cases=len(cases))
        except Exception as exc:
            failures.append(dict(rbf=rbf,config_id=cid,error=repr(exc),warnings=sorted(set(str(w.message) for w in caught))))
            log('fcc_probes_failure',rbf=rbf,config_id=cid,error=repr(exc))
    pd.DataFrame(rows).to_csv(ROOT/'fcc_responses.csv',index=False)
    core.save_json(ROOT/'fcc_failures.json',failures)
    log('fcc_probes_complete',rows=len(rows),failed_configs=len(failures))

def induced_metrics(x,ref,pairs):
    result={}
    for metric,xx,rr in [('raw',x,ref),('norm',core.unit(x),core.unit(ref))]:
        d=cdist(xx,xx);result['pair_'+metric]=d[pairs['i'],pairs['j']]
        np.fill_diagonal(d,np.inf);result['nn_'+metric]=np.argpartition(d,9,axis=1)[:,:10]
        dr=cdist(xx,rr);result['D10_'+metric]=np.partition(dr,9,axis=1)[:,:10].mean(axis=1)
    return result

def compare(a,b,metadata,pairs):
    aggregate=[];details=[]
    for metric in ['raw','norm']:
        aa=a['pair_'+metric];bb=b['pair_'+metric];corrs=[];changes=[]
        for ci,name in enumerate(pairs['names']):
            mask=pairs['category']==ci;x=aa[mask];y=bb[mask]
            corr=core.rho(x,y);change=float(np.median(np.abs(x-y))/max(np.median(y),1e-12))
            corrs.append(corr);changes.append(change)
            details.append(dict(**metadata,metric=metric,category=name,spearman=corr,median_relative_change=change))
        near=float(core.overlap(a['nn_'+metric],b['nn_'+metric]).mean())
        row=dict(**metadata,metric=metric,spearman_min=min(corrs),relative_change_max=max(changes),
            nn_overlap_mean=near,D10_spearman=core.rho(a['D10_'+metric],b['D10_'+metric]))
        row['passed']=row['spearman_min']>=.98 and row['relative_change_max']<=.05 and near>=.8
        aggregate.append(row)
    return aggregate,details

def analyze():
    pairs=dict(np.load(MAIN/'fixed_pairs.npz'));cache={};rows=[];details=[];norm_summaries=[]
    folder=ROOT/'comparison_metrics';folder.mkdir(exist_ok=True)
    for rbf,c in itertools.product(['gto','polynomial'],configurations()):
        try:v=vectors(rbf,c)
        except FileNotFoundError:continue
        for weighting in ['unique_terms','symmetric_frobenius']:
            x=v['evaluation'];ref=v['reference']
            if weighting=='symmetric_frobenius':
                w=diagonal_weights(c);x=x*w;ref=ref*w
            m=induced_metrics(x,ref,pairs)
            cache[rbf,c['config_id'],weighting]=m
            np.savez(folder/f"{rbf}_{c['config_id']}_{weighting}.npz",**m)
        norm_summaries.append(dict(rbf=rbf,config_id=c['config_id'],median_norm=float(np.median(np.linalg.norm(v['evaluation'],axis=1)))))
    for rbf,rc,weighting in itertools.product(['gto','polynomial'],[4.,5.],['unique_terms','symmetric_frobenius']):
        for n1,n2 in [(6,8),(8,10)]:
            lo=core.config_id(rc,.3,n1,10);hi=core.config_id(rc,.3,n2,10)
            if (rbf,lo,weighting) not in cache or (rbf,hi,weighting) not in cache:continue
            meta=dict(comparison='within_basis_radial',rbf=rbf,r_cut=rc,weighting=weighting,lower=lo,upper=hi)
            a,b=compare(cache[rbf,lo,weighting],cache[rbf,hi,weighting],meta,pairs)
            rows.extend(a);details.extend(b)
    for c,weighting in itertools.product(configurations(),['unique_terms','symmetric_frobenius']):
        cid=c['config_id']
        if ('polynomial',cid,weighting) not in cache:continue
        meta=dict(comparison='cross_basis',rbf='gto_to_polynomial',r_cut=c['r_cut'],weighting=weighting,lower='gto_'+cid,upper='polynomial_'+cid)
        a,b=compare(cache['gto',cid,weighting],cache['polynomial',cid,weighting],meta,pairs)
        rows.extend(a);details.extend(b)
    pd.DataFrame(rows).to_csv(ROOT/'distance_comparisons.csv',index=False)
    pd.DataFrame(details).to_csv(ROOT/'distance_comparisons_by_category.csv',index=False)
    pd.DataFrame(norm_summaries).to_csv(ROOT/'norm_scales.csv',index=False)
    p=pd.read_csv(ROOT/'fcc_responses.csv');decomp=[]
    keys=['rbf','config_id','r_cut','n_max','weighting','case_id','kind','amplitude','unit']
    raw=p[p.metric=='raw'];normalized=p[p.metric=='norm']
    joined=raw.merge(normalized,on=keys,suffixes=('_raw','_norm'),validate='one_to_one')
    for record in joined.to_dict('records'):
        t=record['norm_ratio_raw'];dr=record['distance_over_baseline_vector_norm_raw'];dn=record['distance_norm']
        identity=dr*dr-((t-1)**2+t*dn*dn)
        valid=dr>1e-8
        decomp.append({key:record[key] for key in keys}|dict(norm_ratio=t,relative_raw=dr,d_norm=dn,
            magnitude_fraction=(t-1)**2/(dr*dr) if valid else np.nan,
            direction_fraction=t*dn*dn/(dr*dr) if valid else np.nan,identity_residual=identity))
    pd.DataFrame(decomp).to_csv(ROOT/'fcc_decomposition.csv',index=False)
    log('analysis_complete',comparisons=len(rows),category_rows=len(details),fcc_decompositions=len(decomp))

def rotation_audit():
    """Localize the failed scientific acceptance check; do not repair or replace its data."""
    from ase import Atoms
    from ase.io import read
    import urllib.request
    reference=json.loads((ROOT/'fcc_reference.json').read_text())
    original=read(ROOT/'fcc_structures/fcc_baseline.xyz')
    center=int(np.flatnonzero(original.arrays['original_index']==reference['center_original_index'])[0])
    rel=relative_positions(original,center);radius=np.linalg.norm(rel,axis=1)
    indices=np.flatnonzero(radius<7.)
    indices=np.r_[center,indices[indices!=center]]
    cluster=Atoms('Ni'*len(indices),positions=rel[indices],pbc=False)
    rows=[];block_rows=[]
    for rbf,rc,lmax,periodic in itertools.product(['gto','polynomial'],[4.,5.],[0,4,6,8,9,10],[False,True]):
        at=original if periodic else cluster;ci=center if periodic else 0
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            soap=SOAP(species=['Ni'],periodic=periodic,r_cut=rc,sigma=.3,n_max=8,l_max=lmax,
                rbf=rbf,average='off',dtype='float64')
            baseline=soap.create(at,centers=[ci])[0]
            for seed in [core.SEED+801,core.SEED+803,core.SEED+804]:
                q,_=np.linalg.qr(np.random.default_rng(seed).normal(size=(3,3)))
                if np.linalg.det(q)<0:q[:,0]*=-1
                changed=at.copy();changed.positions=at.positions@q.T
                if periodic:changed.set_cell(at.cell.array@q.T,scale_atoms=False)
                rotated=soap.create(changed,centers=[ci])[0]
                relative=float(np.linalg.norm(rotated-baseline)/np.linalg.norm(baseline))
                normalized=float(np.linalg.norm(rotated/np.linalg.norm(rotated)-baseline/np.linalg.norm(baseline)))
                rows.append(dict(rbf=rbf,r_cut=rc,n_max=8,l_max=lmax,periodic=periodic,seed=seed,
                    relative_raw=relative,d_norm=normalized,passed=relative<1e-8 and normalized<1e-8))
                if lmax==10:
                    for l in range(11):
                        a=baseline[l*36:(l+1)*36];b=rotated[l*36:(l+1)*36]
                        block_rows.append(dict(rbf=rbf,r_cut=rc,periodic=periodic,seed=seed,l=l,
                            block_difference_over_full_norm=float(np.linalg.norm(b-a)/np.linalg.norm(baseline))))
    pd.DataFrame(rows).to_csv(ROOT/'rotation_audit.csv',index=False)
    pd.DataFrame(block_rows).to_csv(ROOT/'rotation_audit_by_l.csv',index=False)
    # Preserve a source lead without asserting that it proves the compiled binary's root cause.
    url='https://raw.githubusercontent.com/SINGROUP/dscribe/v2.1.2/dscribe/ext/soapGeneral.cpp'
    try:
        content=urllib.request.urlopen(url,timeout=30).read()
        path=ROOT/'implementation'/'upstream_v2.1.2_soapGeneral.cpp';path.write_bytes(content)
        core.save_json(ROOT/'implementation'/'source_lead.json',dict(url=url,sha256=core.sha(path),
            observation='factorListSet contains c[65]=0 at (l,m)=(10,10); factorY indexes l*(l+1)/2+m and getYlmi uses it.',
            status='mechanism lead, not a patched binary or a proven complete root cause'))
    except Exception as exc:core.save_json(ROOT/'implementation'/'source_lead.json',dict(url=url,error=repr(exc)))
    log('rotation_audit_complete',rows=len(rows),failed=sum(not r['passed'] for r in rows))

def save_figure(fig,name):
    import matplotlib.pyplot as plt
    folder=ROOT/'figures';folder.mkdir(exist_ok=True)
    fig.savefig(folder/(name+'.png'),dpi=170,bbox_inches='tight')
    fig.savefig(folder/(name+'.pdf'),bbox_inches='tight')
    plt.close(fig)

def figures():
    import matplotlib.pyplot as plt
    from ase.io import read
    poly=pd.read_csv(ROOT/'polynomial_summary.csv');health=[]
    for c in configurations():
        record=json.loads((MAIN/'records'/f"{c['config_id']}.json").read_text(encoding='utf-8'))
        health.append(dict(rbf='gto',r_cut=c['r_cut'],n_max=c['n_max'],residual=record['basis_residual_max'],condition=record['basis_condition_max']))
    for row in poly.to_dict('records'):
        health.append(dict(rbf='polynomial',r_cut=row['r_cut'],n_max=row['n_max'],residual=row['basis_residual_max'],condition=row['basis_condition_max']))
    h=pd.DataFrame(health);h.to_csv(ROOT/'basis_health_comparison.csv',index=False)
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for (rbf,rc),part in h.groupby(['rbf','r_cut']):
        color='#0072B2' if rbf=='gto' else '#D55E00';style='-' if rc==4 else '--'
        part=part.sort_values('n_max');label=f'{rbf}, rc={rc:g} A'
        axes[0].semilogy(part.n_max,part.condition,'o'+style,color=color,label=label)
        axes[1].semilogy(part.n_max,part.residual,'o'+style,color=color,label=label)
    axes[0].set_ylabel('max overlap matrix condition number');axes[1].set_ylabel('max radial orthogonality residual')
    for ax in axes:ax.set_xlabel('n_max');ax.set_xticks([6,8,10]);ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Basis health: no warning does not establish accurate orthogonalization')
    fig.tight_layout();save_figure(fig,'basis_health')
    d=pd.read_csv(ROOT/'distance_comparisons.csv')
    sub=d[(d.comparison=='within_basis_radial')&(d.weighting=='unique_terms')]
    fig,axes=plt.subplots(2,3,figsize=(14,8))
    for row,rc in enumerate([4.,5.]):
        for (rbf,metric),part in sub[sub.r_cut==rc].groupby(['rbf','metric']):
            color='#0072B2' if rbf=='gto' else '#D55E00';style='-' if metric=='raw' else '--'
            part=part.sort_values('lower');x=np.arange(len(part));label=rbf+' '+metric
            for col,field in enumerate(['spearman_min','relative_change_max','nn_overlap_mean']):
                axes[row,col].plot(x,part[field],'o'+style,color=color,label=label)
        for col,(label,threshold) in enumerate([('minimum stratified distance rho',.98),('maximum median relative change',.05),('mean top-10 neighbor overlap',.8)]):
            ax=axes[row,col];ax.axhline(threshold,color='grey',ls=':',lw=1)
            ax.set_xticks([0,1],['6 to 8','8 to 10']);ax.set_xlabel('n_max change');ax.set_ylabel(label)
            ax.set_title(f'rc={rc:g} A, sigma=0.3 A, l_max=10');ax.grid(alpha=.15)
    axes[0,0].legend(fontsize=8);fig.suptitle('Original two distances: the three tests must pass together')
    fig.tight_layout();save_figure(fig,'radial_convergence')
    responses=pd.read_csv(ROOT/'fcc_responses.csv')
    sub=responses[(responses.weighting=='unique_terms')&(responses.n_max==8)]
    for norm_only,name in [(False,'fcc_distance_responses'),(True,'fcc_norm_changes')]:
        fig,axes=plt.subplots(2,3,figsize=(14,8))
        for row,rc in enumerate([4.,5.]):
            for col,(kind,xlabel,multiplier) in enumerate([('strain','uniform strain (%)',100),('angle','one-neighbor rotation (degree)',1),('vacancy','first-shell neighbors removed',1)]):
                p=sub[(sub.r_cut==rc)&(sub.kind==kind)]
                for (rbf,metric),part in p.groupby(['rbf','metric']):
                    if norm_only and metric!='raw':continue
                    color='#0072B2' if rbf=='gto' else '#D55E00';style='-' if metric=='raw' else '--'
                    xs=np.r_[part.amplitude.to_numpy()*multiplier,0.]
                    ys=(part.norm_ratio.to_numpy()-1)*100 if norm_only else part.distance_over_baseline_vector_norm.to_numpy()
                    ys=np.r_[ys,0.];order=np.argsort(xs)
                    axes[row,col].plot(xs[order],ys[order],'o'+style,color=color,ms=3,label=rbf if norm_only else rbf+' '+metric)
                ax=axes[row,col];ax.set_xlabel(xlabel);ax.set_title(f'rc={rc:g} A, n_max=8');ax.grid(alpha=.15)
                ax.set_ylabel('SOAP norm change (%)' if norm_only else 'dimensionless distance response')
        axes[0,0].legend(fontsize=8)
        fig.suptitle('Norm changes under controlled geometry' if norm_only else 'Raw distance / original norm versus normalized distance (distinct definitions)')
        fig.tight_layout();save_figure(fig,name)
    deco=pd.read_csv(ROOT/'fcc_decomposition.csv')
    fig,axes=plt.subplots(2,2,figsize=(11,8))
    for ax,(rbf,rc) in zip(axes.flat,itertools.product(['gto','polynomial'],[4.,5.])):
        p=deco[(deco.rbf==rbf)&(deco.r_cut==rc)&(deco.n_max==8)&(deco.weighting=='unique_terms')]
        names=['strain_+0.01','angle_+10','vacancy_1'];p=p.set_index('case_id').loc[names]
        mag=p.magnitude_fraction.to_numpy();direction=p.direction_fraction.to_numpy()
        ax.bar(np.arange(3),mag,color='#0072B2',label='norm-magnitude term')
        ax.bar(np.arange(3),direction,bottom=mag,color='#E69F00',label='normalized-direction term')
        for i,v in enumerate(mag):ax.text(i,.5*v if v>=.08 else v+.04,f'{100*v:.1f}%',ha='center',va='center',fontsize=9,color='white' if v>=.08 else 'black')
        ax.set_ylim(0,1.04);ax.set_xticks(range(3),['+1% strain','10 degree angle','one vacancy'])
        ax.set_ylabel('fraction of raw squared distance');ax.set_title(f'{rbf}, rc={rc:g} A, n_max=8')
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,fontsize=8,loc='upper center',bbox_to_anchor=(.5,.94),ncol=2)
    fig.suptitle('Exact decomposition; polynomial l=10 shown as implementation diagnostics')
    fig.tight_layout(rect=(0,0,1,.91));save_figure(fig,'fcc_distance_decomposition')
    fig=plt.figure(figsize=(14,4))
    reference=json.loads((ROOT/'fcc_reference.json').read_text());cutoff=reference['first_shell_cutoff_A']
    files=['fcc_baseline.xyz','strain_p0.01.xyz','angle_p10.xyz','vacancy_1.xyz']
    titles=['ideal FCC','+1% radial stretch','one neighbor rotated 10 degrees','one first-shell neighbor removed']
    for i,(file,title) in enumerate(zip(files,titles)):
        at=read(ROOT/'fcc_structures'/file);center=int(np.flatnonzero(at.arrays['original_index']==reference['center_original_index'])[0])
        rel=relative_positions(at,center);r=np.linalg.norm(rel,axis=1);mask=(r>1e-8)&(r<cutoff)
        ax=fig.add_subplot(1,4,i+1,projection='3d')
        ax.scatter(*rel[mask].T,s=35,color='#0072B2');ax.scatter([0],[0],[0],color='#D55E00',marker='*',s=100)
        for v in rel[mask]:ax.plot([0,v[0]],[0,v[1]],[0,v[2]],color='grey',lw=.6)
        ax.set_title(title,fontsize=8);ax.set_xlabel('x (A)',fontsize=8);ax.set_ylabel('y (A)',fontsize=8);ax.set_zlabel('z (A)',fontsize=8)
        ax.set_xlim(-2.7,2.7);ax.set_ylim(-2.7,2.7);ax.set_zlim(-2.7,2.7);ax.set_box_aspect([1,1,1]);ax.view_init(20,35)
    fig.suptitle('Controlled local geometry; red star = common central atom')
    fig.tight_layout();save_figure(fig,'fcc_geometry')

def table(frame,columns,labels=None):
    labels=labels or columns
    lines=['| '+' | '.join(labels)+' |','| '+' | '.join(['---']*len(columns))+' |']
    for row in frame.to_dict('records'):
        text=[]
        for key in columns:
            value=row[key]
            if isinstance(value,(bool,np.bool_)):value='通过' if value else '未通过'
            elif isinstance(value,float):value=f'{value:.5g}'
            text.append(str(value))
        lines.append('| '+' | '.join(text)+' |')
    return '\n'.join(lines)

def report():
    figures()
    poly=pd.read_csv(ROOT/'polynomial_summary.csv');d=pd.read_csv(ROOT/'distance_comparisons.csv')
    main_edges=d[(d.comparison=='within_basis_radial')&(d.weighting=='unique_terms')]
    weighted_edges=d[(d.comparison=='within_basis_radial')&(d.weighting=='symmetric_frobenius')]
    deco=pd.read_csv(ROOT/'fcc_decomposition.csv');reference=json.loads((ROOT/'fcc_reference.json').read_text())
    example=deco[(deco.rbf=='gto')&(deco.config_id=='r4_s0p3_n8_l10')&(deco.weighting=='unique_terms')&deco.case_id.isin(['strain_+0.01','angle_+10','vacancy_1'])].copy()
    example['norm_change_percent']=(example.norm_ratio-1)*100
    example['magnitude_percent']=example.magnitude_fraction*100
    protocol=json.loads((ROOT/'protocol.json').read_text())
    lines=['# Ni SOAP 径向基函数与受控结构对照', '',f'生成时间：{core.stamp()}。这是上一轮的独立补充，原 240 组数据、协议和判据保持原样。',
        '', '## 主要结论', '',
        f'- 六组 polynomial 全部生成有限向量，均没有警告；但 n_max=10 的实际径向正交化残差为 {poly[poly.n_max==10].basis_residual_max.min():.4g}—{poly[poly.n_max==10].basis_residual_max.max():.4g}。无警告不能证明正交化准确。',
        '- polynomial 在本轮多个相邻阶数比较中提高了近邻重合率，尤其是 8→10；仍未通过全部预先设定的收敛条件，不能直接替换 GTO 后冻结参数。',
        '- 夹角不变的均匀伸缩仍改变 d_norm；中心到各邻居半径不变的角度扰动仍改变范数。norm 含有径向信息，SOAP 范数不是密度的专用指标。',
        '- 两种距离均对邻居缺失有响应。范数项在 raw 平方距离中的占比随扰动、观察范围和基函数改变，不能直接套用一个固定密度占比。',
        '', '## 固定的试验设置', '',
        '- 新计算：r_cut=4/5 Å、sigma=0.3 Å、l_max=10、n_max=6/8/10，共六组 polynomial；GTO 对照复用原六组向量。',
        '- 真实 Ni 对照继续使用原 1000 初筛中心、400 个不重叠的深体相参考、同一批 20000 分层原子对。没有重新抽样，没有能量过滤。',
        '- 主分析仍是 DScribe 唯一项向量的 d_raw 和归一化后的 d_norm。对称项 √2 权重仅作独立诊断，未覆盖主距离或原选择结果。',
        '- 比较不同基函数诱导的距离、近邻和 D(k)，不直接比较其系数坐标。相同 n_max 不代表相同逼近精度；跨基函数差别也不应全部归结为数值误差。',
        '- 原相邻边标准保持不变：十类原子对最低 Spearman >=0.98，最大中位相对距离变化 <=5%，平均前十近邻重合率 >=80%。',
        '', '## 径向数值诊断', '',
        'polynomial 的诊断直接积分 DScribe 实际生成的径向函数：使用与实现匹配的 100 点 Gauss-Legendre 节点，计算 ∫g_n(r)g_m(r)r²dr 与单位矩阵的无穷范数差。节点一致性已经检查。GTO 对照沿用原重叠矩阵诊断。条件数依赖基函数缩放，应结合正交化残差理解。', '',
        table(poly,['config_id','basis_condition_max','basis_residual_max','warning_count'],['配置','重叠矩阵条件数','径向正交化残差','警告数']), '',
        '![径向基函数数值表现](figures/basis_health.png)',
        '', 'polynomial 的 n_max=6 残差约 2×10⁻⁸，8 为 10⁻⁵ 量级，10 增至约 0.03。增加阶数时逼近能力与数值误差同时变化，因此不能把最高阶结果当作真值。',
        '', '## 真实 Ni 环境的径向比较', '',
        table(main_edges,['rbf','r_cut','lower','upper','metric','spearman_min','relative_change_max','nn_overlap_mean','passed'],
            ['基函数','rc','低阶','高阶','距离','最低ρ','最大相对变化','近邻重合率','全部条件']), '',
        f'主距离相邻边 {int(main_edges.passed.sum())}/{len(main_edges)} 通过；对称权重诊断边 {int(weighted_edges.passed.sum())}/{len(weighted_edges)} 通过。', '',
        '![径向距离收敛](figures/radial_convergence.png)', '',
        'polynomial 8→10 的平均近邻重合率约 0.84—0.86，高于本轮 GTO 对应的约 0.65—0.77。但其 n_max=10 正交化残差也明显增大，这只能说明该分支值得做稳定化后的数值对照，不能据此认定已经收敛。组内细节在 `distance_comparisons_by_category.csv`。',
        '', '跨基函数的全部 24 个距离比较记录在 `distance_comparisons.csv`，其 comparison=cross_basis 部分是表示差异诊断，表中通过标记只表示满足同一数值参考线，不是物理等价检验。',
        '', '## FCC 受控几何', '',
        f'从原深体相参考原子的前 12 近邻距离估计 FCC 晶格常数 a={reference["lattice_constant_A"]:.9f} Å；理想第一壳半径 {reference["nearest_neighbor_A"]:.9f} Å。建立 8×8×8 常规 FCC 超胞，共 2048 个原子，跟踪同一个中心原子的原始身份。', '',
        '22 个构型包含：基准；±0.5%、±1%、±2% 均匀伸缩；一个第一壳邻居绕中心转动 ±0.5°、±1°、±2°、±5°、±10°；移除一个或一对相反方向第一壳邻居；整体平移、旋转和原子置换。', '',
        '- 角度试验检查所有剩余原子到中心的半径不变、第一壳仍为 12。',
        '- 伸缩试验检查第一壳邻居间夹角不变。',
        '- 缺邻居试验检查第一壳配位数依次为 11、10，剩余原子半径不变。',
        '- 平移、旋转、置换用于数值不变性检查。合成几何不代表温度、能垒或真实缺陷平衡态。', '',
        '![受控局域几何](figures/fcc_geometry.png)', '',
        '所有几何验证值见 `fcc_geometry.csv`；完整周期结构在 `fcc_structures/`，可用结构软件查看。缺邻居同时改变配位和相对密度；均匀伸缩同时改变距离和密度，不把它们称为单一信息通道。',
        '', '## 两种距离到底响应什么', '',
        '为展示响应，将 raw 距离除以原始向量范数，记为 relative_raw；norm 距离仍是单位向量间欧氏距离。它们均无量纲，但定义不同，不按谁的数值更大判优。', '',
        '![受控结构的距离响应](figures/fcc_distance_responses.png)', '',
        '![受控结构的范数变化](figures/fcc_norm_changes.png)', '',
        '以数值状况较好的 GTO、rc=4 Å、sigma=0.3 Å、n_max=8、l_max=10 为例：', '',
        table(example,['case_id','norm_change_percent','relative_raw','d_norm','magnitude_percent'],
              ['几何变化','范数变化%','raw/原始范数','d_norm','范数项占raw平方距离%']), '',
        '这些例子直接说明：归一化并没有消除径向结构变化；向量范数也会随纯角度重排变化。删除邻居的范数项较大，是该配置该探针的结果，不代表实际晶界中范数项的唯一物理来源。', '',
        '精确分解为 relative_raw²=(norm_ratio−1)²+norm_ratio·d_norm²，其中 norm_ratio=||p_new||/||p_original||。数值恒等式单独验证；极小的不变性距离不解释其分解百分比。', '',
        '![平方距离分解](figures/fcc_distance_decomposition.png)',
        '', '## 建议的后续处置', '',
        '1. 继续并行报告 d_raw、d_norm 和范数，避免把它们预先贴成“密度距离”与“纯角度距离”。',
        '2. 优先诊断 polynomial 的高阶正交化；若实现稳定化对照，应单独记录实现、精度和基函数变化，保留这次 DScribe 原生结果。稳定化之前不扩大阶数。',
        '3. 同时保留“表示收敛”与“科学结论稳健”两类证据。前者本轮未达到原标准，后者需要针对具体结构比较问题判断，不能事后放宽原标准制造通过结果。',
        '4. 受控结构可帮助解释响应，但未校准可接受的热扰动或应变范围，也未提供独立晶界重复，因此本轮不冻结唯一参数或唯一主距离。',
        '', '## 复现与文件', '',
        '- `protocol.json`：执行前固定的设置、主分析约定、六组 GTO 来源文件哈希。',
        '- `polynomial_summary.csv`、`basis_health_comparison.csv`：数值诊断。',
        '- `distance_comparisons.csv`、`distance_comparisons_by_category.csv`：全部距离比较和十类原子对细节。',
        '- `fcc_responses.csv`、`fcc_decomposition.csv`：两种距离、范数和分解。',
        '- `vectors/`、`metrics/`、`probe_vectors/`：原始向量及派生数据。',
        '- `verification.json`：数值、几何、原始数据未改变等检查。',
        '- `source/`：运行脚本、依赖脚本快照；`implementation/` 保留实际安装版本的 polynomial 基函数函数。',
        '', '从工作区依次执行 `work/ni_soap_sweep_20261101/radial_followup.py prepare`、`scan`、`probes`、`analyze`、`rotation_audit`、`report`、`verify`；已有逐配置结果可复用。',
        '', f'原主扫描协议 SHA256：{next(x["sha256"] for x in protocol["main_input_evidence"] if x["path"].endswith("protocol.json"))}。',
        '', '算法参考：[DScribe SOAP 文档](https://singroup.github.io/dscribe/latest/tutorials/descriptors/soap.html)。实际实现已经保存在本地，所有结论以本次环境的输出为依据。']
    audit=pd.read_csv(ROOT/'rotation_audit.csv')
    inv=pd.read_csv(ROOT/'fcc_responses.csv')
    inv=inv[(inv.kind=='invariance')&(inv.weighting=='unique_terms')]
    inv_by_config=inv.groupby(['rbf','config_id']).distance_over_baseline_vector_norm.max().reset_index()
    strongest=float(inv_by_config[inv_by_config.rbf=='polynomial'].distance_over_baseline_vector_norm.max())
    lines.insert(lines.index('## 主要结论')+2,
        f'- **物理验收未通过：原生 polynomial、l_max=10 在整体旋转后出现最高约 {strongest:.3%} 的相对向量变化。其近邻稳定性提高不能作为物理准确性提高的证据。GTO 对照通过不变性检查。**')
    audit_table=audit.groupby(['rbf','r_cut','l_max']).relative_raw.max().reset_index()
    audit_table=audit_table[audit_table.l_max.isin([9,10])].pivot(index=['r_cut','l_max'],columns='rbf',values='relative_raw').reset_index()
    section=['', '## 验收失败：polynomial 的 l=10 旋转偏差', '',
        '使用原来 1e-8 的相对不变性容差进行验收时，发现六组原生 polynomial、l_max=10 均不能通过整体旋转检查。平移和同元素置换保持数值稳定，GTO 的整体旋转也稳定。该结果作为科学验收失败保留，未放宽阈值。', '',
        '补充定位固定 n_max=8，分别测试 r_cut=4/5 Å，l_max=0/4/6/8/9/10，三种整体旋转，并同时使用周期 FCC 超胞和非周期 7 Å 局域簇，共 144 次比较。12 个失败全部属于 polynomial、l_max=10；低阶和 GTO 在这些旋转下通过。非周期簇重复了周期结构的偏差，排除了“仅由周期边界引起”的解释。', '',
        table(audit_table,['r_cut','l_max','gto','polynomial'],['rc','最高角阶','GTO 最大相对旋转变化','polynomial 最大相对旋转变化']), '',
        '逐 l 通道诊断将偏差集中到 l=10，l≤9 通道保持数值稳定。公开 v2.1.2 `soapGeneral.cpp` 中，factorY(l,m) 使用 l(l+1)/2+m 索引，factorListSet 对 (l,m)=(10,10) 对应的 c[65] 给出 0，并由 getYlmi 使用。这是与观测吻合的实现线索；本轮未修改或重编译依赖库，尚未用修补前后对照证明完整根因。源码与哈希保存在 `implementation/`。', '',
        '证据：`rotation_audit.csv`、`rotation_audit_by_l.csv`、`fcc_responses.csv`。报告里的 polynomial 原生数据用于记录表示与实现诊断；物理响应解释优先采用已经通过不变性检查的 GTO 对照。低阶旋转测试通过，也不等于其径向展开或角向精度已经收敛。', '']
    position=lines.index('## 真实 Ni 环境的径向比较');lines[position:position]=section
    old='2. 优先诊断 polynomial 的高阶正交化；若实现稳定化对照，应单独记录实现、精度和基函数变化，保留这次 DScribe 原生结果。稳定化之前不扩大阶数。'
    if old in lines:lines[lines.index(old)]='2. 先处理原生 polynomial 在 l=10 的旋转偏差，再诊断高阶径向正交化。修改实现或限制角阶均需单独对照和重新验证；当前 l=10 数据不用于冻结物理距离，暂不扩大阶数。'
    (ROOT/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    (ROOT/'README.md').write_text('# SOAP 径向与结构对照\n\n阅读 [REPORT.md](REPORT.md)。\n\n六组 polynomial 与复用的六组 GTO，22 个 FCC 控制构型，两种距离均保留。\n',encoding='utf-8')
    source=ROOT/'source';source.mkdir(exist_ok=True)
    for file in ['radial_followup.py','run_sweep.py']:
        shutil.copy2(core.SOURCE/file,source/file)
    log('report_generated',report=str(ROOT/'REPORT.md'),figures=len(list((ROOT/'figures').glob('*.png'))))

def verify():
    protocol=json.loads((ROOT/'protocol.json').read_text())
    for item in protocol['main_input_evidence']:assert core.sha(item['path'])==item['sha256']
    s=pd.read_csv(ROOT/'polynomial_summary.csv');assert len(s)==6 and (s.status=='ok').all()
    for c in configurations():
        v=vectors('polynomial',c);dimension=c['n_max']*(c['n_max']+1)//2*(c['l_max']+1)
        assert v['evaluation'].shape==(1000,dimension) and v['reference'].shape==(400,dimension)
        core.unit(v['evaluation']);core.unit(v['reference'])
        saved=dict(np.load(ROOT/'metrics'/f"polynomial_{c['config_id']}.npz"))
        repeat=dict(np.load(ROOT/'comparison_metrics'/f"polynomial_{c['config_id']}_unique_terms.npz"))
        for metric in ['raw','norm']:
            np.testing.assert_allclose(saved['pair_'+metric],repeat['pair_'+metric],rtol=1e-13,atol=1e-13)
            np.testing.assert_allclose(saved['D10_'+metric],repeat['D10_'+metric],rtol=1e-13,atol=1e-13)
    p=pd.read_csv(ROOT/'fcc_responses.csv');assert len(p)==12*22*2*2
    assert not json.loads((ROOT/'fcc_failures.json').read_text())
    inv=p[p.kind=='invariance'];maximum=float(inv.distance_over_baseline_vector_norm.max())
    by_rbf=inv.groupby('rbf').distance_over_baseline_vector_norm.max().to_dict()
    assert by_rbf['gto']<1e-8
    decomp=pd.read_csv(ROOT/'fcc_decomposition.csv');max_identity=float(decomp.identity_residual.abs().max());assert max_identity<1e-12
    fractions=decomp[decomp.relative_raw>1e-8]
    np.testing.assert_allclose(fractions.magnitude_fraction+fractions.direction_fraction,1,rtol=1e-8,atol=1e-8)
    geometry=pd.read_csv(ROOT/'fcc_geometry.csv');assert len(geometry)==22
    assert geometry[geometry.kind=='angle'].max_radial_change_A.max()<1e-10
    assert geometry[geometry.kind=='strain'].common_first_shell_cosine_change.max()<1e-10
    assert list(geometry[geometry.kind=='vacancy'].first_shell_count)==[11,10]
    comparisons=pd.read_csv(ROOT/'distance_comparisons.csv');assert len(comparisons)==56
    details=pd.read_csv(ROOT/'distance_comparisons_by_category.csv');assert len(details)==560
    assert len(list((ROOT/'figures').glob('*.png')))==6
    import re
    for link in re.findall(r'!\[[^\]]*\]\(([^)]+)\)',(ROOT/'REPORT.md').read_text(encoding='utf-8')):assert (ROOT/link).exists()
    evidence=dict(verified_at=core.stamp(),main_inputs_unchanged=True,polynomial_success=6,reused_gto=6,
        fcc_cases=22,response_rows=len(p),distance_comparisons=56,category_rows=560,
        max_invariance_relative_distance=maximum,max_decomposition_identity_residual=max_identity,
        geometry_checks_passed=True,main_metric_recomputation_matches=True,
        invariance_max_by_rbf=by_rbf,artifact_consistency_checks_passed=True,
        scientific_invariance_all_passed=maximum<1e-8,all_checks_passed=maximum<1e-8,
        scientific_acceptance='failed_rotation_invariance_for_tested_native_polynomial_l10' if maximum>=1e-8 else 'passed_invariance')
    core.save_json(ROOT/'verification.json',evidence)
    source=ROOT/'source';shutil.copy2(__file__,source/Path(__file__).name)
    manifest=[dict(file=str(p.relative_to(ROOT)),sha256=core.sha(p),bytes=p.stat().st_size)
              for folder in [source,ROOT/'implementation'] for p in folder.glob('*.py')]
    core.save_json(ROOT/'source_manifest.json',manifest)
    log('verification_complete',**evidence)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','scan','probes','analyze','rotation_audit','report','verify'])
    args=parser.parse_args();globals()[args.stage]()

if __name__=='__main__':main()
