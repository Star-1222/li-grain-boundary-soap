import os
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'):os.environ[k]='1'
import sys,json,warnings,time,itertools
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation
from scipy.stats import spearmanr
from dscribe.descriptors import SOAP
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[1]/'tmp/soap_timing'))
from cfg_reader import read_cfg,get_auxiliary
manifest=json.loads((ROOT/'manifest.json').read_text());sel=np.load(ROOT/'selections.npz');centers=[sel[f'centers_{i}'] for i in range(2)];refs=[sel[f'references_{i}'] for i in range(2)];pops=[sel[f'population_{i}'] for i in range(2)]
atoms=[read_cfg(ROOT.parents[1]/f'tmp/soap_timing/ni_{s}.cfg',coordinate_convention='scaled') for s in ['stable','metastable']]
# Diagnostic representatives, deliberately not designated selected/optimal: no radial convergence edge passed.
parameters=[(3.25,.5,12,12),(4,.3,10,8),(5,.3,10,8),(6,.3,10,8),(6,.7,10,8)]
rng=np.random.default_rng(84120);probe=[]
for ci in centers:probe.append(np.concatenate([rng.choice(ci[:300],30,False),rng.choice(ci[300:],30,False)]))
allsummary=[]
for r,s,n,l in parameters:
 key=f'r{r:g}_s{s:g}_n{n}_l{l}';out=ROOT/'diagnostics'/key;out.mkdir(parents=True,exist_ok=True);t=time.perf_counter()
 if (out/'summary.json').exists():allsummary.append(json.loads((out/'summary.json').read_text()));continue
 with warnings.catch_warnings(record=True) as ws:
  warnings.simplefilter('always');soap=SOAP(species=['Ni'],periodic=True,r_cut=r,sigma=s,n_max=n,l_max=l,rbf='gto',average='off',dtype='float64')
  base=[soap.create(at,centers=ci,n_jobs=1) for at,ci in zip(atoms,probe)]
 x0=np.vstack(base);n0=np.linalg.norm(x0,axis=1);metrics=[]
 R=Rotation.from_rotvec(np.array([.37,.61,.23])).as_matrix()
 variants=[('translation',None,None),('rotation',None,None),('permutation',None,None)]
 variants += [('displacement',amp,seed) for amp,seed in itertools.product([.01,.03,.05],[7101,7102,7103])]
 variants += [('strain',amp,None) for amp in [-.01,-.005,.005,.01]]
 for kind,amp,seed in variants:
  xx=[]
  for st,(at,ci) in enumerate(zip(atoms,probe)):
   at=at.copy();ci=ci.copy()
   if kind=='translation':at.positions+=np.array([2.13,-1.37,4.29]);at.wrap()
   elif kind=='rotation':at.positions=at.positions@R.T;at.set_cell(at.cell.array@R.T)
   elif kind=='permutation':
    perm=np.random.default_rng(987+st).permutation(len(at));inv=np.argsort(perm);at=at[perm];ci=inv[ci]
   elif kind=='displacement':
    dx=np.random.default_rng(seed+st*100).normal(size=(len(at),3));dx*=amp/np.sqrt(np.mean(np.sum(dx*dx,axis=1)));at.positions+=dx;at.wrap()
   elif kind=='strain':at.set_cell(at.cell.array*(1+amp),scale_atoms=True)
   xx.append(soap.create(at,centers=ci,n_jobs=1))
  x=np.vstack(xx);nn=np.linalg.norm(x,axis=1)
  row={'kind':kind,'amplitude':amp,'seed':seed,'raw_median':float(np.median(np.linalg.norm(x-x0,axis=1))),'raw_max':float(np.max(np.linalg.norm(x-x0,axis=1))),'norm_median':float(np.median(np.linalg.norm(x/nn[:,None]-x0/n0[:,None],axis=1))),'norm_max':float(np.max(np.linalg.norm(x/nn[:,None]-x0/n0[:,None],axis=1)))}
  metrics.append(row)
 # Full 4659 original evaluation population, 400 disjoint references, both periodic interfaces.
 blocks=[soap.create(at,centers=np.concatenate([pop,ref]),n_jobs=1) for at,pop,ref in zip(atoms,pops,refs)]
 x=np.vstack([b[:len(pop)] for b,pop in zip(blocks,pops)]);rr=np.vstack([b[len(pop):] for b,pop in zip(blocks,pops)]);norm=np.linalg.norm(x,axis=1);nr=np.linalg.norm(rr,axis=1);data={'norm':norm}
 ref_sensitivity={}
 for kind,xx,yy in [('raw',x,rr),('norm',x/norm[:,None],rr/nr[:,None])]:
  ds=cdist(xx,yy);near=np.partition(ds,19,axis=1)[:,:20];near.sort(axis=1)
  for k in (5,10,20):data[f'{kind}_d{k}']=near[:,:k].mean(axis=1)
  for j,name in enumerate(['stable_ref','metastable_ref']):
   dd=np.partition(ds[:,j*200:(j+1)*200],9,axis=1)[:,:10].mean(axis=1);data[f'{kind}_{name}']=dd
   ref_sensitivity[f'{kind}_{name}_rho_vs_joint']=float(spearmanr(dd,data[f'{kind}_d10']).statistic)
 np.savez_compressed(out/'population.npz',**data)
 pop_info=[]
 for st,(at,pop) in enumerate(zip(atoms,pops)):
  ids=get_auxiliary(at,'id').astype(int);ps=set(centers[st]);region=[rr['region'] for rr in __import__('csv').DictReader(Path(r'C:/Users/1/Documents/Codex/2026-10-31/bang/work/li-grain-boundary-soap/soap/Stage02_SOAP表示与审计/03_输出/stage2_75_atom_metadata.csv').open(encoding='utf-8-sig')) if rr['structure']==['stable','metastable'][st]]
  for idx,reg in zip(pop,region):pop_info.append({'state':['stable','metastable'][st],'atom_id':int(ids[idx]),'region':reg,'z_A':float(at.positions[idx,2]),'z_scaled':float(at.get_scaled_positions(wrap=True)[idx,2]),'held_out':bool(idx not in ps)})
 (out/'population_metadata.json').write_text(json.dumps(pop_info))
 summary={'key':key,'role':'diagnostic representative; not convergence-qualified, not an optimal selection','r_cut':r,'sigma':s,'n_max':n,'l_max':l,'dimension':int(x.shape[1]),'probe_count':len(x0),'population_count':len(x),'heldout_count':sum(z['held_out'] for z in pop_info),'perturbation':metrics,'reference_sensitivity':ref_sensitivity,'warnings':sorted(set(str(w.message) for w in ws)),'total_s':time.perf_counter()-t}
 (out/'summary.json').write_text(json.dumps(summary,indent=2));allsummary.append(summary);print(key+' complete '+str(round(summary['total_s'],1))+'s',flush=True)
(ROOT/'diagnostics_summary.json').write_text(json.dumps(allsummary,indent=2));print('DIAGNOSTICS_DONE',flush=True)
