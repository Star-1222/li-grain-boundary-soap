import os
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[k]='1'
import sys,json,time,warnings,hashlib,itertools,csv
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr
from dscribe.descriptors import SOAP
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tmp'/'soap_timing'))
from cfg_reader import read_cfg,get_auxiliary
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT.parents[1]/'tmp'/'soap_timing'
META=Path(r'C:/Users/1/Documents/Codex/2026-10-31/bang/work/li-grain-boundary-soap/soap/Stage02_SOAP表示与审计/03_输出/stage2_75_atom_metadata.csv')
SEED=20261101
rng=np.random.default_rng(SEED)
labels=['stable','metastable']
atoms=[read_cfg(SOURCE/('ni_'+s+'.cfg'),coordinate_convention='scaled') for s in labels]
rows=list(csv.DictReader(META.open(encoding='utf-8-sig')))

# Round-robin sampling over fixed interface identity, spatial quartile and CNA/PTM.
def stratified(pool,at,n):
 z=at.get_scaled_positions(wrap=True)[pool,2]; gb=(np.abs(z-.5)<np.minimum(z,1-z)).astype(int)
 cna=get_auxiliary(at,'c_CNA').astype(int)[pool]; ptm=get_auxiliary(at,'c_PTM_1').astype(int)[pool]
 # Position quartiles within each GB band; include all structural classes without PE filtering.
 q=np.zeros(len(pool),int)
 for g in (0,1):
  ix=np.flatnonzero(gb==g); zz=z[ix].copy(); zz[zz>.75]-=1
  order=ix[np.argsort(zz,kind='stable')]; q[order]=np.minimum(np.arange(len(ix))*4//max(len(ix),1),3)
 groups={}
 for j,key in enumerate(zip(gb,q,cna,ptm)): groups.setdefault(tuple(map(int,key)),[]).append(int(pool[j]))
 stacks=[list(rng.permutation(v)) for k,v in sorted(groups.items())]
 out=[]
 while len(out)<n:
  for st in stacks:
   if st and len(out)<n: out.append(st.pop())
 return np.array(out,int)

manifest={'seed':SEED,'created':'2026-11-01','purpose':'formal 240 configuration fixed-center Ni SOAP sweep','grid':{'r_cut':[3.25,4,5,6],'sigma':[.3,.5,.7],'n_max':[6,8,10,12],'l_max':[4,6,8,10,12]},'descriptor':{'species':['Ni'],'periodic':True,'rbf':'gto','average':'off','dtype':'float64','weighting':None},'reference_selection':'geometrically >9 Angstrom from every previously identified interface atom along periodic z; no PE filtering; excluded entire prior 4659 population','bulk_margin_A':9,'structures':[],'convergence_thresholds':{'pair_spearman_min':.98,'median_relative_change_max':.05,'top10_overlap_min':.8},'limitations':['two pre-existing Ni structures, not an independent GB ensemble','thresholds are engineering criteria, not physical proof of metastability','normalized and raw distances audited in parallel']}
centers=[];refs=[];population=[];center_info=[]
for state,(s,at) in enumerate(zip(labels,atoms)):
 ids=get_auxiliary(at,'id').astype(int); idmap={v:i for i,v in enumerate(ids)}
 sr=[r for r in rows if r['structure']==s]
 iface=np.array([idmap[int(r['atom_id'])] for r in sr if r['region']=='interface'])
 oldbulk=np.array([idmap[int(r['atom_id'])] for r in sr if r['region']=='bulk'])
 pop=np.array([idmap[int(r['atom_id'])] for r in sr])
 z=at.get_scaled_positions(wrap=True)[:,2]; L=at.cell.lengths()[2]
 iz=z[iface]
 # 47k x 2k is avoided with sorted interface positions and periodic search.
 zz=np.sort(np.concatenate([iz-1,iz,iz+1])); j=np.searchsorted(zz,z)
 dist=np.minimum(abs(z-zz[np.maximum(j-1,0)]),abs(z-zz[np.minimum(j,len(zz)-1)]))*L
 far=np.flatnonzero(dist>9)
 eligible=np.setdiff1d(far,pop)
 bulkpool=np.intersect1d(oldbulk,far)
 if len(bulkpool)<200: raise RuntimeError(f'{s} far original bulk only {len(bulkpool)}')
 ci=np.concatenate([stratified(iface,at,300),rng.choice(bulkpool,200,False)])
 ri=rng.choice(eligible,200,False)
 assert len(np.intersect1d(ci,ri))==0
 centers.append(ci);refs.append(ri);population.append(pop)
 for k,idx in enumerate(ci):
  center_info.append({'state':s,'region':'interface' if k<300 else 'bulk','atom_id':int(ids[idx]),'atom_index':int(idx),'z_A':float(at.positions[idx,2]),'z_scaled':float(z[idx]),'cna':int(get_auxiliary(at,'c_CNA')[idx]),'ptm':int(get_auxiliary(at,'c_PTM_1')[idx])})
 manifest['structures'].append({'state':s,'natoms':len(at),'sha256':hashlib.sha256((SOURCE/('ni_'+s+'.cfg')).read_bytes()).hexdigest(),'cell_A':at.cell.tolist(),'interface_population':len(iface),'old_bulk_population':len(oldbulk),'far_original_bulk':len(bulkpool),'reference_eligible':len(eligible),'center_ids':ids[ci].tolist(),'reference_ids':ids[ri].tolist(),'population_ids':ids[pop].tolist()})
(ROOT/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
(ROOT/'centers.json').write_text(json.dumps(center_info,indent=2),encoding='utf-8')
np.savez(ROOT/'selections.npz',**{f'{key}_{i}':arr for key,seq in [('centers',centers),('references',refs),('population',population)] for i,arr in enumerate(seq)})
# Fixed distance pairs, 600 per class and state. Same pairs for every descriptor.
pair_a=[];pair_b=[];pair_group=[]
for st in range(2):
 off=st*500; I=np.arange(off,off+300); B=np.arange(off+300,off+500)
 for kind,aa,bb in [('BB',B,B),('BI',B,I),('II',I,I)]:
  pairs=set()
  while len(pairs)<600:
   a=int(rng.choice(aa));b=int(rng.choice(bb))
   if a==b: continue
   if kind!='BI' and a>b:a,b=b,a
   pairs.add((a,b))
  for a,b in sorted(pairs): pair_a.append(a);pair_b.append(b);pair_group.append(f'{labels[st]}_{kind}')
pair_a=np.array(pair_a);pair_b=np.array(pair_b)
np.savez(ROOT/'fixed_pairs.npz',a=pair_a,b=pair_b,group=np.array(pair_group))

def compare(a,b):
 out={}
 for k in ('raw','norm'):
  x=a[k+'_pairs'];y=b[k+'_pairs']; scale=max(float(np.median(abs(y))),1e-12)
  rho=float(spearmanr(x,y).statistic)
  rel=float(np.median(abs(x-y)/np.maximum(abs(y),scale*1e-6)))
  overlap=float(np.mean([len(set(xi)&set(yi))/10 for xi,yi in zip(a[k+'_knn'],b[k+'_knn'])]))
  out[k]={'rho':rho,'median_relative_change':rel,'top10_overlap':overlap,'pass':bool(rho>=.98 and rel<=.05 and overlap>=.8)}
 return out
allrows=[];convergence=[]
for r,sigma in itertools.product(manifest['grid']['r_cut'],manifest['grid']['sigma']):
 cache={}
 for n,l in itertools.product(manifest['grid']['n_max'],manifest['grid']['l_max']):
  key=f'r{r:g}_s{sigma:g}_n{n}_l{l}'; jp=ROOT/(key+'.json');npz=ROOT/(key+'.npz')
  if jp.exists() and npz.exists():
   row=json.loads(jp.read_text());d=dict(np.load(npz));print('resume '+key,flush=True)
  else:
   t=time.perf_counter();row={'key':key,'r_cut':r,'sigma':sigma,'n_max':n,'l_max':l};d={}
   try:
    with warnings.catch_warnings(record=True) as ws:
     warnings.simplefilter('always')
     soap=SOAP(species=['Ni'],periodic=True,r_cut=r,sigma=sigma,n_max=n,l_max=l,rbf='gto',average='off',dtype='float64')
     blocks=[soap.create(at,centers=np.concatenate([ci,ri]),n_jobs=1) for at,ci,ri in zip(atoms,centers,refs)]
     row['warnings']=sorted(set(str(w.message) for w in ws))
    row['create_s']=time.perf_counter()-t
    x=np.vstack([b[:500] for b in blocks]);ref=np.vstack([b[500:] for b in blocks]);nx=np.linalg.norm(x,axis=1);nr=np.linalg.norm(ref,axis=1)
    row.update(dimension=x.shape[1],finite=bool(np.isfinite(x).all() and np.isfinite(ref).all()),zero_norms=int((nx<1e-14).sum()+(nr<1e-14).sum()),norm_min=float(nx.min()),norm_max=float(nx.max()))
    if not row['finite'] or row['zero_norms']:raise ValueError('nonfinite or zero SOAP norms')
    d['soap_norm']=nx
    for kind,xx,rr in [('raw',x,ref),('norm',x/nx[:,None],ref/nr[:,None])]:
     dd=cdist(xx,rr);near=np.argpartition(dd,19,axis=1)[:,:20];ds=np.take_along_axis(dd,near,axis=1);ds.sort(axis=1)
     for kk in (5,10,20): d[f'{kind}_d{kk}']=ds[:,:kk].mean(axis=1)
     d[kind+'_pairs']=np.linalg.norm(xx[pair_a]-xx[pair_b],axis=1)
     nn=cdist(xx,xx);np.fill_diagonal(nn,np.inf);d[kind+'_knn']=np.argpartition(nn,9,axis=1)[:,:10]
     for st in range(2):
      for region,sl in [('interface',slice(st*500,st*500+300)),('bulk',slice(st*500+300,st*500+500))]:
       vals=d[kind+'_d10'][sl];row[f'{kind}_{labels[st]}_{region}_d10_median']=float(np.median(vals))
    row['ok']=True
    np.savez_compressed(npz,**d)
   except Exception as e:
    row.update(ok=False,error=repr(e));np.savez_compressed(npz,failed=np.array([1]))
   row['total_s']=time.perf_counter()-t;jp.write_text(json.dumps(row,indent=2),encoding='utf-8')
   print(json.dumps({k:row.get(k) for k in ['key','ok','dimension','create_s','total_s','error']}),flush=True)
  allrows.append(row)
  if row.get('ok'):cache[n,l]=d
  (ROOT/'progress.json').write_text(json.dumps({'complete':len(allrows),'total':240,'last':key,'successful':sum(bool(a.get('ok')) for a in allrows)}),encoding='utf-8')
 # Separate n/l increments: retain paired metrics to demonstrate convergence, not descriptor size.
 for (n,l),d in cache.items():
  for axis,n2,l2 in [('n',n+2,l),('l',n,l+2)]:
   if (n2,l2) in cache:
    convergence.append({'r_cut':r,'sigma':sigma,'axis':axis,'n_from':n,'l_from':l,'n_to':n2,'l_to':l2,**compare(d,cache[n2,l2])})
 (ROOT/'convergence.json').write_text(json.dumps(convergence,indent=2),encoding='utf-8')
 (ROOT/'results.json').write_text(json.dumps(allrows,indent=2),encoding='utf-8')
keys=sorted(set().union(*(a.keys() for a in allrows)))
with (ROOT/'results.csv').open('w',newline='',encoding='utf-8-sig') as f:
 w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(allrows)
print('DONE '+str(len(allrows)),flush=True)
