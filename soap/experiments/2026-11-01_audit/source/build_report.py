"""Generate inspected static scientific figures and a Chinese evidence report."""
import run_sweep as core
import json
import math
import itertools
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.spatial import cKDTree

ROOT=core.ROOT
FIG=ROOT/'figures'
COLORS={'stable_bulk':'#0072B2','stable_interface':'#D55E00',
        'metastable_bulk':'#009E73','metastable_interface':'#CC79A7'}

def save(fig,name):
    fig.savefig(FIG/(name+'.png'),dpi=170,bbox_inches='tight')
    fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight')
    plt.close(fig)

def caption(c):
    return f"rc={c['r_cut']:g}, sigma={c['sigma']:g}, n/l={c['n_max']}/{c['l_max']}"

def md_table(frame,columns,labels=None,formats=None):
    labels=labels or columns;formats=formats or {}
    lines=['| '+' | '.join(labels)+' |','| '+' | '.join(['---']*len(columns))+' |']
    for row in frame.to_dict('records'):
        values=[]
        for key in columns:
            value=row[key]
            if isinstance(value,(bool,np.bool_)):text='是' if value else '否'
            elif key in formats:text=formats[key](value)
            elif isinstance(value,float):text=f'{value:.4g}'
            else:text=str(value)
            values.append(text.replace('|','/'))
        lines.append('| '+' | '.join(values)+' |')
    return '\n'.join(lines)

def grid_figures(convergence,candidates):
    cmap=ListedColormap(['#eeeeee','#56B4E9','#E69F00','#009E73'])
    norm=BoundaryNorm([-.5,.5,1.5,2.5,3.5],4)
    chosen={c['config_id'] for c in candidates}
    for metric in ['raw','norm']:
        fig,axes=plt.subplots(3,4,figsize=(14,10),sharex=True,sharey=True)
        for ri,sigma in enumerate(core.SIGMA):
            for ci,rc in enumerate(core.RC):
                ax=axes[ri,ci]
                part=convergence[(convergence.metric==metric)&(convergence.r_cut==rc)&(convergence.sigma==sigma)]
                values=np.full((4,5),np.nan)
                for row in part.to_dict('records'):
                    y=core.NMAX.index(row['n_max']);x=core.LMAX.index(row['l_max'])
                    values[y,x]=int(row['plateau_n'])+2*int(row['plateau_l'])
                    if row['warning_count']:
                        ax.plot(x,y,'x',color='black',ms=8,mew=1.5)
                    if row['config_id'] in chosen:ax.plot(x,y,'*',ms=14,mfc='white',mec='black')
                ax.imshow(values,origin='lower',cmap=cmap,norm=norm,aspect='auto',zorder=-1)
                ax.set_xticks(range(5),core.LMAX);ax.set_yticks(range(4),core.NMAX)
                ax.set_title(f'rc={rc:g} A, sigma={sigma:g} A',fontsize=10)
                if ri==2:ax.set_xlabel('l_max')
                if ci==0:ax.set_ylabel('n_max')
        fig.suptitle(f'{metric}: two-edge basis convergence; x = warning, star = review candidate',fontsize=12)
        fig.subplots_adjust(bottom=.13,top=.92,hspace=.28)
        cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=axes.ravel().tolist(),orientation='horizontal',fraction=.025,pad=.08)
        cb.set_ticks([0,1,2,3],labels=['neither axis','radial only','angular only','both axes'])
        save(fig,'convergence_'+metric)

def basis_figure(sweep):
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for rc in core.RC:
        sub=sweep[(sweep.r_cut==rc)&(sweep.sigma==.5)&(sweep.l_max==12)].sort_values('n_max')
        axes[0].semilogy(sub.n_max,sub.basis_condition_max,'o-',label=f'rc={rc:g} A')
        axes[1].semilogy(sub.n_max,sub.basis_residual_max,'o-',label=f'rc={rc:g} A')
    axes[0].set_ylabel('max cond(S) over angular channels')
    axes[1].set_ylabel('max ||B S B^T - I||_infinity')
    for ax in axes:ax.set_xlabel('n_max');ax.set_xticks(core.NMAX);ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('GTO basis diagnostics at l_max=12 (no sample-count dependence)')
    fig.tight_layout();save(fig,'basis_health')

def candidate_figures(candidates):
    fig,axes=plt.subplots(2,3,figsize=(14,8))
    for ax,c in zip(axes.flat,candidates):
        d=pd.read_csv(ROOT/'review'/f"{c['config_id']}.csv")
        for group in core.GROUPS:
            p=d[d.group==group]
            ax.scatter(p.D10_raw,p.D10_norm,s=4,alpha=.35,color=COLORS[group],label=group)
        ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel('D10 raw');ax.set_ylabel('D10 norm')
        ax.set_title(caption(c)+f"\nrho(all)={core.rho(d.D10_raw,d.D10_norm):.3f}",fontsize=9)
        ax.grid(alpha=.15)
    axes[0,0].legend(fontsize=6,markerscale=3)
    fig.suptitle('Same environments, two metrics; raw/norm amplitudes are not comparable')
    fig.tight_layout();save(fig,'raw_vs_norm')
    for metric in ['raw','norm']:
        fig,axes=plt.subplots(2,3,figsize=(14,8))
        for ax,c in zip(axes.flat,candidates):
            d=pd.read_csv(ROOT/'review'/f"{c['config_id']}.csv")
            for group in core.GROUPS:
                p=np.sort(d.loc[d.group==group,'D10_'+metric].to_numpy())
                ax.plot(p,np.arange(1,len(p)+1)/len(p),color=COLORS[group],label=group,lw=1.2)
            ax.set_xscale('log');ax.set_xlabel('D10 '+metric);ax.set_ylabel('empirical CDF')
            ax.set_title(caption(c),fontsize=9);ax.grid(alpha=.15)
        axes[0,0].legend(fontsize=7)
        fig.suptitle(f'Full evaluation population: {metric} (sampling labels, not independent replicas)')
        fig.tight_layout();save(fig,'distributions_'+metric)

def perturb_figures(candidates,perturb):
    fig,axes=plt.subplots(2,3,figsize=(14,8))
    for ax,c in zip(axes.flat,candidates):
        sub=perturb[(perturb.config_id==c['config_id'])&(perturb.kind=='random')&(perturb.region=='interface')]
        for metric,color in [('raw','#0072B2'),('norm','#D55E00')]:
            p=sub[sub.metric==metric]
            values=p.groupby('amplitude').distance_over_interface_pair_median
            med=values.median();lo=values.quantile(.05);hi=values.quantile(.95)
            ax.plot(med.index,med.values,'o-',color=color,label=metric)
            ax.fill_between(med.index,lo.values,hi.values,color=color,alpha=.12)
        ax.set_title(caption(c),fontsize=9);ax.set_xlabel('synthetic displacement RMS (A)')
        ax.set_ylabel('self-distance / interface-pair median');ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Interface response to synthetic perturbations; band = 5th-95th percentile')
    fig.tight_layout();save(fig,'perturbation_response')
    fig,axes=plt.subplots(2,3,figsize=(14,8))
    for ax,c in zip(axes.flat,candidates):
        sub=perturb[(perturb.config_id==c['config_id'])&(perturb.kind=='strain')&(perturb.region=='interface')]
        for metric,color in [('raw','#0072B2'),('norm','#D55E00')]:
            p=sub[sub.metric==metric].groupby('amplitude').distance_over_interface_pair_median.median()
            ax.plot(p.index*100,p.values,'o-',color=color,label=metric)
        ax.set_title(caption(c),fontsize=9);ax.set_xlabel('isotropic strain (%)')
        ax.set_ylabel('self-distance / interface-pair median');ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8);fig.suptitle('Strain response is reported, not minimized as a nuisance')
    fig.tight_layout();save(fig,'strain_response')

def spatial_figures(candidates):
    geometry=json.loads((ROOT/'geometry.json').read_text(encoding='utf-8'))
    fig,axes=plt.subplots(6,2,figsize=(13,17))
    for row,c in enumerate(candidates):
        d=pd.read_csv(ROOT/'review'/f"{c['config_id']}.csv")
        for col,metric in enumerate(['raw','norm']):
            ax=axes[row,col]
            for state,color in [('stable','#0072B2'),('metastable','#D55E00')]:
                part=d[d.structure==state]
                ax.scatter(part.z,part['D10_'+metric],s=3,alpha=.3,color=color,label=state)
            for lo,hi in geometry['stable']['runs_A']:ax.axvspan(lo,hi,color='grey',alpha=.15)
            ax.set_ylabel('D10 '+metric);ax.set_title(caption(c),fontsize=9)
            if row==5:ax.set_xlabel('z (A)')
            if row==0:ax.legend(fontsize=8)
    fig.suptitle('Spatial population check; grey = historical stable interface mask',fontsize=13)
    fig.tight_layout();save(fig,'spatial_profiles')

def representatives(candidates):
    """Export actual neighborhoods and plot their radial profiles and projections."""
    atoms=core.structures();out=ROOT/'representatives';out.mkdir(exist_ok=True)
    rows=[]
    for c in candidates:
        d=pd.read_csv(ROOT/'review'/f"{c['config_id']}.csv")
        chosen=[]
        for state in ['stable','metastable']:
            p=d[(d.structure==state)&(d.region=='interface')&(~d.pilot)].copy()
            p['rank_raw']=p.D10_raw.rank(pct=True);p['rank_norm']=p.D10_norm.rank(pct=True)
            p['rank_gap']=p.rank_raw-p.rank_norm
            for reason,idx in [('raw_high_norm_low',p.rank_gap.idxmax()),('norm_high_raw_low',p.rank_gap.idxmin()),
                               ('norm_low',p.D10_norm.idxmin()),('norm_high',p.D10_norm.idxmax())]:
                chosen.append((state,reason,d.loc[idx]))
        fig,axes=plt.subplots(2,4,figsize=(15,7))
        for ax,(state,reason,row) in zip(axes.flat,chosen):
            at=atoms[state];delta=at.positions-at.positions[int(row.atom_index)]
            lengths=at.cell.lengths();delta-=np.round(delta/lengths)*lengths
            radius=c['r_cut']+c['sigma']*np.sqrt(-2*np.log(.001))
            dist=np.linalg.norm(delta,axis=1);idx=np.flatnonzero(dist<=radius)
            order=np.argsort(dist[idx]);idx=idx[order]
            file=out/f"{c['config_id']}_{state}_{reason}_{int(row.atom_id)}.xyz"
            with file.open('w',encoding='utf-8') as f:
                f.write(f'{len(idx)}\ncenter_id={int(row.atom_id)}; radius={radius:.6f}; relative minimum-image coordinates\n')
                for i in idx:f.write('Ni '+' '.join(f'{v:.10f}' for v in delta[i])+'\n')
            ax.scatter(delta[idx,0],delta[idx,2],c=dist[idx],s=10,cmap='viridis',vmin=0,vmax=radius)
            ax.scatter([0],[0],marker='*',s=70,color='red');ax.set_aspect('equal')
            ax.set_xlabel('relative x (A)');ax.set_ylabel('relative z (A)')
            ax.set_title(f'{state}, {reason}\nid={int(row.atom_id)}, CNA/PTM={int(row.cna)}/{int(row.ptm)}',fontsize=8)
            rows.append(dict(config_id=c['config_id'],structure=state,reason=reason,atom_id=int(row.atom_id),
                D10_raw=float(row.D10_raw),D10_norm=float(row.D10_norm),norm=float(row.soap_norm),
                coordination=int(row.coordination_nominal_rc),radius_A=radius,file=str(file.relative_to(ROOT))))
        fig.suptitle(caption(c)+'; held-out examples selected to inspect metric disagreement',fontsize=11)
        fig.tight_layout();save(fig,'neighborhoods_'+c['config_id'])
    pd.DataFrame(rows).to_csv(ROOT/'representative_environments.csv',index=False)

def decisions(candidates,convergence,perturb,robust):
    rows=[]
    invariant=json.loads((ROOT/'invariance.json').read_text())
    for c in candidates:
        cid=c['config_id']
        for metric in ['raw','norm']:
            conv=convergence[(convergence.config_id==cid)&(convergence.metric==metric)].iloc[0]
            p=perturb[(perturb.config_id==cid)&(perturb.metric==metric)]
            inv=p[p.kind.isin(['rotation','translation','permutation'])]
            rand=p[(p.kind=='random')&(p.amplitude==.01)&(p.region=='interface')]
            retention=rand.nn10_retention.mean()
            r=robust[(robust.config_id==cid)&(robust.metric==metric)&(robust.subset=='heldout')]
            ref=r[r.variant=='alternate_reference']
            kv=r[r.variant.isin(['k5','k20'])]
            rows.append(dict(config_id=cid,metric=metric,r_cut=c['r_cut'],sigma=c['sigma'],
                n_max=c['n_max'],l_max=c['l_max'],dimension=c['dimension'],
                basis_converged=bool(conv.converged),warning_count=c['warning_count'],
                invariance_pass=bool((inv.relative_vector_change<1e-8).all()),
                noise_001_ratio_median=float(rand.distance_over_interface_pair_median.median()),
                noise_001_nn_retention=float(retention),
                alternate_ref_heldout_rho_min=float(ref.spearman.min()),
                k5_20_heldout_rho_min=float(kv.spearman.min()),
                selection_reason=c['selection_reason']))
    table=pd.DataFrame(rows);table.to_csv(ROOT/'candidate_decision_table.csv',index=False)
    decisions={}
    for metric in ['raw','norm']:
        # No post-hoc noise threshold: acceptable perturbation scale is not physically calibrated.
        eligible=table[(table.metric==metric)&table.basis_converged&table.invariance_pass&(table.warning_count==0)]
        if len(eligible):
            candidate=eligible.sort_values(['dimension','r_cut','sigma']).iloc[0]
            decisions[metric]=dict(status='provisional_basis_converged_candidate',config_id=candidate.config_id,
                note='Chosen for representation economy among audited candidates; does not establish a uniquely correct physical length scale or noise tolerance.')
        else:
            decisions[metric]=dict(status='no_audited_candidate_meets_predeclared_basis_convergence',config_id=None)
    core.save_json(ROOT/'decisions.json',dict(metrics=decisions,
        primary_distance='not_frozen',reason='raw and norm discard different information; current two configurations do not establish a universal choice.',
        invariance=invariant,limitations='No independent GB replication or physically calibrated perturbation tolerance. Held-out atoms are an internal check only.'))
    return table,decisions

def write_report(sweep,convergence,candidates,perturb,robust,table,decision):
    protocol=json.loads((ROOT/'protocol.json').read_text(encoding='utf-8'))
    inv=json.loads((ROOT/'invariance.json').read_text())
    edges=pd.read_csv(ROOT/'convergence_edges.csv')
    passed=convergence.groupby('metric').converged.sum().to_dict()
    warning=sweep[sweep.warning_count>0]
    summary=table[['config_id','metric','basis_converged','warning_count','noise_001_ratio_median','noise_001_nn_retention','alternate_ref_heldout_rho_min']]
    converged_grid=[]
    for (metric,rc,sigma),part in convergence.groupby(['metric','r_cut','sigma']):
        converged_grid.append(dict(metric=metric,r_cut=rc,sigma=sigma,converged_count=int(part.converged.sum())))
    pd.DataFrame(converged_grid).to_csv(ROOT/'converged_counts_by_scale.csv',index=False)
    elapsed=sweep.total_s.sum()
    lines=[
        '# Ni 双晶 SOAP 参数与距离联合遍历报告',
        '',f'生成时间：{core.stamp()}。所有“通过”均指本文预先登记的内部检查，不等于物理普适性。',
        '', '## 结果摘要',
        '',f'- 240 组描述符配置：成功 {int((sweep.status=="ok").sum())}，失败 {int((sweep.status!="ok").sum())}；含警告配置 {len(warning)}。',
        f'- 两种距离均完成评价；按两条连续相邻边检查，raw 有 {passed.get("raw",0)} 组、norm 有 {passed.get("norm",0)} 组达到径向和角向共同收敛参考线。',
        f'- 已对 6 组候选做扰动/不变性和 4659 个评价中心复核；其中 3659 个中心未参与初筛，另有共同参考 400 个、备用参考 400 个。',
        f'- 平移、整体旋转、同元素置换的不变性检查：{"通过" if inv["passed"] else "未通过"}；最大相对向量变化 {inv["max_relative_vector_change"]:.3g}。',
        '- 主距离未冻结。raw 保留幅度与方向，norm 去除向量整体幅度；本轮未把幅度直接解释成密度，也未按 stable/metastable 分离程度选参数。',
        '', '## 预先固定的方法',
        '', '- 数据为现有 stable、metastable 两个 Ni CFG，每个 47120 原子；输入哈希、环境版本和随机种子见 `protocol.json`。',
        '- 240 组 = r_cut [3.25,4,5,6] × sigma [0.3,0.5,0.7] × n_max [6,8,10,12] × l_max [4,6,8,10,12]。',
        '- 固定 DScribe GTO、Ni、PBC、逐中心、float64、无额外径向加权、无压缩；固定 DScribe 唯一项输出上的欧氏权重，未改成对称项重加权的原始 SOAP 核。',
        '- 初筛 1000 个中心：每个构型 300 界面和 200 深体相。抽样按界面分箱、CNA/PTM 和空间层次进行，未使用 PEATOM 或 SOAP 值。',
        '- 复现历史 200 分箱掩膜（非主导比例 >0.1，每箱 >=50），得到 1600/2059 个界面中心。深体相距离界面分箱边缘 >10 Å，且 CNA/PTM=1/1。',
        '- 完整结构始终参与邻居搜索；抽样只减少中心。最大补齐搜索半径约 8.60 Å，小于最短盒长的一半。',
        '- 主参考每个构型随机分层选 200 个深体相；备用每个构型另选 200 个。参考、备用参考、评价原子 ID 不重叠，不用能量过滤。',
        '- 每组保存原始向量，归一化向量由同一产物派生；raw 与 norm 分别计算 D(k)、近邻和距离排序。',
        '- 固定 20000 对分层原子对，包含四个人口组的十种组对类别。收敛逐类别检查最低 Spearman，不用总体相关掩盖组内变化。',
        '', '## 展开收敛',
        '', '固定 r_cut、sigma，分别增加 n_max 或 l_max。每条相邻边的条件：十种组对中最低距离 Spearman >=0.98、最大中位距离变化/高阶距离中位数 <=5%、1000 中心的平均前十近邻重合率 >=80%。',
        '', '每个轴需要两条连续边均通过：内部节点用前后两边，最高阶用最后两边；最低阶只有一边，因此本轮不能被认定为已完成两边收敛验证。边界不通过时明确保留“网格内未收敛”，未事后放宽阈值。近邻并列与小距离排序可使此判据偏严格，因此同时报告三个原始指标。',
        '', '![raw 收敛](figures/convergence_raw.png)',
        '', '![norm 收敛](figures/convergence_norm.png)',
        '', '网格颜色：灰=两个轴均未满足两边判据，蓝=仅径向，橙=仅角向，绿=两个轴；叉号=警告；星号=复核候选。详细边与类别汇总见 `convergence_edges.csv`、`convergence_summary.csv`。',
        '', '## 数值警告',
        '', '重叠矩阵来自 GTO 径向基函数，与原子样本数量无关。报告 cond(S) 及正交化残差 ||B S Bᵀ−I||∞，不把警告直接等同于所有距离失效。',
        '', '![基函数数值检查](figures/basis_health.png)',
        '', md_table(warning,['config_id','basis_condition_max','basis_residual_max'],['含警告配置','最大条件数','最大正交化残差']),
        '', '## 候选与扰动检查',
        '', '初筛后按收敛证据给每个 r_cut 保留一组，另保留一个展宽尺度和旧基准。所有候选在计算扰动和完整人口前即已写入 `candidates.json`；不是最终最优参数名单。',
        '', md_table(summary,list(summary.columns),['配置','距离','两轴收敛','警告数','0.01Å扰动/界面对距离','0.01Å近邻保留率','备用参考留出组最低ρ']),
        '', '120 个固定中心（四组各 30）；对完整结构施加三维 RMS 0.01/0.03/0.05 Å 的随机位移，各三个固定种子；同一随机模式随幅度缩放。另施加 ±0.5%、±1% 均匀应变，以及平移/旋转/置换。所有位移均为合成探针，不代表温度或真实轨迹。',
        '', '扰动图的分母是同构型界面—界面固定原子对距离的中位数，分子是扰动前后同中心的距离。它用于比较信号和扰动尺度，不能当作已校准的正确/错误分类阈值。阴影为中心和种子混合分布的 5%—95% 分位，不是独立样本置信区间。',
        '', '![随机扰动](figures/perturbation_response.png)',
        '', '![均匀应变](figures/strain_response.png)',
        '', '## 两种距离的比较',
        '', 'd_raw²=(||p_i||−||p_j||)²+||p_i||||p_j|| d_norm²。每个网格配置都验证这一恒等式。归一化改变距离排序与近邻关系，不能将二者简单视为单位换算。不同参数的距离绝对值也不能直接按大小判优。',
        '', '![raw 与 norm](figures/raw_vs_norm.png)',
        '', '![raw 分布](figures/distributions_raw.png)',
        '', '![norm 分布](figures/distributions_norm.png)',
        '', '范数与名义 r_cut 内配位数的相关、raw/norm 的组内相关，以及对共同体相均值的平方距离范数贡献均单独保存。相关不证明范数只反映密度，也不证明 norm 是纯键角度量。见 `physical_audit.csv` 和 `review/*_decomposition.csv`。',
        '', '## 完整人口与参考敏感性',
        '', '固定候选复算所有 4659 个中心；分别报告初筛/留出人口、两个周期界面，检查 k=1/5/10/20、备用共同参考和按构型参考。参考或 k 改变可能影响细小的组内排序；不能只用混合人口的高相关宣称稳健。',
        '', '![空间剖面](figures/spatial_profiles.png)',
        '', '详见 `review_summary.csv`、`reference_k_robustness.csv` 和 `candidate_agreement.csv`。近邻原子和两个周期界面均有相关性，3659 个留出原子不是独立晶界验证集。',
        '', '## 参数处置与适用范围',
        '', f'- raw：{decision["raw"]["status"]}；候选 {decision["raw"].get("config_id") or "无"}。',
        f'- norm：{decision["norm"]["status"]}；候选 {decision["norm"].get("config_id") or "无"}。',
        '- 若有满足基函数收敛、无警告和不变性检查的候选，只把维数较小者列为经济性优先的暂定候选，不宣称其物理观察范围最优。',
        '- 本轮没有物理校准的位移容差，也没有独立 GB 重复；不据此冻结唯一 r_cut、sigma 或主距离。没有通过的配置如实标为未收敛，未用最高阶代替真值。',
        '- 已导出每组候选的 8 个留出界面环境：raw/norm 排名分歧最大的两类，以及 norm 的高/低例子。XYZ 位于 `representatives/`，索引在 `representative_environments.csv`，投影见 `figures/neighborhoods_*.png`。这些是用于检查解释的例子，不是独立物理真值标签。',
        '', '## 计算与复现',
        '', f'- 240 组单组计时之和 {elapsed/60:.2f} 分钟（包含向量保存和初筛后处理，不含脚本准备）。',
        f'- Python {protocol["python"]}；DScribe {protocol["versions"]["dscribe"]}；NumPy {protocol["versions"]["numpy"]}；SciPy {protocol["versions"]["scipy"]}。单进程单线程。',
        '- `progress.jsonl` 保留阶段时间；`protocol.json` 保留事先固定的方法、输入哈希和环境；源代码快照在 `source/`。',
        '- 运行顺序：run_sweep.py prepare → sweep → convergence；validate_candidates.py perturb → review；build_report.py。从本工作区 `work/ni_soap_sweep_20261101/` 执行对应脚本。',
        '- `vectors/` 保留全部 240 组原始初筛向量；`metrics/` 保留两种距离、kNN 和 D(k)；`review/` 保留候选的完整人口向量及逐原子结果。',
        '', '方法参考：[DScribe SOAP 文档](https://singroup.github.io/dscribe/latest/tutorials/descriptors/soap.html)、[GTO 与截断扩展实现](https://singroup.github.io/dscribe/latest/_modules/dscribe/descriptors/soap.html)。输入来源为用户仓库，版本和实际文件哈希见协议。',
    ]
    # Explain the diagnostic failure mechanism rather than just reporting zero passes.
    radial=edges[edges.axis=='n'];angular=edges[edges.axis=='l']
    edge_counts=[]
    for (metric,axis),part in edges.groupby(['metric','axis']):
        edge_counts.append(dict(metric=metric,axis=axis,passed=int(part.passed.sum()),total=len(part),
            rho_pass=int((part.spearman_min>=.98).sum()),change_pass=int((part.relative_change_max<=.05).sum()),
            neighbor_pass=int((part.nn_overlap_mean>=.8).sum())))
    detail=['', '主要瓶颈是径向。以下为相邻边通过数量；三项必须同时满足才算边通过。', '',
        md_table(pd.DataFrame(edge_counts),['metric','axis','passed','total','rho_pass','change_pass','neighbor_pass'],
                 ['距离','轴','三项均通过','总边数','排序通过','幅度变化通过','近邻通过']),
        '', '六组候选 n_max=10→12 的每个原子对类别结果见 `radial_category_diagnostics.csv`。这能区分体相小距离排序与界面结构距离的变化，避免只凭最差汇总值解释机制。', '']
    lines[lines.index('## 数值警告'):lines.index('## 数值警告')]=detail
    extension=pd.read_csv(ROOT/'radial_extension/summary.csv')
    extension_edges=pd.read_csv(ROOT/'radial_extension/convergence_edges.csv')
    section=['', '## 径向边界补查：n_max=14、16', '',
        '因主网格没有径向边通过，在六组已选候选的同一 r_cut、sigma、l_max 上额外检查 14/16；补查协议在计算前单独登记，原 240 组结果和判据保持原样。', '',
        f'12 个补查配置中 {int((extension.status=="ok").sum())} 个生成向量、{int((extension.status=="failed").sum())} 个构造失败；可比较的两种距离相邻边中 {int(extension_edges.passed.sum())}/{len(extension_edges)} 通过原判据。', '',
        md_table(extension,['config_id','status','warning_count','basis_residual_max'],['补查配置','状态','警告数','正交化残差']), '',
        '例如 r_cut=6、sigma=0.3、l_max=8 的 n_max=16 虽能生成有限向量，其最大正交化残差约 0.361；不能把高阶有限结果当成准确真值。成功高阶配置和失败异常均保留在 `radial_extension/`。继续提高 n_max 未在这些补查配置中建立可信收敛。', '',
        '另用四个观察尺度、n_max=10→12 做了对称项权重诊断：按 DScribe 单元素 l-major 顺序，对 n≠n′ 分量乘 √2 后分别计算 raw/norm。其 8 条诊断边仍未通过原判据，因此本轮所检验的径向变化不能仅归结为是否恢复对称项权重。该诊断不是主距离的替换，也未扩展为另一轮完整参数选择。见 `symmetry_weight_diagnostic.csv`、`metric_diagnostics.json`。分量顺序与唯一项来自 [DScribe v2.1.2 源码](https://github.com/SINGROUP/dscribe/blob/v2.1.2/dscribe/ext/soapGTO.cpp#L2240)。', '']
    lines[lines.index('## 数值警告'):lines.index('## 数值警告')]=section
    audit=pd.read_csv(ROOT/'physical_audit.csv')
    same=audit[audit.group.str.endswith('_interface') & audit.rho_Draw_Dnorm.notna()]
    pivot=same.pivot(index='config_id',columns='group',values='rho_Draw_Dnorm').reset_index()
    distance_detail=['', '下表是完整界面人口内部的 D_raw/D_norm 排序相关；比混合全部体相和界面的总体相关更能显示方法差别。', '',
        md_table(pivot,['config_id','stable_interface','metastable_interface'],['配置','stable 界面ρ','metastable 界面ρ']), '']
    lines[lines.index('## 完整人口与参考敏感性'):lines.index('## 完整人口与参考敏感性')]=distance_detail
    (ROOT/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    (ROOT/'README.md').write_text('# Ni SOAP 扫描结果\n\n完整中文报告：[REPORT.md](REPORT.md)。\n\n关键表：`sweep_summary.csv`、`convergence_summary.csv`、`candidate_decision_table.csv`、`reference_k_robustness.csv`。\n\n本轮同时验证 d_raw 和 d_norm，未冻结唯一主距离。\n',encoding='utf-8')

def main():
    FIG.mkdir(exist_ok=True)
    candidates=json.loads((ROOT/'candidates.json').read_text(encoding='utf-8'))
    sweep=pd.read_csv(ROOT/'sweep_summary.csv');conv=pd.read_csv(ROOT/'convergence_summary.csv')
    perturb=pd.read_csv(ROOT/'perturbation_atomwise.csv');robust=pd.read_csv(ROOT/'reference_k_robustness.csv')
    grid_figures(conv,candidates);basis_figure(sweep);candidate_figures(candidates)
    perturb_figures(candidates,perturb);spatial_figures(candidates);representatives(candidates)
    table,decision=decisions(candidates,conv,perturb,robust)
    write_report(sweep,conv,candidates,perturb,robust,table,decision)
    source=ROOT/'source';source.mkdir(exist_ok=True)
    for p in core.SOURCE.glob('*.py'):__import__('shutil').copy2(p,source/p.name)
    manifest=[]
    for folder in [source,ROOT/'inputs']:
        for path in folder.glob('*'):
            if path.is_file():manifest.append(dict(path=str(path.relative_to(ROOT)),sha256=core.sha(path),bytes=path.stat().st_size))
    core.save_json(ROOT/'source_input_manifest.json',manifest)
    core.log('report_generated',report=str(ROOT/'REPORT.md'),figures=len(list(FIG.glob('*.png'))),decisions=decision)

if __name__=='__main__':main()
