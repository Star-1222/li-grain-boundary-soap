"""Meaningful end-to-end consistency checks on scientific artifacts."""
import run_sweep as core
import json
import numpy as np
import pandas as pd

ROOT=core.ROOT
def main():
    frame=core.population()
    assert not frame.duplicated(['structure','atom_index']).any()
    assert frame.pilot.sum()==1000 and frame.perturb.sum()==120
    assert (frame.role=='evaluation').sum()==4659
    assert (frame.role=='reference').sum()==400
    assert (frame.role=='alternate_reference').sum()==400
    for role in ['reference','alternate_reference']:
        assert frame.loc[frame.role==role,'distance_to_interface_A'].min()>10
    rows=pd.read_csv(ROOT/'sweep_summary.csv');assert len(rows)==240
    assert (rows.status=='ok').all()
    max_abs=0.;max_rel=0.
    for c in core.configurations():
        cid=c['config_id'];vectors=np.load(ROOT/'vectors'/f'{cid}.npz')
        dimension=c['n_max']*(c['n_max']+1)//2*(c['l_max']+1)
        assert vectors['evaluation'].shape==(1000,dimension)
        assert vectors['reference'].shape==(400,dimension)
        assert np.isfinite(vectors['evaluation']).all() and np.isfinite(vectors['reference']).all()
        m=np.load(ROOT/'metrics'/f'{cid}.npz')
        for metric in ['raw','norm']:
            assert np.isfinite(m['pair_'+metric]).all()
            assert not (m['nn_'+metric]==np.arange(1000)[:,None]).any()
            assert np.all(m['D1_'+metric]<=m['D5_'+metric]+1e-13)
            assert np.all(m['D5_'+metric]<=m['D10_'+metric]+1e-13)
            assert np.all(m['D10_'+metric]<=m['D20_'+metric]+1e-13)
        residual=float(m['identity_max_absolute_residual'])
        relative=residual/max(float((m['pair_raw']**2).max()),1e-30)
        max_abs=max(max_abs,residual);max_rel=max(max_rel,relative)
        assert relative<1e-10
    candidates=json.loads((ROOT/'candidates.json').read_text(encoding='utf-8'))
    assert len(candidates)==6
    maximum_recompute_difference=0.
    for c in candidates:
        cid=c['config_id'];review=pd.read_csv(ROOT/'review'/f'{cid}.csv')
        assert len(review)==4659 and (~review.pilot).sum()==3659
        full=np.load(ROOT/'review'/f'{cid}_vectors.npz')['evaluation']
        pilot=np.load(ROOT/'vectors'/f'{cid}.npz')['evaluation']
        maximum_recompute_difference=max(maximum_recompute_difference,float(np.max(np.abs(full[review.pilot.to_numpy()]-pilot))))
        np.testing.assert_allclose(full[review.pilot.to_numpy()],pilot,rtol=1e-12,atol=1e-12)
        for metric in ['raw','norm']:
            np.testing.assert_allclose(review.loc[review.pilot,'D10_'+metric],
                np.load(ROOT/'metrics'/f'{cid}.npz')['D10_'+metric],rtol=1e-10,atol=1e-12)
    perturb=pd.read_csv(ROOT/'perturbation_atomwise.csv')
    assert len(perturb)==6*16*120*2
    inv=json.loads((ROOT/'invariance.json').read_text());assert inv['passed']
    assert len(list((ROOT/'representatives').glob('*.xyz')))==48
    figures=list((ROOT/'figures').glob('*.png'));assert len(figures)>=10
    report=(ROOT/'REPORT.md').read_text(encoding='utf-8')
    import re
    for path in re.findall(r'!\[[^\]]*\]\(([^)]+)\)',report):assert (ROOT/path).exists()
    checks=dict(verified_at=core.stamp(),grid_success=240,population=4659,heldout=3659,
        paired_raw_norm=True,max_distance_identity_absolute_residual=max_abs,
        max_distance_identity_relative_residual=max_rel,
        max_full_vs_pilot_vector_difference=maximum_recompute_difference,
        invariance=inv,perturbation_rows=len(perturb),representative_files=48,
        figure_count=len(figures),all_checks_passed=True)
    core.save_json(ROOT/'verification.json',checks)
    core.log('verification_complete',**checks)

if __name__=='__main__':main()
