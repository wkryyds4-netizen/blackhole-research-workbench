"""Paired search-policy and artificial-thinning experiment, fixed before run.

Reuse the 512 SRSWOR draws, seed 20261025, from thinning_benchmark.py.
All six models and both search policies share exactly the same retained rows.
Calibrates a declared artificial design conditional on this finite bank only.
"""
from pathlib import Path
import json
import hashlib
import numpy as np
from scipy.integrate import cumulative_trapezoid
from run_selection_study import physical,init_model,weights,coordinates
from thinning_benchmark import srs_covariance

BASE=Path(__file__).resolve().parent


def main():
    source=json.loads((BASE/'o3-provenance.json').read_text())
    study=json.loads((BASE/'selection-results.json').read_text())
    policy=json.loads((BASE/'pipeline-results.json').read_text())
    assert hashlib.sha256((BASE/'o3-bank.npz').read_bytes()).hexdigest()==source['cache_sha256']
    with np.load(BASE/'o3-bank.npz',allow_pickle=False) as bank:
        raw={k[5:]:bank[k] for k in bank.files if k.startswith('full_')}
        actual=bank['combined_o3_to_full']
    with np.load(BASE/'thinning-checkpoints.npz',allow_pickle=False) as bank:
        previous_xi=bank['selections'];previous_q=bank['quantiles']
    init_model();w=weights(physical(raw))
    matched=np.maximum.reduce([raw['ifar_gstlal'],raw['ifar_mbta'],raw['ifar_pycbc_bbh'],raw['ifar_pycbc_hyperbank']])>1
    union=matched|(raw['ifar_cwb']>1)
    # Policy is an additional axis; row selection is shared, including zero contributions.
    bank_w=np.stack([w*union[:,None,None],w*matched[:,None,None]],axis=1)
    m=len(w);n=len(actual);assert (m,n)==(284672,128454)
    factor=source['attributes']['analysis_time_s']/source['combined_attributes']['analysis_time_s']/source['attributes']['total_generated']
    fixed=np.array(study['selection_components']['o1'])+np.array(study['selection_components']['o2'])
    full=fixed[None,:,:]+factor*bank_w.sum(axis=0)
    historical=fixed[None,:,:]+factor*m/n*bank_w[actual].sum(axis=0)
    with np.load(BASE/'pipeline-checkpoints.npz',allow_pickle=False) as bank:means=bank['means']
    r,fs=coordinates()
    lognumerator=np.log(fs[None,:,:]*means[:,:,0,None]+(1-fs[None,:,:])*means[:,:,1,None]).sum(axis=0)
    def quantiles(xi):
        logp=lognumerator[None,:,:]-len(means)*np.log(fs[None,:,:]*xi[:,:,0,None]+(1-fs[None,:,:])*xi[:,:,1,None])
        p=np.exp(logp-logp.max(axis=2)[:,:,None]);cdf=cumulative_trapezoid(p,r,axis=2,initial=0);cdf/=cdf[:,:,-1,None]
        return np.array([[np.interp([.05,.5,.95],c,r) for c in policy_c] for policy_c in cdf])
    full_q=quantiles(full);historical_q=quantiles(historical)
    assert np.max(abs(full_q[1]-np.array(policy['matched_filter_quantiles'])))<1e-11
    rng=np.random.default_rng(20261025);reps=512;xis=[];qs=[]
    for i in range(reps):
        take=rng.choice(m,n,replace=False)
        xi=fixed[None,:,:]+factor*m/n*bank_w[take].sum(axis=0)
        xis.append(xi);qs.append(quantiles(xi))
        if (i+1)%64==0:print(f'Paired policy/thinning: {i+1}/{reps}',flush=True)
    xis=np.array(xis);qs=np.array(qs)
    regression_xi=float(np.max(abs(xis[:,0]/previous_xi-1)))
    regression_q=float(np.max(abs(qs[:,0]-previous_q)))
    assert regression_xi<1e-11 and regression_q<1e-11
    analytic=(factor*m)**2*srs_covariance(bank_w.reshape(m,-1),n)
    empirical=np.cov(xis.reshape(reps,-1),rowvar=False,ddof=1)
    variance_ratio=np.diag(empirical)/np.diag(analytic)
    # Shapes: repeat, policy, amplitude family.
    gaps=100*np.stack([qs[:,:,0,1]-qs[:,:,2,1],qs[:,:,3,1]-qs[:,:,5,1]],axis=-1)
    full_gaps=100*np.stack([full_q[:,0,1]-full_q[:,2,1],full_q[:,3,1]-full_q[:,5,1]],axis=-1)
    historical_gaps=100*np.stack([historical_q[:,0,1]-historical_q[:,2,1],historical_q[:,3,1]-historical_q[:,5,1]],axis=-1)
    delta=gaps[:,1]-gaps[:,0];full_delta=full_gaps[1]-full_gaps[0]
    interaction=delta-full_delta
    diagnostics=[]
    for k,label in enumerate(['Beta(1.5,8.5)','Beta(2,5)']):
        diagnostics.append(dict(amplitude=label,full_policy_gap_change_pp=float(full_delta[k]),
            paired_policy_gap_change_5_50_95_pp=np.quantile(delta[:,k],[.05,.5,.95]).tolist(),
            interaction_5_50_95_pp=np.quantile(interaction[:,k],[.05,.5,.95]).tolist(),
            union_matched_gap_correlation=float(np.corrcoef(gaps[:,0,k],gaps[:,1,k])[0,1]),
            policy_change_std_pp=float(np.std(delta[:,k],ddof=1)),
            unpaired_variance_sum_std_pp=float(np.sqrt(np.var(gaps[:,0,k],ddof=1)+np.var(gaps[:,1,k],ddof=1))),
            matched_gap_nonpositive_fraction=float(np.mean(gaps[:,1,k]<=0))))
    result=dict(schema='blackhole-paired-policy-thinning/1',seed=20261025,replicates=reps,
        protocol='SRS without replacement of 128454 from 284672 accepted O3 rows. Shared indices across two policies, six models and both components; failed detections contribute zero. Fixed events and O1/O2.',
        policies=['union_including_cwb','matched_filter_only'],full_gaps_pp=full_gaps.tolist(),historical_subset_gaps_pp=historical_gaps.tolist(),
        gap_5_50_95_pp=np.moveaxis(np.quantile(gaps,[.05,.5,.95],axis=0),0,-1).tolist(),
        matched_historical_gap_empirical_percentile=(100*np.mean(gaps[:,1]<=historical_gaps[1],axis=0)).tolist(),
        diagnostics=diagnostics,selection_variance_empirical_over_analytic=variance_ratio.tolist(),
        original_union_selection_relative_regression_error=regression_xi,original_union_quantile_regression_error=regression_q,
        o3_bank_sha256=source['cache_sha256'],
        limits='These are conditional artificial-thinning distributions, not astrophysical credible intervals, p-values, or calibrated errors of the actual historical injection campaign. Event/search compatibility remains documentary. No physical spin-reference correction.')
    (BASE/'paired-policy-results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    np.savez_compressed(BASE/'paired-policy-checkpoints.npz',selections=xis,quantiles=qs,gaps=gaps,analytic_covariance=analytic,empirical_covariance=empirical)
    print(json.dumps({k:result[k] for k in ['full_gaps_pp','gap_5_50_95_pp','diagnostics','original_union_quantile_regression_error']},indent=2),flush=True)


if __name__=='__main__':main()
