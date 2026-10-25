"""Paired SRS-without-replacement benchmark conditional on the saved O3 bank."""
from pathlib import Path
import hashlib
import itertools
import json
import sys
import numpy as np
from scipy.integrate import cumulative_trapezoid

BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE))
from run_selection_study import physical,init_model,weights,coordinates


def srs_covariance(values,n):
    m=len(values)
    return (1-n/m)/n*np.cov(values,rowvar=False,ddof=1)


def main():
    # Exhaustively verify the finite-population formula on a nontrivial tiny bank.
    tiny=np.array([[0.,0.],[1.,4.],[3.,2.],[9.,7.],[2.,11.]])
    exact=np.array([tiny[list(ix)].mean(axis=0) for ix in itertools.combinations(range(5),3)])
    covariance_error=float(np.max(abs(np.cov(exact,rowvar=False,ddof=0)-srs_covariance(tiny,3))))
    assert covariance_error<1e-12
    source=json.loads((BASE/'o3-provenance.json').read_text());study=json.loads((BASE/'selection-results.json').read_text())
    assert hashlib.sha256((BASE/'o3-bank.npz').read_bytes()).hexdigest()==source['cache_sha256']
    with np.load(BASE/'o3-bank.npz',allow_pickle=False) as bank:
        raw={k[5:]:bank[k] for k in bank.files if k.startswith('full_')}
        actual=bank['combined_o3_to_full']
    init_model();w=weights(physical(raw))
    found=np.maximum.reduce([raw[k] for k in raw if k.startswith('ifar_')])>1
    w*=found[:,None,None];m=len(w);n=len(actual)
    assert m==284672 and n==128454
    c=source['attributes']['analysis_time_s']/source['combined_attributes']['analysis_time_s']/source['attributes']['total_generated']
    fixed=np.array(study['selection_components']['o1'])+np.array(study['selection_components']['o2'])
    total=fixed+c*w.sum(axis=0)
    actual_xi=fixed+c*m/n*w[actual].sum(axis=0)
    assert np.max(abs(actual_xi/np.array(study['selection_components']['exact_mixture'])-1))<1e-11
    with np.load(BASE/'selection-checkpoints.npz',allow_pickle=False) as cp:
        means=cp['means'][cp['names'].tolist().index('five_full_o3')]
    r,fs=coordinates();numerator=np.log(fs[None,:,:]*means[:,:,0,None]+(1-fs[None,:,:])*means[:,:,1,None]).sum(axis=0)
    def quantiles(selection):
        logs=numerator-len(means)*np.log(fs*selection[:,0,None]+(1-fs)*selection[:,1,None])
        p=np.exp(logs-logs.max(axis=1)[:,None]);cdf=cumulative_trapezoid(p,r,axis=1,initial=0);cdf/=cdf[:,-1,None]
        return np.array([np.interp([.05,.5,.95],v,r) for v in cdf])
    reference=quantiles(total);actual_q=quantiles(actual_xi)
    rng=np.random.default_rng(20261025);reps=512;xis=[];qs=[]
    for i in range(reps):
        take=rng.choice(m,size=n,replace=False)
        xi=fixed+c*m/n*w[take].sum(axis=0)
        xis.append(xi);qs.append(quantiles(xi))
        if (i+1)%64==0:print(f'Conditional thinning: {i+1}/{reps}',flush=True)
    xis=np.array(xis);qs=np.array(qs)
    expected_cov=(c*m)**2*srs_covariance(w.reshape(m,-1),n)
    simulated_cov=np.cov(xis.reshape(reps,-1),rowvar=False,ddof=1)
    diagonal_ratio=np.diag(simulated_cov)/np.diag(expected_cov)
    standardized_mean=(xis.mean(axis=0)-total).ravel()/np.sqrt(np.diag(expected_cov)/reps)
    gaps=100*np.stack([qs[:,0,1]-qs[:,2,1],qs[:,3,1]-qs[:,5,1]],axis=1)
    full_gaps=100*np.array([reference[0,1]-reference[2,1],reference[3,1]-reference[5,1]])
    actual_gaps=100*np.array([actual_q[0,1]-actual_q[2,1],actual_q[3,1]-actual_q[5,1]])
    result=dict(schema='blackhole-conditional-o3-thinning/1',seed=20261025,replicates=reps,full_accepted_rows=m,retained_rows=n,
        protocol='Uniform simple random sampling without replacement on the saved accepted bank, all failed-detection contributions zero, shared indices across six models and two components.',
        full_bank_gaps_pp=full_gaps.tolist(),actual_combined_subset_gaps_pp=actual_gaps.tolist(),
        simulated_gap_5_50_95_pp=np.quantile(gaps,[.05,.5,.95],axis=0).T.tolist(),
        actual_gap_empirical_percentile=(100*np.mean(gaps<=actual_gaps,axis=0)).tolist(),
        simulated_median_shift_5_50_95_pp=np.quantile(100*(qs[:,:,1]-reference[None,:,1]),[.05,.5,.95],axis=0).T.tolist(),
        actual_median_shift_from_full_pp=(100*(actual_q[:,1]-reference[:,1])).tolist(),
        exact_small_bank_covariance_error=covariance_error,
        selection_covariance_diagonal_simulated_over_analytic=diagonal_ratio.tolist(),
        selection_mean_standardized_by_mc_standard_error=standardized_mean.tolist(),
        max_quantile_regression_error=float(np.max(abs(reference-np.array(study['scenarios']['five_full_o3']['quantiles'])))),
        o3_bank_sha256=source['cache_sha256'],
        limits='Calibrates only the specified artificial thinning design conditional on the full saved bank. Actual historical thinning law is unverified. Empirical percentile is not a p-value; intervals are neither astrophysical credible intervals nor calibrated uncertainty of the historical campaign. O1/O2, events, noise, spin reference, and population models are fixed.')
    assert result['max_quantile_regression_error']<1e-11
    (BASE/'thinning-results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    np.savez_compressed(BASE/'thinning-checkpoints.npz',selections=xis,quantiles=qs,analytic_selection_covariance=expected_cov,simulated_selection_covariance=simulated_cov)
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
