"""Numerical checks for the declared retention hypothesis, not source proof."""
from pathlib import Path
import itertools
import json
import numpy as np
from scipy.special import gammaincc,gammaln
from scipy.stats import ncx2

BASE=Path(__file__).resolve().parent


def main():
    checks=[]
    for df in [4,6]:
        for snr in [5.,8.5,10.,20.]:
            n=np.arange(1000);mu=snr**2/2
            poisson=np.exp(n*np.log(mu)-mu-gammaln(n+1))
            assert abs(poisson.sum()-1)<1e-12
            for threshold in [8.5,10.]:
                mixture=np.sum(poisson*gammaincc(df/2+n,threshold**2/2))
                ref=ncx2.sf(threshold**2,df,snr**2)
                checks.append(abs(mixture/ref-1))
    assert max(checks)<1e-10
    # All eight outcomes of a three-row Bernoulli weighted-sum problem.
    p=np.array([.07,.38,.91]);w=np.array([[.2,.3],[1.1,.1],[.5,.7]])
    rows=np.array(list(itertools.product([0.,1.],repeat=3)))
    probabilities=np.prod(np.where(rows,p,1-p),axis=1)
    values=rows@w;mean=probabilities@values
    cov=np.einsum('n,ni,nj->ij',probabilities,values-mean,values-mean)
    exact=(w.T*(p*(1-p)))@w
    error=float(np.max(abs(cov-exact)));assert error<1e-14
    # Verify stored simulation covariance is recomputed from the actual draws.
    with np.load(BASE/'cumulative-retention-checkpoints.npz',allow_pickle=False) as cp:
        simulations=cp['simulated_selections'].reshape((-1,12));analytic=cp['analytic_covariance']
    ratio=np.diag(np.cov(simulations,rowvar=False,ddof=1))/np.diag(analytic)
    saved=json.loads((BASE/'cumulative-retention-results.json').read_text())['conditional_noise_simulation']['empirical_to_analytic_selection_variance']
    assert np.max(abs(ratio-saved))<1e-12
    result=dict(schema='blackhole-retention-numerical-validation/1',noncentral_survival_poisson_mixture_max_relative_error=float(max(checks)),
        exact_bernoulli_covariance_max_abs_error=error,empirical_variance_ratio_range=[float(ratio.min()),float(ratio.max())],
        scope='Distribution arithmetic and Bernoulli-sum covariance only. The empirically inferred 8.5 cutoff and historical noise-generator identity are NOT established by these checks. 256 empirical variances have finite simulation scatter.')
    (BASE/'cumulative-retention-validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
