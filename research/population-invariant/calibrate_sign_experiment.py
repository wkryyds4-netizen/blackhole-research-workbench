"""Exact calibration in a fully specified selected/noisy-sign experiment.

This tests the likelihood/selection implementation against a generative model.
It is NOT calibration of real GW strain analyses, PE chains, or injection MC.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import binom
from run_control import posterior

BASE=Path(__file__).resolve().parent


def main():
    n=69;rows=[]
    for en,ep,eta in [(.5,.5,0.),(.25,.75,0.),(.25,.75,.2),(.25,.75,.35)]:
        # A source is negative with probability g. Detection is Bernoulli
        # with sign-dependent efficiency. Its observed sign is flipped with
        # probability eta, independently conditional on the true sign.
        probability=lambda g:(g*en*(1-eta)+(1-g)*ep*eta)/(g*en+(1-g)*ep)
        intervals=[];predictive=[];covered=[];errors=[]
        for k in range(n+1):
            negative=np.tile([en*(1-eta),ep*eta],(k,6,1))
            positive=np.tile([en*eta,ep*(1-eta)],(n-k,6,1))
            means=np.concatenate([negative,positive],axis=0)
            _,q,_=posterior(means,np.tile([en,ep],(6,1)),8001)
            assert np.max(abs(q-q[0]))==0
            scale=float(binom.pmf(k,n,np.clip(k/n,eta,1-eta)))
            density=lambda g:float(binom.pmf(k,n,probability(g)))/scale
            norm=quad(density,0,1,epsabs=1e-11,epsrel=1e-11)[0]
            cdf=lambda g:quad(density,0,g,epsabs=1e-11,epsrel=1e-11)[0]/norm
            exact=[brentq(lambda g:cdf(g)-p,0,1,xtol=1e-12) for p in [.05,.5,.95]]
            errors.append(float(np.max(abs(q[0]-exact))))
            intervals.append(q[0].tolist());predictive.append(scale*norm)
            covered.append(cdf(q[0,2])-cdf(q[0,0]))
        prior_predictive=float(sum(predictive))
        coverage=float(np.dot(predictive,covered))
        error=max(errors)
        assert abs(prior_predictive-1)<1e-10 and error<2e-6 and abs(coverage-.9)<3e-5
        intervals=np.array(intervals);fixed=[]
        for g in [.05,.1,.2,.35,.5,.8,.95]:
            inside=(intervals[:,0]<=g)&(intervals[:,2]>=g)
            rate=float(np.sum(binom.pmf(np.arange(n+1),n,probability(g))*inside))
            fixed.append(dict(true_g=g,coverage=rate))
        row=dict(negative_detection_efficiency=en,positive_detection_efficiency=ep,sign_flip_probability=eta,
            n_selected=n,all_possible_negative_counts=n+1,prior_predictive_normalization=prior_predictive,
            prior_predictive_90pct_interval_coverage=coverage,maximum_quantile_error=error,fixed_parameter_coverages=fixed,
            intervals_by_observed_negative_count=intervals.tolist())
        rows.append(row);print(json.dumps({k:row[k] for k in ['negative_detection_efficiency','positive_detection_efficiency','sign_flip_probability','prior_predictive_90pct_interval_coverage','maximum_quantile_error']}),flush=True)
    report=dict(schema='blackhole-selected-noisy-sign-exact-calibration/1',experiments=rows,
        method='Enumerate all70 outcomes at N69. Compare the population posterior routine against independent binomial probabilities, adaptive quadrature and root-finding. Integrate coverage over the specified uniform prior and selected-catalog generative distribution.',
        inference_code_sha256=hashlib.sha256((BASE/'run_control.py').read_bytes()).hexdigest(),
        limits='Algorithm calibration in a discrete known likelihood with known selection and sign errors. No GW strain simulation, realistic PE posterior approximation, uncertainty in historical detection efficiency, chain correlations, or astrophysical formation-channel coverage is established. Fixed-g coverages are reported, not asserted to be uniformly90%.')
    (BASE/'sign-experiment-calibration.json').write_text(json.dumps(report,indent=2),encoding='utf-8')


if __name__=='__main__':main()
