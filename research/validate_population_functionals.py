"""Independent quadrature and convergence audit of population_functionals.py."""
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('study',ROOT/'research/population_functionals.py')
study=importlib.util.module_from_spec(spec)
spec.loader.exec_module(study)


def run():
    source=ROOT/'results/model.json'
    data=json.loads(source.read_text(encoding='utf-8'))
    a=study.run_study(data,4001)
    b=study.run_study(data,8001)
    maximum={'quantiles_4001_8001':0.,'retained_mass_4001_8001':0.,
             'mc_passing_mass_4001_8001':0.,'tv_4001_8001':0.,
             'common_vs_conditioned_native_quad_quantiles':0.,
             'native_vs_quad_quantiles':0.,'retained_mass_vs_quad':0.,
             'old501_cdf_convention_vs_continuous':0.,'old501_density_vs_independent':0.}
    per_scenario={}
    for key,rec in data['scenarios'].items():
        left,right=a['scenarios'][key],b['scenarios'][key]
        for region in ['native','common_r']:
            pairs=[(left[region][v],right[region][v]) for v in ['f','r','a']] if region=='native' else [(left[region],right[region])]
            for x,y in pairs:
                maximum['quantiles_4001_8001']=max(maximum['quantiles_4001_8001'],np.max(np.abs(np.array(x['quantiles'])-y['quantiles'])))
            for quantity in ['prior','posterior']:
                maximum['mc_passing_mass_4001_8001']=max(maximum['mc_passing_mass_4001_8001'],abs(left[region]['mc_passing_mass'][quantity]-right[region]['mc_passing_mass'][quantity]))
        maximum['retained_mass_4001_8001']=max(maximum['retained_mass_4001_8001'],abs(left['common_r']['retained_original_posterior_mass']-right['common_r']['retained_original_posterior_mass']))
        # Construct likelihood independently using scalar loops and saved moments.
        em=np.asarray(rec['event_moments']); sel=rec['selection']; n=len(em)
        def logl(f):
            return sum(np.log(f*row[0]+(1-f)*row[1]) for row in em)-n*np.log(f*sel[0]+(1-f)*sel[1])
        oldf=np.array(data['f_grid'])
        oldlogs=np.array([logl(f) for f in oldf])
        logscale=float(np.max(oldlogs))
        fun=lambda f: np.exp(logl(f)-logscale)
        z=quad(fun,0,1,epsabs=1e-11,epsrel=1e-11)[0]
        quantiles=[brentq(lambda f: quad(fun,0,f,epsabs=1e-11,epsrel=1e-11)[0]/z-q,0,1,xtol=1e-12) for q in [.05,.5,.95]]
        maximum['native_vs_quad_quantiles']=max(maximum['native_vs_quad_quantiles'],np.max(np.abs(np.array(quantiles)-right['native']['f']['quantiles'])))
        # Exact conditioning of the original uniform-f prior on the shared r support.
        sigma=[.3,.6,1.2][int(key.split('-')[0])]
        c0,c1=study.functional_coefficients(sigma)['r']
        limit=(a['common_r_support'][0]-c0)/c1
        zc=quad(fun,0,limit,epsabs=1e-11,epsrel=1e-11)[0]
        condq=[c0+c1*brentq(lambda f: quad(fun,0,f,epsabs=1e-11,epsrel=1e-11)[0]/zc-q,0,limit,xtol=1e-12) for q in [.95,.5,.05]]
        error=float(np.max(np.abs(np.array(condq)-right['common_r']['quantiles'])))
        maximum['common_vs_conditioned_native_quad_quantiles']=max(maximum['common_vs_conditioned_native_quad_quantiles'],error)
        maximum['retained_mass_vs_quad']=max(maximum['retained_mass_vs_quad'],abs(zc/z-right['common_r']['retained_original_posterior_mass']))
        # Legacy model.interval uses cumulative midpoint masses without the half-cell correction.
        oldp=np.exp(oldlogs-np.max(oldlogs));oldp/=oldp.sum()
        oldq=np.interp([.05,.5,.95],np.cumsum(oldp),oldf)
        maximum['old501_cdf_convention_vs_continuous']=max(maximum['old501_cdf_convention_vs_continuous'],np.max(abs(oldq-quantiles)))
        assert np.max(abs(oldq-rec['corrected']))<1e-12
        newlogs=study.log_likelihood(rec,oldf)
        newp=np.exp(newlogs-np.max(newlogs));newp/=newp.sum()
        maximum['old501_density_vs_independent']=max(maximum['old501_density_vs_independent'],np.max(abs(newp-oldp)))
        per_scenario[key]={'conditioned_native_vs_common_r_quantile_max_error':error,'native_continuous_quantiles':quantiles,'legacy_quantiles':oldq.tolist()}
    for x,y in zip(a['groups'],b['groups']):
        for k in x['common_r_tv']:
            maximum['tv_4001_8001']=max(maximum['tv_4001_8001'],abs(x['common_r_tv'][k]-y['common_r_tv'][k]))
    # Constant likelihood, including selection: each event and selection are f independent.
    mock={'scenarios':{},'default':data['default'],'events':[{'name':'synthetic-flat-likelihood'}]}
    for key in data['scenarios']:
        mock['scenarios'][key]={'event_moments':[[1,1,1,1,1]],'selection':[1,1,.0001,.0001,.0001], 'sample_counts':[1000],'total_generated':10000}
    flat=study.run_study(mock,4001)
    uniform_errors=[]
    for row in flat['scenarios'].values():
        for region in [row['native']['f'],row['native']['r'],row['native']['a'],row['common_r']]:
            lo,hi=region['prior_support']
            expected=lo+(hi-lo)*np.array([.05,.5,.95])
            uniform_errors.append(max(abs(np.array(region['quantiles'])-expected)))
            assert abs(region['kl_to_uniform_prior_nats'])<1e-12
    maximum['constant_likelihood_uniform_quantile_error']=max(uniform_errors)
    assert maximum['quantiles_4001_8001']<2e-6
    assert maximum['native_vs_quad_quantiles']<1e-6
    assert maximum['common_vs_conditioned_native_quad_quantiles']<1e-6
    assert maximum['retained_mass_vs_quad']<1e-6
    assert maximum['tv_4001_8001']<2e-6
    assert maximum['old501_density_vs_independent']<1e-12
    assert maximum['constant_likelihood_uniform_quantile_error']<1e-12
    assert len(a['scenarios'])==81 and len(a['groups'])==27
    result={'status':'pass','date':'2026-10-24','input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'implementation_sha256':hashlib.sha256((ROOT/'research/population_functionals.py').read_bytes()).hexdigest(),
            'validation_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scenarios':81,'groups':27,'max_errors':{k:float(v) for k,v in maximum.items()},
            'notes':['Legacy 501 midpoint CDF is not identical to continuous endpoint/trapezoid CDF; likelihood agrees.',
                     'MC pass-mass has a step-function boundary; its grid convergence is slower than smooth posterior quantiles.',
                     'Common-r prior equals conditioning native uniform-f prior on its corresponding f interval; it changes the scientific prior.',
                     'These are mathematical and quadrature tests, not event-selection or Monte Carlo bias validation.'],
            'per_scenario':per_scenario}
    dest=ROOT/'research/population-stage1/independent-functional-validation.json'
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps({'status':result['status'],'max_errors':result['max_errors']},indent=2))


if __name__=='__main__':
    run()
