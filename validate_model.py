"""Statistical and numerical checks; includes a selection-biased recovery test."""
import json
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from scipy.integrate import quad
from astropy.cosmology import Planck15
import model

ROOT=Path(__file__).resolve().parent
data=json.loads((ROOT/'results/model.json').read_text(encoding='utf-8'))
baseline=data['scenarios'][data['default']]
checks={}

# Check normalized spin components over the entire supported parameter space.
checks['spin_norm_max_error']=max(abs(quad(lambda x:model.aligned_tilt(x,s),-1,1)[0]-1) for s in model.SIGMAS)
checks['magnitude_norm_max_error']=max(abs(quad(lambda x:model.magnitude_pdf(x,a,b),0,1)[0]-1) for a,b in model.SPINS)
assert checks['spin_norm_max_error']<1e-8 and checks['magnitude_norm_max_error']<1e-8

# Recover known fraction from a deliberately biased, deterministic mock catalog.
# A(x)=2x, B(x)=2(1-x), p_detect(x)=0.2+0.6x; xi_A=.6, xi_B=.4.
truth=.3
x=np.linspace(0,1,100001)
pop=(truth*2*x+(1-truth)*2*(1-x))*(.2+.6*x)
cdf=np.cumsum(pop);cdf/=cdf[-1]
observed=np.interp((np.arange(1000)+.5)/1000,cdf,x)
f=(np.arange(2001)+.5)/2001
log_un=np.log(f[:,None]*2*observed+(1-f[:,None])*2*(1-observed)).sum(axis=1)
log_co=log_un-len(observed)*np.log(f*.6+(1-f)*.4)
checks['synthetic_truth']=truth
checks['synthetic_corrected_map']=float(f[np.argmax(log_co)])
checks['synthetic_uncorrected_map']=float(f[np.argmax(log_un)])
assert abs(checks['synthetic_corrected_map']-truth)<.003
assert abs(checks['synthetic_uncorrected_map']-truth)>.1

# Identical component populations must leave the mixture fraction unidentified.
same={**baseline,'event_moments':[[1,1,1,1,1]]*69,'selection':[.2,.2,.001,.001,.001]}
assert np.max(np.abs(model.density(same)-1))<1e-10
checks['identical_components_flat']=True

original_f=model.F.copy()
model.F=(np.arange(2001)+.5)/2001
fine=model.interval(model.density(baseline))
model.F=original_f
coarse=model.interval(model.density(baseline))
checks['grid_501_vs_2001_max_quantile_difference']=float(np.max(np.abs(np.array(fine)-coarse)))
assert checks['grid_501_vs_2001_max_quantile_difference']<.002

for record in data['scenarios'].values():
    for corrected in [True,False]:
        pdf=model.density(record,corrected)
        assert np.isfinite(pdf).all() and abs(pdf.mean()-1)<1e-10
checks['all_81_scenarios_finite_and_normalized']=True

# Bootstrap the actual Monte Carlo integration inputs, NOT astrophysical events.
# This isolates Monte Carlo variation from population sampling uncertainty.
model.ZGRID=np.linspace(0,2.3,3001)
model.ZDV=Planck15.differential_comoving_volume(model.ZGRID).value
model.ZNORM={lam:np.trapezoid(model.ZDV*(1+model.ZGRID)**(lam-1),model.ZGRID) for lam in model.LAMBDAS}
npz=np.load(ROOT/'data/posteriors.npz')
keys=['mass_1','mass_ratio','redshift','a_1','a_2','cos_tilt_1','cos_tilt_2','prior']
weights=[model.component_weights({k:npz[f'{s["index"]}_{k}'] for k in keys},1,0,2,.6) for s in data['events']]
inj=dict(np.load(ROOT/'data/injections.npz'))
ia,ib=model.component_weights(inj,1,0,2,.6)
total=float(inj['total_generated'])
rng=np.random.default_rng(20261022)
medians=[]
for replicate in range(40):
    em=[]
    for a,b in weights:
        idx=rng.integers(0,len(a),len(a));em.append([float(a[idx].mean()),float(b[idx].mean())])
    idx=rng.integers(0,len(ia),len(ia))
    rec={'event_moments':em,'selection':[float(ia[idx].sum()/total),float(ib[idx].sum()/total)]}
    medians.append(model.interval(model.density(rec))[1])
checks['mc_bootstrap_replicates']=40
checks['mc_bootstrap_median_sd']=float(np.std(medians,ddof=1))
checks['mc_bootstrap_median_range']=np.quantile(medians,[.05,.95]).tolist()
checks['baseline_conditional_interval']=coarse
checks['baseline_mc_diagnostics']=baseline['diagnostics']
(ROOT/'results/validation.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
print(json.dumps(checks,indent=2))
