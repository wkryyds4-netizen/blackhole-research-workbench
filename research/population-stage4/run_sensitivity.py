"""Conditional two-event replacement audit. Not calibrated waveform uncertainty."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import json
import sys
import hashlib
import numpy as np
from scipy.integrate import cumulative_trapezoid
from astropy.cosmology import Planck15
from recover_prior import redshift_grid, density_function, prior, LAL

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'research'))
import model
from population_functionals import aligned_cdf

KEYS=[f'{s}-{a}-0-1' for a in [0,1] for s in range(3)]
INDICES=[38,48]
PARAMS=['mass_1','mass_ratio','redshift','a_1','a_2','cos_tilt_1','cos_tilt_2']

def init_model():
    model.ZGRID=np.linspace(0,2.3,3001)
    model.ZDV=Planck15.differential_comoving_volume(model.ZGRID).value
    model.ZNORM={lam:np.trapezoid(model.ZDV*(1+model.ZGRID)**(lam-1),model.ZGRID) for lam in model.LAMBDAS}

def weights(d):
    values=[]
    for key in KEYS:
        s,a,m,z=map(int,key.split('-'))
        left,right=model.component_weights(d,a,m,model.LAMBDAS[z],model.SIGMAS[s])
        values.append(np.column_stack([left,right]))
    return np.stack(values,axis=1)

def coordinates(points=1001):
    r=np.linspace(float(aligned_cdf(0,1.2)),.5,points)
    fs=np.array([(.5-r)/(.5-float(aligned_cdf(0,model.SIGMAS[int(k[0])]))) for k in KEYS])
    return r,fs

def baseline_log(data,fs):
    logs=[]
    for ki,key in enumerate(KEYS):
        rec=data['scenarios'][key]
        em=np.array(rec['event_moments'])[:,:2]
        f=fs[ki]
        mix=f[None,:]*em[:,0,None]+(1-f[None,:])*em[:,1,None]
        xi=f*rec['selection'][0]+(1-f)*rec['selection'][1]
        logs.append(np.log(mix).sum(axis=0)-len(em)*np.log(xi))
    return np.array(logs)

def evaluate(data,old,new,points=1001):
    """new has shape [replicate, pilot, model, component]."""
    r,fs=coordinates(points)
    logs=np.broadcast_to(baseline_log(data,fs),(len(new),len(KEYS),len(r))).copy()
    for ei in range(2):
        before=fs*old[ei,:,0,None]+(1-fs)*old[ei,:,1,None]
        after=fs[None,:,:]*new[:,ei,:,0,None]+(1-fs[None,:,:])*new[:,ei,:,1,None]
        logs+=np.log(after)-np.log(before)[None,:,:]
    p=np.exp(logs-logs.max(axis=-1,keepdims=True))
    p/=np.trapezoid(p,r,axis=-1)[...,None]
    cdf=cumulative_trapezoid(p,r,axis=-1,initial=0)
    cdf/=cdf[:,:,-1,None]
    qs=np.array([[np.interp([.05,.5,.95],cdf[b,k],r) for k in range(len(KEYS))] for b in range(len(new))])
    return qs

def ess_info(w):
    r,fs=coordinates()
    mean=w.mean(axis=0);second=np.einsum('nki,nkj->kij',w,w)/len(w)
    e=fs*mean[:,0,None]+(1-fs)*mean[:,1,None]
    e2=fs**2*second[:,0,0,None]+(1-fs)**2*second[:,1,1,None]+2*fs*(1-fs)*second[:,0,1,None]
    ess=len(w)*e**2/e2
    minimum=ess.argmin(axis=1)
    maxfractions=[]
    for k,g in enumerate(minimum):
        one=fs[k,g]*w[:,k,0]+(1-fs[k,g])*w[:,k,1]
        maxfractions.append(float(one.max()/one.sum()))
    return dict(n=len(w),min_ess_over_common_support=ess.min(axis=1).tolist(),
                r_at_min_ess=r[minimum].tolist(),max_weight_fraction_at_min_ess=maxfractions)

def main():
    data=json.loads((ROOT/'results/model.json').read_text())
    init_model()
    local={};banks={};checks=[];counts={}
    with np.load(ROOT/'data/posteriors.npz',allow_pickle=False) as a, np.load(BASE/'pilot-bank.npz',allow_pickle=False) as bank:
        grid=redshift_grid(a);interp,norm=density_function(grid)
        for ix in INDICES:
            d={k:a[f'{ix}_{k}'] for k in PARAMS};d['prior']=a[f'{ix}_prior']
            assert np.allclose(prior(d,interp),d['prior'],rtol=1e-12,atol=0)
            local[ix]=weights(d)
            pos=next(i for i,e in enumerate(data['events']) if e['index']==ix)
            for ki,key in enumerate(KEYS):
                got=model.reduce_weights(local[ix][:,ki,0],local[ix][:,ki,1])
                ref=np.array(data['scenarios'][key]['event_moments'][pos])
                checks.append(float(np.max(abs(got/ref-1))))
            for label in ['phenom','mixed']:
                q={k:bank[f'{ix}_{label}_{k}'] for k in PARAMS}
                q['prior']=prior(q,interp)
                banks[ix,label]=weights(q)
                counts[f'{ix}_{label}']=ess_info(banks[ix,label])
            for name,cosmo in [('direct_planck',Planck15),('direct_lal',LAL)]:
                z=d['redshift']
                q=dict(d)
                # Use the same constant normalization: it cancels from conditional shape.
                q['prior']=cosmo.differential_comoving_volume(z).value*4*np.pi/(1+z)/norm*(1+z)**2*d['mass_1']/4
                banks[ix,name]=weights(q)
            counts[f'{ix}_original']=ess_info(local[ix])
    assert max(checks)<1e-12
    # Inspect the largest weights in the previously identified stress setting.
    _,fs=coordinates()
    k=0
    rmin=counts['38_phenom']['r_at_min_ess'][k]
    r,_=coordinates();g=int(np.argmin(abs(r-rmin)))
    one=fs[k,g]*banks[38,'phenom'][:,k,0]+(1-fs[k,g])*banks[38,'phenom'][:,k,1]
    with np.load(ROOT/'research/population-stage3/pilot-row-map.npz',allow_pickle=False) as rowmap:
        original_rows=set(rowmap['38_original_row_indices'].tolist())
    top=np.argsort(one)[-5:][::-1]
    tail=dict(event='GW190805_211137',key=KEYS[k],r=float(r[g]),
              full_release_n=len(one),original_n=len(local[38]),
              top_five=[dict(source_row=int(i),fraction_of_integral=float(one[i]/one.sum()),
                             included_in_original_subset=int(i) in original_rows) for i in top])
    old=np.stack([local[ix].mean(axis=0) for ix in INDICES])
    replacements={'original':old}
    for label in ['phenom','mixed']:
        for which in ['38','48','both']:
            new=old.copy()
            for ei,ix in enumerate(INDICES):
                if which=='both' or which==str(ix):new[ei]=banks[ix,label].mean(axis=0)
            replacements[f'{label}_{which}']=new
    for label in ['direct_planck','direct_lal']:
        replacements[f'{label}_both']=np.stack([banks[ix,label].mean(axis=0) for ix in INDICES])
    names=list(replacements);new=np.stack(list(replacements.values()))
    qs=evaluate(data,old,new)
    median_r=float(qs[0,0,1])
    median_f=(.5-median_r)/(.5-float(aligned_cdf(0,.3)))
    def at_median(w):
        w=median_f*w[:,0,0]+(1-median_f)*w[:,0,1]
        return w, float(w.sum()**2/np.sum(w*w))
    medweights,fulless=at_median(banks[38,'phenom'])
    _,oldess=at_median(local[38])
    tail['at_original_posterior_median']=dict(r=median_r,full_release_ess=fulless,original_subset_ess=oldess,
        largest_weight_fraction=float(medweights.max()/medweights.sum()),
        largest_source_row=int(medweights.argmax()),
        largest_in_original_subset=int(medweights.argmax()) in original_rows)
    refined=evaluate(data,old,new,2001)
    grid_error=float(np.max(abs(qs-refined)))
    assert grid_error<1e-6
    rng=np.random.default_rng(20261025)
    subs=np.empty((512,2,6,2))
    for b in range(512):
        for ei,ix in enumerate(INDICES):
            w=banks[ix,'phenom'];ids=rng.choice(len(w),len(local[ix]),replace=False)
            subs[b,ei]=w[ids].mean(axis=0)
    sub_q=evaluate(data,old,subs)
    full_q=qs[names.index('phenom_both')]
    diff=100*(sub_q[:,:,1]-full_q[:,1])
    rows={}
    for j,name in enumerate(names):
        rows[name]=dict(quantiles=qs[j].tolist(),median_shift_from_original_pp=(100*(qs[j,:,1]-qs[0,:,1])).tolist(),
            width_sigma03_minus12_pp=[float(100*(qs[j,0,1]-qs[j,2,1])),float(100*(qs[j,3,1]-qs[j,5,1]))])
    result=dict(schema='blackhole-population-two-event-sensitivity/1',keys=KEYS,indices=INDICES,
        prior_assumption='All released labels evaluated using recovered tutorial prior shape; Mixed/SEOB prior supports not fully established.',
        reference_protocol='Native stored tilt columns unchanged; pilot configurations at 10 and 20 Hz; no infinity substitution.',
        frozen_model_sha256=hashlib.sha256((ROOT/'results/model.json').read_bytes()).hexdigest(),
        moment_reconstruction_max_relative_error=max(checks),grid_1001_to_2001_max_quantile_error=grid_error,
        scenarios=rows,weight_diagnostics=counts,weight_tail_audit=tail,
        finite_release_thinning=dict(seed=20261025,repeats=512,replace=False,draws_per_event=6913,
            quantile_levels=[.05,.5,.95],median_minus_full_phenom_pp_quantiles=np.quantile(diff,[.05,.5,.95],axis=0).T.tolist(),
            sd_pp=diff.std(axis=0,ddof=1).tolist(),
            scope='Conditional subset variability of fixed finite IMRPhenomXPHM releases; not full waveform/PE uncertainty.'),
        unchanged='Other 67 events, all injections, selection integrals, model settings and common r prior.')
    (BASE/'sensitivity-results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    np.savez_compressed(BASE/'sensitivity-checkpoints.npz',old=old,new=new,names=np.array(names),quantiles=qs,
                        thinning_median_minus_full_pp=diff,subsample_means=subs)
    print(json.dumps(dict(shifts_pp={k:v['median_shift_from_original_pp'] for k,v in rows.items()},
        thinning=result['finite_release_thinning'],grid_error=grid_error),indent=2))

if __name__=='__main__':main()
