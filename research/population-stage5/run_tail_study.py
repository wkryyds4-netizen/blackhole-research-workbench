"""Same-waveform, two-prior representation and influential-draw diagnostics."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import sys
import json
import hashlib
import numpy as np
from scipy.integrate import cumulative_trapezoid

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
sys.path.insert(0,str(ROOT/'research/population-stage4'))
from recover_prior import LAL
from run_sensitivity import init_model,weights,coordinates,baseline_log,KEYS,PARAMS

def density_pair(d):
    z=d['redshift'];m=d['mass_1']
    dl=LAL.luminosity_distance(z).value
    ddz=dl/(1+z)+(1+z)*LAL.hubble_distance.value/LAL.efunc(z)
    jac=m*(1+z)**2/4
    return (LAL.differential_comoving_volume(z).value*4*np.pi/(1+z)*jac,
            dl**2*ddz*jac,dl)

def summarize(data,means,points=1001):
    r,fs=coordinates(points)
    pos=next(i for i,e in enumerate(data['events']) if e['index']==38)
    logs=baseline_log(data,fs)
    for k,key in enumerate(KEYS):
        old=np.array(data['scenarios'][key]['event_moments'][pos])[:2]
        before=fs[k]*old[0]+(1-fs[k])*old[1]
        after=fs[k]*means[k,0]+(1-fs[k])*means[k,1]
        logs[k]+=np.log(after/before)
    p=np.exp(logs-logs.max(axis=1)[:,None]);p/=np.trapezoid(p,r,axis=1)[:,None]
    c=cumulative_trapezoid(p,r,axis=1,initial=0);c/=c[:,-1,None]
    return np.array([np.interp([.05,.5,.95],c[k],r) for k in range(6)]),p

def weighted_diagnostic(w,d,reference_r,tracked_row=None):
    from population_functionals import aligned_cdf
    rows=[]
    low=(d['a_1']<.2)&(d['a_2']<.2)
    for k,key in enumerate(KEYS):
        sigma=[.3,.6,1.2][int(key[0])]
        f=(.5-reference_r[k])/(.5-float(aligned_cdf(0,sigma)))
        a=f*w[:,k,0]+(1-f)*w[:,k,1]
        frac=a/a.sum()
        rows.append(dict(key=key,r=float(reference_r[k]),ess=float(a.sum()**2/np.sum(a*a)),
            max_weight_fraction=float(frac.max()),largest_row=int(frac.argmax()),
            n_for_half_integral=int(np.searchsorted(np.cumsum(np.sort(frac)[::-1]),.5)+1),
            low_spin_count=int(low.sum()),low_spin_integral_fraction=float(frac[low].sum()),
            tracked_row_fraction=None if tracked_row is None else float(frac[tracked_row])))
    return rows

def main():
    init_model();data=json.loads((ROOT/'results/model.json').read_text())
    arrays={}
    with np.load(BASE/'representations.npz',allow_pickle=False) as bank:
        for label in ['cosmo','nocosmo']:
            arrays[label]={k:bank[f'{label}_{k}'] for k in PARAMS+['luminosity_distance']}
        mapping=bank['cosmo_to_nocosmo_rows']
    with np.load(ROOT/'data/posteriors.npz',allow_pickle=False) as a:
        local={k:a[f'38_{k}'] for k in PARAMS+['prior']}
    old_w=weights(local)
    means={'original':old_w.mean(axis=0)};ws={};checks={}
    for label,d in arrays.items():
        pc,pn,dl=density_pair(d)
        checks[label+'_dl_max_relative_error']=float(np.max(abs(dl/d['luminosity_distance']-1)))
        d['prior']=pc if label=='cosmo' else pn
        ws[label]=weights(d);means[label+'_lal']=ws[label].mean(axis=0)
        if label=='nocosmo':
            ratio=pc/pn
            wpc=weights(dict(d,prior=pc))
            identity=np.mean(wpc*ratio[:,None,None],axis=0)
            checks['change_of_measure_max_relative_error']=float(np.max(abs(identity/means[label+'_lal']-1)))
            assert checks['change_of_measure_max_relative_error']<1e-12
            # Weighted cosmo representation: no rejection/resampling of the original draws.
            means['nocosmo_weighted_cosmo']=identity/ratio.mean()
    # Stored redshift may have been numerically inverted; don't silently change coordinates.
    assert max(checks[k] for k in checks if '_dl_' in k)<1e-3
    local_pc,_,_=density_pair(local)
    means['original_lal']=weights(dict(local,prior=local_pc)).mean(axis=0)
    tracked_no=int(mapping[2691])
    for label,index in [('cosmo',2691),('nocosmo',tracked_no)]:
        w=ws[label]
        means[label+'_drop_tracked']=(w.sum(axis=0)-w[index])/(len(w)-1)
    summaries={};curves={};grid_errors=[]
    baseline_q,_=summarize(data,means['original'])
    for label,value in means.items():
        q,p=summarize(data,value);q2,_=summarize(data,value,2001)
        grid_errors.append(float(np.max(abs(q-q2))))
        summaries[label]=dict(quantiles=q.tolist(),median_shift_from_original_pp=(100*(q[:,1]-baseline_q[:,1])).tolist(),
            width_sigma03_minus12_pp=[float(100*(q[0,1]-q[2,1])),float(100*(q[3,1]-q[5,1]))])
        curves[label]=p
    assert max(grid_errors)<1e-6
    identity_q=np.array(summaries['nocosmo_weighted_cosmo']['quantiles'])
    exact_q=np.array(summaries['nocosmo_lal']['quantiles'])
    checks['weighted_representation_max_quantile_error']=float(np.max(abs(identity_q-exact_q)))
    assert checks['weighted_representation_max_quantile_error']<1e-12
    ref=baseline_q[:,1]
    diagnostics={'original':weighted_diagnostic(old_w,local,ref),
        'cosmo_lal':weighted_diagnostic(ws['cosmo'],arrays['cosmo'],ref,2691),
        'nocosmo_lal':weighted_diagnostic(ws['nocosmo'],arrays['nocosmo'],ref,tracked_no)}
    r,fs=coordinates()
    shape={}
    for label,w in means.items():
        shape[label]=(fs*w[:,0,None]+(1-fs)*w[:,1,None])/w[:,1,None]
    result=dict(schema='blackhole-population-tail-representation/1',event='GW190805_211137',keys=KEYS,
        reference_r=ref.tolist(),n={k:len(v['mass_1']) for k,v in arrays.items()},
        cosmology='Planck15_lal; fixed existing source coordinates and fixed population model',
        densities='cosmo: dVc/dz/(1+z); nocosmo: dL^2 * ddL/dz; both times m1*(1+z)^2/4; constants omitted',
        source_match_sha256=hashlib.sha256((BASE/'representation-match.json').read_bytes()).hexdigest(),
        scenarios=summaries,weight_diagnostics=diagnostics,checks=checks,max_grid_quantile_error=max(grid_errors),
        tracked_cosmo_row=2691,tracked_nocosmo_row=tracked_no,
        drop_warning='Single-point influence diagnostic only; the main analysis retains all rows.',
        interpretation='Dependent finite representations of one PE analysis; no independent observations or calibrated convergence claim.',
        unchanged='Other 68 events, selection integrals, six population settings, native 10 Hz tilt definition.')
    (BASE/'tail-results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    np.savez_compressed(BASE/'tail-checkpoints.npz',r=r,fs=fs,names=np.array(list(means)),
        means=np.stack(list(means.values())),densities=np.stack(list(curves.values())),
        anchored_event_curves=np.stack(list(shape.values())),
        prior_ratio_nocosmo=ratio,
        **{f'sorted_fraction_{label}':np.sort((f*w[:,0,0]+(1-f)*w[:,0,1])/(f*w[:,0,0]+(1-f)*w[:,0,1]).sum())[::-1]
           for label,w in [('original',old_w),('cosmo',ws['cosmo']),('nocosmo',ws['nocosmo'])]
           for f in [float(np.interp(ref[0],r,fs[0]))]})
    print(json.dumps(dict(n=result['n'],checks=checks,
        shifts={k:v['median_shift_from_original_pp'] for k,v in summaries.items()},
        low_sigma03={k:v[0] for k,v in diagnostics.items()}),indent=2))

if __name__=='__main__':main()
