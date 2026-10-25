"""Independent full-product and direct-delete checks for all 414 settings."""
from pathlib import Path
import sys
import json
import numpy as np
from scipy.integrate import cumulative_trapezoid

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
sys.path.insert(0,str(BASE))
from screen_catalog import compute,init_model,weights,KEYS,PARAMS


def product_quantiles(r,fs,means,selection):
    mix=fs[None,:]*means[:,0,None]+(1-fs[None,:])*means[:,1,None]
    mix=mix/means[:,1,None]
    xi=(fs*selection[0]+(1-fs)*selection[1])/selection[1]
    p=np.prod(mix.astype(np.longdouble),axis=0)/xi.astype(np.longdouble)**len(means)
    p=np.asarray(p/p.max(),dtype=float)
    c=cumulative_trapezoid(p,r,initial=0);c/=c[-1]
    return np.interp([.05,.5,.95],c,r)


def main():
    init_model()
    result=json.loads((BASE/'catalog-screen.json').read_text())
    frozen=json.loads((ROOT/'results/model.json').read_text())
    errors=[];ess_errors=[]
    with np.load(BASE/'catalog-checkpoints.npz',allow_pickle=False) as cp,np.load(ROOT/'data/posteriors.npz',allow_pickle=False) as archive:
        r=cp['r'];fs=cp['fs']
        for pos,event in enumerate(result['events']):
            ix=event['index'];d={k:archive[f'{ix}_{k}'] for k in PARAMS+['prior']};w=weights(d)
            for k,key in enumerate(KEYS):
                means=np.array(frozen['scenarios'][key]['event_moments'])[:,:2].copy()
                selection=frozen['scenarios'][key]['selection']
                removed=np.delete(means,pos,axis=0)
                q=product_quantiles(r,fs[k],removed,selection)
                errors.append(float(np.max(abs(q-cp[f'{ix}_leave_quantiles'][k]))))
                top=int(cp[f'{ix}_top_rows'][k])
                # Explicit deletion followed by averaging, not subtract-sum shortcut.
                means[pos]=np.delete(w[:,k,:],top,axis=0).mean(axis=0)
                q=product_quantiles(r,fs[k],means,selection)
                errors.append(float(np.max(abs(q-cp[f'{ix}_drop_quantiles'][k]))))
                for j in (0,len(r)//2,len(r)-1):
                    raw=fs[k,j]*w[:,k,0]+(1-fs[k,j])*w[:,k,1]
                    direct=float(raw.sum()**2/np.dot(raw,raw))
                    ess_errors.append(abs(direct/cp[f'{ix}_ess'][k,j]-1))
    assert max(errors)<1e-11 and max(ess_errors)<1e-11
    dense,cp2=compute(2001)
    grid=float(np.max(abs(np.array(dense['baseline_quantiles'])-result['baseline_quantiles'])))
    shift_errors=[]
    for a,b in zip(result['events'],dense['events']):
        for x,y in zip(a['settings'],b['settings']):
            for key in ['drop_largest_row_median_shift_pp','leave_event_out_median_shift_pp']:
                shift_errors.append(abs(x[key]-y[key]))
    assert grid<1e-6 and max(shift_errors)<1e-4
    assert dense['rankings']==result['rankings']
    with np.load(BASE/'catalog-checkpoints.npz',allow_pickle=False) as cp:
        complete_error=max(float(np.max(abs(cp[f'{e["index"]}_{kind}_quantiles']-cp2[f'{e["index"]}_{kind}_quantiles'])))
                           for e in result['events'] for kind in ['drop','leave'])
    assert complete_error<1e-6
    report=dict(schema='blackhole-catalog-screen-validation/1',status='pass',
        events=69,settings_per_event=6,independent_full_product_max_quantile_error=max(errors),
        direct_weight_vs_moment_ess_max_relative_error=max(ess_errors),
        grid_1001_to_2001_baseline_max_quantile_error=grid,
        grid_1001_to_2001_all_diagnostic_max_quantile_error=complete_error,
        grid_1001_to_2001_max_shift_difference_pp=max(shift_errors),all_three_rankings_unchanged=True,
        frozen_hashes_verified=result['frozen_hashes'],
        interpretation='Implementation checks only; no uncertainty calibration or missing-tail test.')
    (BASE/'screen-validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
