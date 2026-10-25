"""All 69 stored events: numerical concentration and separate influence diagnostics."""
from pathlib import Path
import sys
import json
import hashlib
import numpy as np
from scipy.integrate import cumulative_trapezoid

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
sys.path.insert(0,str(ROOT/'research/population-stage4'))
from run_sensitivity import init_model,weights,coordinates,baseline_log,KEYS,PARAMS,model


def quantiles(r,logs):
    p=np.exp(logs-np.max(logs,axis=-1,keepdims=True))
    c=cumulative_trapezoid(p,r,axis=-1,initial=0);c/=c[...,-1,None]
    return np.array([np.interp([.05,.5,.95],row,r) for row in c])


def compute(points=1001):
    init_model()
    data=json.loads((ROOT/'results/model.json').read_text())
    expected=json.loads((ROOT/'research/population-stage3/validation.json').read_text())['frozen_hashes']
    for name,sha in expected.items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha
    r,fs=coordinates(points);logs=baseline_log(data,fs);q=quantiles(r,logs)
    entries=[];max_moment_error=0.;cp={}
    with np.load(ROOT/'data/posteriors.npz',allow_pickle=False) as archive:
        for pos,event in enumerate(data['events']):
            ix=event['index'];d={k:archive[f'{ix}_{k}'] for k in PARAMS+['prior']}
            w=weights(d);n=len(w);rows=[]
            means=w.mean(axis=0)
            sumsq=np.einsum('nki,nkj->kij',w,w)
            sumw=w.sum(axis=0)
            numerator=(fs*sumw[:,0,None]+(1-fs)*sumw[:,1,None])**2
            denominator=fs**2*sumsq[:,0,0,None]+(1-fs)**2*sumsq[:,1,1,None]+2*fs*(1-fs)*sumsq[:,0,1,None]
            ess=numerator/denominator
            drop_logs=[];leave_logs=[];top_rows=[]
            for k,key in enumerate(KEYS):
                frozen=np.array(data['scenarios'][key]['event_moments'][pos])
                got=model.reduce_weights(w[:,k,0],w[:,k,1])
                err=float(np.max(abs(got/frozen-1)));assert err<1e-10
                max_moment_error=max(max_moment_error,err)
                f=float(np.interp(q[k,1],r,fs[k]))
                contribution=f*w[:,k,0]+(1-f)*w[:,k,1]
                fraction=contribution/contribution.sum();top=int(np.argmax(fraction));top_rows.append(top)
                before=fs[k]*means[k,0]+(1-fs[k])*means[k,1]
                after_means=(sumw[k]-w[top,k])/(n-1)
                after=fs[k]*after_means[0]+(1-fs[k])*after_means[1]
                drop_logs.append(logs[k]+np.log(after/before))
                selection=data['scenarios'][key]['selection']
                xi=fs[k]*selection[0]+(1-fs[k])*selection[1]
                # The leave-event-out diagnostic uses N=68, including xi^(-68).
                leave_logs.append(logs[k]-np.log(before)+np.log(xi))
                central=(r>=q[k,0])&(r<=q[k,2])
                rows.append(dict(key=key,r_median=float(q[k,1]),
                    ess_at_median=float(1/np.sum(fraction**2)),
                    min_ess_in_original_90pct_interval=float(ess[k,central].min()),
                    min_ess_over_common_support=float(ess[k].min()),
                    max_row_fraction_at_median=float(fraction[top]),largest_tutorial_row=top,
                    n_for_half_integral=int(np.searchsorted(np.cumsum(np.sort(fraction)[::-1]),.5)+1)))
            drop_q=quantiles(r,np.array(drop_logs));leave_q=quantiles(r,np.array(leave_logs))
            for k,row in enumerate(rows):
                row['drop_largest_row_median_shift_pp']=float(100*(drop_q[k,1]-q[k,1]))
                row['leave_event_out_median_shift_pp']=float(100*(leave_q[k,1]-q[k,1]))
            entries.append(dict(index=ix,event=event['name'],n_rows=n,
                identity_status='exact seven-coordinate match in stage three' if ix in [38,48] else 'proposed tutorial identity; summary match only',
                min_median_ess=min(x['ess_at_median'] for x in rows),
                max_row_fraction=max(x['max_row_fraction_at_median'] for x in rows),
                max_abs_drop_row_shift_pp=max(abs(x['drop_largest_row_median_shift_pp']) for x in rows),
                max_abs_leave_event_shift_pp=max(abs(x['leave_event_out_median_shift_pp']) for x in rows),settings=rows))
            cp[f'{ix}_ess']=ess;cp[f'{ix}_means']=means
            cp[f'{ix}_drop_quantiles']=drop_q;cp[f'{ix}_leave_quantiles']=leave_q
            cp[f'{ix}_top_rows']=np.array(top_rows)
    order={
        'lowest_median_ess':[x['index'] for x in sorted(entries,key=lambda x:(x['min_median_ess'],x['index']))],
        'largest_single_row_shift':[x['index'] for x in sorted(entries,key=lambda x:(-x['max_abs_drop_row_shift_pp'],x['index']))],
        'largest_event_shift':[x['index'] for x in sorted(entries,key=lambda x:(-x['max_abs_leave_event_shift_pp'],x['index']))]}
    candidates=sorted(set(i for ranked in order.values() for i in ranked[:3]))
    result=dict(schema='blackhole-catalog-numerical-screen/1',keys=KEYS,points=points,
        n_events=len(entries),n_settings=len(KEYS),baseline_quantiles=q.tolist(),events=entries,
        rankings=order,followup_indices=candidates,frozen_hashes=expected,
        max_raw_moment_relative_error=max_moment_error,
        selection_policy='Union of top 3 per ranking, original-index tie break; follow-up priority only, no exclusions.',
        interpretation='Conditional stored-array diagnostics; ESS ignores chain correlation and absent tails; neither shift is a calibrated uncertainty.')
    return result,dict(r=r,fs=fs,baseline_quantiles=q,**cp)


def main():
    result,cp=compute()
    (BASE/'catalog-screen.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    np.savez_compressed(BASE/'catalog-checkpoints.npz',**cp)
    names={e['index']:e for e in result['events']}
    print(json.dumps({key:[{k:names[ix][k] for k in ['index','event','min_median_ess','max_row_fraction','max_abs_drop_row_shift_pp','max_abs_leave_event_shift_pp']} for ix in order[:5]] for key,order in result['rankings'].items()},indent=2))


if __name__=='__main__':main()
