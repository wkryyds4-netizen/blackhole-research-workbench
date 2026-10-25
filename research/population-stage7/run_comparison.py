"""Fixed-model, fixed-selection comparison of original posterior representations."""
from pathlib import Path
import argparse
import sys
import hashlib
import json
import numpy as np
from scipy.integrate import cumulative_trapezoid

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
sys.path.insert(0,str(ROOT/'research'))
sys.path.insert(0,str(ROOT/'research/population-stage5'))
from population_event import WeightedEvent
from run_tail_study import density_pair,init_model,weights,coordinates,baseline_log,KEYS,PARAMS


def summary(data,replacements,points=1001):
    r,fs=coordinates(points);logs=baseline_log(data,fs)
    for index,means in replacements.items():
        pos=next(i for i,e in enumerate(data['events']) if e['index']==index)
        old=np.array([data['scenarios'][key]['event_moments'][pos][:2] for key in KEYS])
        before=fs*old[:,0,None]+(1-fs)*old[:,1,None]
        after=fs*means[:,0,None]+(1-fs)*means[:,1,None]
        logs+=np.log(after/before)
    p=np.exp(logs-logs.max(axis=1)[:,None]);p/=np.trapezoid(p,r,axis=1)[:,None]
    c=cumulative_trapezoid(p,r,axis=1,initial=0);c/=c[:,-1,None]
    qs=np.array([np.interp([.05,.5,.95],row,r) for row in c])
    return qs,p


def describe(w,d,ref,fs,r,tracked_rows,row_ids):
    mean=w.mean(axis=0);second=np.einsum('nki,nkj->kij',w,w)/len(w)
    mixed=fs*mean[:,0,None]+(1-fs)*mean[:,1,None]
    square=fs**2*second[:,0,0,None]+(1-fs)**2*second[:,1,1,None]+2*fs*(1-fs)*second[:,0,1,None]
    ess=len(w)*mixed**2/square;out=[]
    low=(d['a_1']<.2)&(d['a_2']<.2)
    for k,key in enumerate(KEYS):
        f=float(np.interp(ref[k,1],r,fs[k]));raw=f*w[:,k,0]+(1-f)*w[:,k,1]
        frac=raw/raw.sum();top=int(frac.argmax());tracked=int(tracked_rows[k])
        central=(r>=ref[k,0])&(r<=ref[k,2])
        out.append(dict(key=key,n_rows=len(w),r=float(ref[k,1]),ess=float(1/np.sum(frac**2)),
            minimum_ess_in_original_90pct=float(ess[k,central].min()),
            maximum_contribution=float(frac[top]),largest_source_row=int(row_ids[top]),
            tracked_tutorial_source_row=int(row_ids[tracked]),tracked_tutorial_contribution=float(frac[tracked]),
            n_for_half_integral=int(np.searchsorted(np.cumsum(np.sort(frac)[::-1]),.5)+1),
            low_spin_rows=int(low.sum()),low_spin_integral_fraction=float(frac[low].sum())))
    return out,ess


def metadata(record,representation):
    return dict(event=record['event'],analysis_label=record['analysis_label'],source=record['sources'][representation],
        row_basis='zero-based row in specified official HDF analysis',
        measure=dict(coordinates=PARAMS,units={k:'solar_mass' if k=='mass_1' else 'dimensionless' for k in PARAMS}),
        reference=dict(spin_frequency_hz=record['metadata']['meta_data/meta_data/f_ref'][0],frame='native PE reference'),
        cosmology='Planck15_lal',prior=dict(representation=representation,
            formula='dVc/dz/(1+z)*m1*(1+z)^2/4' if representation=='cosmo' else 'DL^2*dDL/dz*m1*(1+z)^2/4',
            normalization='hyperparameter-independent constants omitted'),
        support_status=dict(declared_analytic_bounds=record['prior_audit']['bounds'],
            caveat='retained-row checks only; unobserved tails and full historical constraints not established'),
        identity_status='exact seven-coordinate match of tutorial and cosmo rows to original-prior release')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--indices',type=int,nargs='+',default=[19,20,25,50]);args=parser.parse_args()
    indices=args.indices;assert set(indices)<=set([19,20,25,50]) and len(set(indices))==len(indices)
    init_model();data=json.loads((ROOT/'results/model.json').read_text())
    screen=json.loads((ROOT/'research/population-stage6/catalog-screen.json').read_text())
    for file,sha in screen['frozen_hashes'].items():assert hashlib.sha256((ROOT/file).read_bytes()).hexdigest()==sha
    provenance=json.loads((BASE/'event-provenance.json').read_text())
    reference,_=summary(data,{})
    r,fs=coordinates();banks={name:{} for name in ['original','original_lal','cosmo','nocosmo','weighted_cosmo']}
    diagnostics=[];checks=[];cp={}
    with np.load(ROOT/'data/posteriors.npz',allow_pickle=False) as archive:
        for index in indices:
            record=next(x for x in provenance if x['index']==index)
            path=BASE/f'event-{index}.npz';assert hashlib.sha256(path.read_bytes()).hexdigest()==record['cache_sha256']
            original={k:archive[f'{index}_{k}'] for k in PARAMS+['prior']}
            top=np.array([x['largest_tutorial_row'] for x in next(e for e in screen['events'] if e['index']==index)['settings']])
            diagnostics_event=dict(index=index,event=record['event'],reference_frequency_hz=record['metadata']['meta_data/meta_data/f_ref'][0],representations={})
            with np.load(path,allow_pickle=False) as a:
                for label in ['original','original_lal','cosmo','nocosmo']:
                    d=original if label.startswith('original') else {k:a[f'{label}_{k}'] for k in PARAMS}
                    pc,pn,dl=density_pair(d)
                    prior=d['prior'] if label=='original' else pn if label=='nocosmo' else pc
                    w=weights(dict(d,prior=prior))
                    if label.startswith('original'):
                        means=w.mean(axis=0);tracked=top;row_ids=np.arange(len(prior))
                    else:
                        error=float(np.max(abs(dl/a[f'{label}_luminosity_distance']-1)));assert error<1e-3
                        event=WeightedEvent(d,np.arange(len(prior)),np.log(prior),np.zeros(len(prior)),metadata(record,label))
                        population=weights(dict(d,prior=np.ones(len(prior))))
                        with np.errstate(divide='ignore'):lp=np.log(population)
                        means=np.exp(event.log_integral(lp))
                        integral_error=float(np.max(abs(means/w.mean(axis=0)-1)));assert integral_error<1e-10
                        tracked=a[f'tutorial_to_{label}'][top];row_ids=event.row_ids
                        checks.append(dict(index=index,label=label,distance_max_relative_error=error,interface_integral_max_relative_error=integral_error))
                        if label=='nocosmo':
                            changed=event.reexpress_prior(np.log(pc),dict(formula='cosmo density on original rows'),event.metadata['support_status'])
                            lm=changed.log_integral(lp);delta=lm-event.log_integral(lp)
                            assert np.ptp(delta)<1e-10
                            banks['weighted_cosmo'][index]=np.exp(lm)
                            checks[-1]['reexpression_log_constant_range']=float(np.ptp(delta))
                    banks[label][index]=means
                    diag,ess=describe(w,d,reference,fs,r,tracked,row_ids)
                    diagnostics_event['representations'][label]=diag
                    cp[f'{index}_{label}_means']=means;cp[f'{index}_{label}_ess']=ess
            diagnostics.append(diagnostics_event)
            print(json.dumps(dict(event=record['event'],low_spin_sigma03={k:{field:v[0][field] for field in ['n_rows','ess','maximum_contribution','tracked_tutorial_contribution','low_spin_rows']} for k,v in diagnostics_event['representations'].items()}),indent=2),flush=True)
    scenarios={'original':{}}
    for index in indices:
        for label in ['original_lal','cosmo','nocosmo']:
            scenarios[f'single_{index}_{label}']={index:banks[label][index]}
    for label in ['original_lal','cosmo','nocosmo','weighted_cosmo']:
        scenarios[f'all_new_{label}']=banks[label]
    # Add the previously validated GW190805 as a separate five-event comparison.
    with np.load(ROOT/'research/population-stage5/tail-checkpoints.npz',allow_pickle=False) as p5:
        for label in ['cosmo','nocosmo']:
            means=p5['means'][p5['names'].tolist().index(label+'_lal')]
            scenarios[f'new_plus_38_{label}']=dict(banks[label],**{})
            scenarios[f'new_plus_38_{label}'][38]=means
    summaries={};curves=[];grid_error=0.
    for name,replacements in scenarios.items():
        q,p=summary(data,replacements);fine,_=summary(data,replacements,2001)
        grid_error=max(grid_error,float(np.max(abs(q-fine))));assert grid_error<1e-6
        summaries[name]=dict(replaced_indices=sorted(replacements),quantiles=q.tolist(),
            median_shift_pp=(100*(q[:,1]-reference[:,1])).tolist(),
            sigma03_minus12_gap_pp=[float(100*(q[0,1]-q[2,1])),float(100*(q[3,1]-q[5,1]))])
        curves.append(p)
    identity=float(np.max(abs(np.array(summaries['all_new_nocosmo']['quantiles'])-summaries['all_new_weighted_cosmo']['quantiles'])))
    assert identity<1e-11
    suffix='' if indices==[19,20,25,50] else 'partial-'+'-'.join(map(str,indices))+'-'
    result=dict(schema='blackhole-multiple-event-representation/1',indices=indices,keys=KEYS,
        n_events_fixed=69,reference_quantiles=reference.tolist(),diagnostics=diagnostics,scenarios=summaries,
        checks=checks,max_grid_quantile_error=grid_error,weighted_representation_max_quantile_difference=identity,
        frozen_hashes=screen['frozen_hashes'],
        limits='Targeted event subset; dependent released samples; native spin references retained. No calibrated PE uncertainty or formation-channel measurement.')
    (BASE/(suffix+'comparison-results.json')).write_text(json.dumps(result,indent=2),encoding='utf-8')
    all_means=np.array([np.array([replacements.get(e['index'],np.array([data['scenarios'][k]['event_moments'][pos][:2] for k in KEYS])) for pos,e in enumerate(data['events'])]) for replacements in scenarios.values()])
    np.savez_compressed(BASE/(suffix+'comparison-checkpoints.npz'),r=r,fs=fs,names=np.array(list(scenarios)),
                        all_event_means=all_means,densities=np.array(curves),**cp)
    print(json.dumps(dict(scopes={k:{field:v[field] for field in ['replaced_indices','median_shift_pp','sigma03_minus12_gap_pp']} for k,v in summaries.items() if k.startswith(('all_','new_'))},grid_error=grid_error),indent=2))


if __name__=='__main__':main()
