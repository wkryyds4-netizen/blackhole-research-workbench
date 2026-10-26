"""Full69,54-model conditional controls on both official O3 banks.

Reads existing verified sources; writes only this study's versioned artifacts.
Run directly, never by executing root model.py. No ESS-based sample rejection.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import hashlib
import json
import sys
from importlib.metadata import version
import h5py
import numpy as np
from scipy.integrate import quad
from scipy.special import betainc
from grid_model import MODELS,weights,mass_pdf

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
CAT=ROOT/'research/population-catalog'
INV=ROOT/'research/population-invariant'
SEL=ROOT/'research/population-selection-audit'
MIX=ROOT/'research/population-mixture'
sys.path.insert(0,str(SEL))
from run_selection_study import init_model,density_pair,physical,PARAMS
from run_control import posterior
import model


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def relative(actual,expected):
    assert np.array_equal(actual==0,expected==0),'Zero components differ'
    take=expected!=0
    return float(np.max(abs(actual[take]/expected[take]-1))) if take.any() else 0.


def write_json(path,value):
    temporary=path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')
    temporary.replace(path)


def moments(w):
    return w.sum(axis=0),np.einsum('nki,nkj->kij',w,w),w.max(axis=0)


def new_stats():
    return dict(sum=np.zeros((54,2)),outer=np.zeros((54,2,2)),maximum=np.zeros((54,2)),rows=0)


def accumulate(stats,w):
    if not len(w):return
    total,outer,maximum=moments(w)
    stats['sum']+=total;stats['outer']+=outer
    stats['maximum']=np.maximum(stats['maximum'],maximum);stats['rows']+=len(w)


def events():
    frozen=json.loads((ROOT/'results/model.json').read_text())
    sources=json.loads((CAT/'event-provenance.json').read_text())
    old=json.loads((INV/'event-integrals.json').read_text())
    pieces=[BASE/'grid_model.py',BASE/'design.json',ROOT/'model.py',INV/'run_control.py',INV/'invariant_model.py',ROOT/'research/population-stage5/run_tail_study.py',ROOT/'research/population-stage4/recover_prior.py']
    protocol=hashlib.sha256(b'mass-redshift-event-integrals/1'+b''.join(p.read_bytes() for p in pieces)).hexdigest()
    cache_path=BASE/'event-integrals.json'
    cache=json.loads(cache_path.read_text()) if cache_path.exists() else dict(protocol=protocol,events=[])
    assert cache['protocol']==protocol,'Numerical protocol changed: version the cache before reuse.'
    used=[]
    for position,event in enumerate(frozen['events']):
        index=event['index'];info=next(r for r in sources if r['index']==index)
        path=CAT/f'event-{index}.npz'
        assert sha(path)==info['cache_sha256']
        saved=next((r for r in cache['events'] if r['index']==index and r['source_cache_sha256']==info['cache_sha256']),None)
        if saved is None:
            with np.load(path,allow_pickle=False) as bank:
                raw={k:bank[k] for k in PARAMS}
            rows=len(raw['mass_1']);stats=new_stats()
            for start in range(0,rows,32768):
                d={k:v[start:start+32768] for k,v in raw.items()}
                _,prior,_=density_pair(d);d['prior']=prior
                accumulate(stats,weights(d))
            means=stats['sum']/rows
            expected=np.array(next(r for r in old['events'] if r['index']==index)['mean'])
            bridge=relative(means[6:12],expected);assert bridge<1e-10
            saved=dict(index=index,event=info['event'],rows=rows,source_cache_sha256=info['cache_sha256'],
                       mean=means.tolist(),second=(stats['outer']/rows).tolist(),maximum=stats['maximum'].tolist(),
                       baseline_relative_error=bridge)
            cache['events']=[r for r in cache['events'] if r['index']!=index]+[saved]
            write_json(cache_path,cache)
            print(f'Full69 grid event {position+1}/69: {info["event"]}, {rows} rows',flush=True)
        used.append(saved)
    assert len(used)==69 and sum(r['rows'] for r in used)==11180578
    return used,protocol


def injection_stats(raw,found,masks,factors):
    # Keep every declared detected row. Targets outside their mass/z supports
    # have zero weights, rather than being silently converted into detections.
    indices=np.flatnonzero(found)
    output={name:new_stats() for name in masks}
    for start in range(0,len(indices),32768):
        take=indices[start:start+32768]
        d=physical({k:v[take] for k,v in raw.items()})
        w=weights(d)*factors[take,None,None]
        for name,mask in masks.items():accumulate(output[name],w[mask[take]])
    return output


def selections():
    source=json.loads((SEL/'o3-provenance.json').read_text())
    assert sha(SEL/'o3-bank.npz')==source['cache_sha256']
    with np.load(SEL/'o3-bank.npz',allow_pickle=False) as bank:
        full={k[5:]:bank[k] for k in bank.files if k.startswith('full_')}
        old={k[9:]:bank[k] for k in bank.files if k.startswith('combined_') and k!='combined_o3_to_full'}
    search=sorted(k for k in full if k.startswith('ifar_'))
    early=(old['name']!=b'o3')&(np.any([old[k]>1 for k in search],axis=0)|(old['optimal_snr_net']>10))
    early_stats=injection_stats(old,early,{'early':early},old['mixture_weight']/source['combined_attributes']['total_generated'])['early']
    time_fraction=source['attributes']['analysis_time_s']/source['combined_attributes']['analysis_time_s']
    out={};bank_info={}
    for bank_name in ['pure_bbh','official_mixture']:
        if bank_name=='pure_bbh':
            raw=full;attrs=source['attributes']
        else:
            record=json.loads((MIX/'download.json').read_text());path=MIX/'sources'/record['file']
            assert sha(path)==record['sha256']
            keys=['mass1_source','mass2_source','redshift','sampling_pdf','mixture_weight']
            keys+=[f'spin{i}{axis}' for i in [1,2] for axis in 'xyz']+search
            with h5py.File(path,'r') as handle:
                raw={k:handle['injections'][k][:] for k in keys}
                attrs={k:float(handle.attrs[k]) for k in ['analysis_time_s','total_generated']}
            assert attrs['analysis_time_s']==source['attributes']['analysis_time_s']
        union=np.any([raw[k]>1 for k in search],axis=0)
        matched=np.any([raw[k]>1 for k in search if k!='ifar_cwb'],axis=0)
        stats=injection_stats(raw,union,dict(matched_filter=matched,with_cwb=union),raw['mixture_weight']*time_fraction/attrs['total_generated'])
        bank_info[bank_name]=dict(accepted_rows=len(union),union_found=int(union.sum()),matched_found=int(matched.sum()),total_generated=attrs['total_generated'])
        for policy,summary in stats.items():out[bank_name+':'+policy]=dict(early=early_stats,o3=summary)
    return out,time_fraction,bank_info


def support_summary(time_fraction,all_selections):
    zgrid=model.ZGRID
    def z_integral(low,high,kappa):
        f=lambda z:np.interp(z,zgrid,model.ZDV)*(1+z)**(kappa-1)
        return quad(f,low,high,points=zgrid[(zgrid>low)&(zgrid<high)],limit=4000,epsabs=.1,epsrel=1e-10)[0]
    rows=[]
    for m in MODELS:
        kappa=m['kappa'];total=z_integral(0,2.3,kappa);tail=z_integral(1.9,2.3,kappa)/total
        pm=lambda x:float(mass_pdf(np.asarray(x),m['alpha'],m['peak_fraction']))
        added=quad(lambda x:pm(x)*(1-(10/x)**2)/(1-(3/x)**2),90,100,epsabs=1e-13)[0]
        aa,bb=m['magnitude'];spin=betainc(aa,bb,.998)**2
        pure_missing=1-spin*(1-tail)
        mixture_missing=1-spin*(1-tail+tail*added)
        per_bank={}
        for bank,missing in [('pure_bbh',pure_missing),('official_mixture',mixture_missing)]:
            # Same missing source mass in either normalized sign component.
            upper=time_fraction*missing
            minimum=min(float(all_selections[bank+':'+policy][m['index']].min()) for policy in ['matched_filter','with_cwb'])
            per_bank[bank]=dict(uncovered_target_probability=float(missing),o3_detection_contribution_upper=float(upper),
                               upper_over_minimum_represented_component=float(upper/minimum))
        rows.append(dict(index=m['index'],redshift_tail_probability=float(tail),high_z_imbh_mass_fraction=float(added),
                         coded_redshift_integral_over_trapezoid_normalization=float(total/model.ZNORM[kappa]),banks=per_bank))
    return dict(models=rows,scope='Normalized target proposal support before hopeless cuts. Selection upper bounds use only p_det<=1; they are not actual omitted contribution measurements.')


def main():
    init_model()
    design=json.loads((BASE/'design.json').read_text())
    assert len(MODELS)==design['models']==54
    for validation,result_path in [(INV/'control-validation.json',INV/'control-results.json'),
                                   (CAT/'native-catalog-validation.json',CAT/'native-catalog-results.json'),
                                   (MIX/'validation.json',MIX/'results.json')]:
        assert json.loads(validation.read_text())['result_sha256']==sha(result_path),validation
    events_used,protocol=events()
    means=np.array([r['mean'] for r in events_used]);second=np.array([r['second'] for r in events_used])
    parts,time_fraction,bank_info=selections()
    selections_by_case={name:rec['early']['sum']+rec['o3']['sum'] for name,rec in parts.items()}
    g=np.linspace(0,1,4001);cases={};saved=dict(g=g,means=means,second=second)
    baseline=json.loads((MIX/'results.json').read_text())
    selection_bridges=[]
    for name,selection in selections_by_case.items():
        _,q,p=posterior(means,selection,len(g))
        _,fine,_=posterior(means,selection,2*len(g)-1)
        error=float(np.max(abs(q-fine)));assert error<2e-6
        bank,policy=name.split(':')
        old_case='pure_bbh_baseline' if bank=='pure_bbh' else 'mixture_full_union'
        expected=np.array(baseline['metrics']['g']['scenarios'][policy][old_case]['quantiles'])
        bridge=float(np.max(abs(q[6:12]-expected)));assert bridge<2e-6
        with np.load(MIX/'study-checkpoints.npz',allow_pickle=False) as oldcp:
            oldsel=oldcp['g_'+policy+'_'+old_case+'_selection']
        selection_error=relative(selection[6:12],oldsel);assert selection_error<1e-10
        selection_bridges.append(selection_error)
        contrasts=[]
        for mi in range(3):
            for zi in range(3):
                start=mi*18+zi*6
                for li,shape in enumerate([0,4,16]):
                    contrasts.append(dict(mass_id=MODELS[start]['mass_id'],kappa=MODELS[start]['kappa'],shape=shape,
                                          value_pp=float(100*(q[start+3+li,1]-q[start+li,1]))))
        mean_mixture=q[None,:,1]*means[:,:,0]+(1-q[None,:,1])*means[:,:,1]
        moment2=q[None,:,1]**2*second[:,:,0,0]+2*q[None,:,1]*(1-q[None,:,1])*second[:,:,0,1]+(1-q[None,:,1])**2*second[:,:,1,1]
        n=np.array([r['rows'] for r in events_used])[:,None]
        event_ess=np.divide(n*mean_mixture**2,moment2,out=np.zeros_like(moment2),where=moment2>0)
        pooled_outer=parts[name]['early']['outer']+parts[name]['o3']['outer']
        contribution=q[:,1]*selection[:,0]+(1-q[:,1])*selection[:,1]
        weighted_second=q[:,1]**2*pooled_outer[:,0,0]+2*q[:,1]*(1-q[:,1])*pooled_outer[:,0,1]+(1-q[:,1])**2*pooled_outer[:,1,1]
        cases[name]=dict(quantiles=q.tolist(),grid_doubling_maximum_quantile_error=error,baseline_quantile_bridge=bridge,
                         magnitude_contrasts=contrasts,smallest_magnitude_contrast_pp=min(x['value_pp'] for x in contrasts),
                         largest_magnitude_contrast_pp=max(x['value_pp'] for x in contrasts),
                         nonpositive_magnitude_contrast_count=sum(x['value_pp']<=0 for x in contrasts),
                         event_minimum_weight_ess_by_model=event_ess.min(axis=0).tolist(),
                         pooled_selection_weight_concentration_ess=(contribution**2/weighted_second).tolist())
        centers=q[:,1].reshape((3,3,2,3))
        cases[name]['median_sensitivity_pp']=dict(mass_only_maximum_span=float(100*np.ptp(centers,axis=0).max()),
            redshift_only_maximum_span=float(100*np.ptp(centers,axis=1).max()),
            joint_mass_redshift_maximum_span=float(100*(centers.max(axis=(0,1))-centers.min(axis=(0,1))).max()),
            scope='Maxima over declared conditional model choices, not uncertainty estimates or nuisance-parameter posteriors.')
        saved[name+'_selection']=selection;saved[name+'_density']=p
        saved[name+'_early']=parts[name]['early']['sum'];saved[name+'_o3']=parts[name]['o3']['sum']
    bank_comparisons={}
    for policy in ['matched_filter','with_cwb']:
        pure='pure_bbh:'+policy;mixture='official_mixture:'+policy
        q1=np.array(cases[pure]['quantiles']);q2=np.array(cases[mixture]['quantiles'])
        tv=.5*np.trapezoid(abs(saved[mixture+'_density']-saved[pure+'_density']),g,axis=1)
        bank_comparisons[policy]=dict(quantile_shift_pp=(100*(q2-q1)).tolist(),posterior_total_variation=tv.tolist(),
            maximum_median_shift_pp=float(np.max(abs(100*(q2[:,1]-q1[:,1])))),maximum_tv=float(tv.max()))
    envelopes=[]
    for name,case in cases.items():
        q=np.array(case['quantiles']).reshape((3,3,2,3,3))
        for si,mag in enumerate([[1.5,8.5],[2,5]]):
            for li,shape in enumerate([0,4,16]):
                allq=q[:,:,si,li,:]
                envelopes.append(dict(case=name,magnitude=mag,shape=shape,median_minimum=float(allq[:,:,1].min()),
                                      median_maximum=float(allq[:,:,1].max()),quantile5_minimum=float(allq[:,:,0].min()),
                                      quantile95_maximum=float(allq[:,:,2].max()),scope='Declared3x3 model envelope; not a credible interval, fitted nuisance posterior or Bayesian model average.'))
    report=dict(schema='blackhole-full69-mass-redshift-grid/1',event_count=69,source_rows=sum(r['rows'] for r in events_used),
                model_count=54,conditional_fits=216,functional=design['functional'],models=MODELS,prior=design['prior'],
                event_protocol=protocol,cases=cases,bank_comparisons=bank_comparisons,model_envelopes=envelopes,
                o3_time_fraction=time_fraction,bank_metadata=bank_info,support=support_summary(time_fraction,selections_by_case),
                maximum_baseline_event_relative_error=max(r['baseline_relative_error'] for r in events_used),
                maximum_baseline_selection_relative_error=max(selection_bridges),
                exact_zero_event_components=[dict(index=r['index'],event=r['event'],model_components=np.argwhere(np.array(r['mean'])==0).tolist()) for r in events_used if np.any(np.array(r['mean'])==0)],
                weight_diagnostics_scope='Contribution concentration only; no chain/campaign correlation calibration, no ESS-based rejection.',
                limits=design['limits'],input_sha256={})
    report['runtime']=dict(python=sys.version.split()[0],packages={name:version(name) for name in ['numpy','scipy','astropy','h5py']})
    paths=[BASE/'design.json',BASE/'grid_model.py',Path(__file__),BASE/'event-integrals.json',ROOT/'model.py',INV/'run_control.py',INV/'invariant_model.py',INV/'event-integrals.json',INV/'control-checkpoints.npz',INV/'control-results.json',INV/'control-validation.json',CAT/'event-provenance.json',CAT/'native-catalog-results.json',CAT/'native-catalog-validation.json',ROOT/'research/population-stage5/run_tail_study.py',ROOT/'research/population-stage4/recover_prior.py',SEL/'run_selection_study.py',SEL/'o3-provenance.json',SEL/'o3-bank.npz',MIX/'download.json',MIX/'results.json',MIX/'validation.json']
    for path in paths:report['input_sha256'][path.relative_to(ROOT).as_posix()]=sha(path)
    np.savez_compressed(BASE/'checkpoints.npz',**saved)
    report['checkpoint_sha256']=sha(BASE/'checkpoints.npz')
    write_json(BASE/'results.json',report)
    print(json.dumps(dict(conditional_fits=216,contrast_ranges={name:[case['smallest_magnitude_contrast_pp'],case['largest_magnitude_contrast_pp']] for name,case in cases.items()},bank_comparisons={k:{a:b for a,b in v.items() if a.startswith('maximum')} for k,v in bank_comparisons.items()}),indent=2))


if __name__=='__main__':main()
