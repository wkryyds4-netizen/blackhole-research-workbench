"""Full-source, fixed69 control for the negative-effective-spin fraction.

No fit is permitted until all original-prior event sources are restored and
validated. Native tilt results and their frozen inputs are never overwritten.
"""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
from scipy.integrate import cumulative_trapezoid
from invariant_model import KEYS,MAGNITUDES,LAMBDAS,component_factors,effective_spin,symmetric_ratio

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
CAT=ROOT/'research/population-catalog'
SEL=ROOT/'research/population-selection-audit'
sys.path.insert(0,str(SEL))
from run_selection_study import init_model,density_pair,physical,PARAMS
import model


def digest(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def weights(d):
    result=[]
    for spin,_ in enumerate(MAGNITUDES):
        base=model.base_density(d,spin,0,2)/d['prior']/4
        for lam in LAMBDAS:
            result.append(base[:,None]*component_factors(d['mass_ratio'],d['a_1'],d['a_2'],d['cos_tilt_1'],d['cos_tilt_2'],lam))
    out=np.stack(result,axis=1)
    assert np.isfinite(out).all() and np.all(out>=0)
    return out


def posterior(means,selection,points):
    g=np.linspace(0,1,points)
    mix=g*means[:,:,0,None]+(1-g)*means[:,:,1,None]
    efficiency=g*selection[:,0,None]+(1-g)*selection[:,1,None]
    assert np.all(efficiency>0)
    with np.errstate(divide='ignore'):
        logs=np.log(mix).sum(axis=0)-len(means)*np.log(efficiency)
    p=np.exp(logs-logs.max(axis=1)[:,None]);p/=np.trapezoid(p,g,axis=1)[:,None]
    c=cumulative_trapezoid(p,g,axis=1,initial=0);c/=c[:,-1,None]
    q=np.array([np.interp([.05,.5,.95],v,g) for v in c])
    return g,q,p


def main():
    init_model()
    native_path=CAT/'native-catalog-results.json';native=json.loads(native_path.read_text())
    validated=json.loads((CAT/'native-catalog-validation.json').read_text())
    assert native['completed_full_source_events']==69,'Restore and validate all69 original-prior sources before fitting this control.'
    assert validated['result_sha256']==digest(native_path)
    frozen=json.loads((ROOT/'results/model.json').read_text())
    provenance=json.loads((CAT/'event-provenance.json').read_text())
    old_integrals=json.loads((CAT/'native-event-integrals.json').read_text())
    assert {r['index'] for r in provenance if r.get('native_source_identity_resolved')}==set(native['full_source_indices_in_likelihood'])
    angular=json.loads((BASE/'angular-validation.json').read_text());assert angular['model_sha256']==digest(BASE/'invariant_model.py')
    protocol=hashlib.sha256(b'invariant-control/1'+(BASE/'invariant_model.py').read_bytes()+(BASE/'run_control.py').read_bytes()+(ROOT/'model.py').read_bytes()+(ROOT/'research/population-stage5/run_tail_study.py').read_bytes()).hexdigest()
    cache_path=BASE/'event-integrals.json'
    cache=json.loads(cache_path.read_text()) if cache_path.exists() else dict(protocol=protocol,events=[])
    assert cache['protocol']==protocol,'Preserve/version existing integrals before changing their protocol.'
    means=[];details=[]
    for event in frozen['events']:
        ix=event['index'];info=next(r for r in provenance if r['index']==ix)
        path=CAT/f'event-{ix}.npz';assert digest(path)==info['cache_sha256']
        saved=next((r for r in cache['events'] if r['index']==ix and r['source_cache_sha256']==info['cache_sha256']),None)
        if saved is None:
            with np.load(path,allow_pickle=False) as archive:d={k:archive[k] for k in archive.files}
            _,prior,dl=density_pair(d);d['prior']=prior
            w=weights(d);mean=w.mean(axis=0);second=np.einsum('nki,nkj->kij',w,w)/len(w)
            old=next(r for r in old_integrals['events'] if r['index']==ix)
            iso=(mean[[0,3],0]+mean[[0,3],1])/2
            bridge=float(np.max(abs(iso/np.array(old['mean'])[[0,3],1]-1)));assert bridge<1e-10
            x=effective_spin(d['mass_ratio'],d['a_1'],d['a_2'],d['cos_tilt_1'],d['cos_tilt_2'])
            xinf=effective_spin(d['mass_ratio'],d['a_1'],d['a_2'],d['cos_tilt_1_infinity_only_prec_avg'],d['cos_tilt_2_infinity_only_prec_avg'])
            finite=np.isfinite(xinf);transport=[];full_transport=[]
            for lam in LAMBDAS:
                a=symmetric_ratio(d['mass_ratio'][finite],d['a_1'][finite],d['a_2'][finite],x[finite],lam)
                b=symmetric_ratio(d['mass_ratio'][finite],d['a_1'][finite],d['a_2'][finite],xinf[finite],lam)
                transport.append(float(np.max(abs(a-b))) if finite.any() else None)
                errors=[]
                for test_g in [.1,.5,.9]:
                    split=lambda t:2*test_g*(t<0)+2*(1-test_g)*(t>0)+(t==0)
                    errors.append(float(np.max(abs(a*split(x[finite])-b*split(xinf[finite])))) if finite.any() else None)
                full_transport.append(errors)
            saved=dict(index=ix,event=info['event'],rows=len(w),source_cache_sha256=info['cache_sha256'],mean=mean.tolist(),second=second.tolist(),
                component_maximum=w.max(axis=0).tolist(),isotropic_bridge_relative_error=bridge,
                exact_zero_effective_spin_rows=int(np.count_nonzero(x==0)),
                released_infinity_finite_rows=int(finite.sum()),released_infinity_nonfinite_rows=int((~finite).sum()),
                maximum_released_effective_spin_difference=float(np.max(abs(x[finite]-xinf[finite]))) if finite.any() else None,
                maximum_symmetric_ratio_transport_difference=transport,
                maximum_full_ratio_transport_difference_at_g01_05_09=full_transport,
                released_infinity_sign_disagreement_rows=int(np.count_nonzero(np.sign(x[finite])!=np.sign(xinf[finite]))),
                limits='Released-derived-column conservation check only; not a proof of invariant PN angular measure or reconstructed historical injection reference.')
            cache['events']=[r for r in cache['events'] if r['index']!=ix]+[saved]
            cache_path.write_text(json.dumps(cache,indent=2),encoding='utf-8')
            print(f'Recomputed invariant event {ix}: {len(w)} rows',flush=True)
        means.append(saved['mean']);details.append(saved)
    means=np.array(means);assert means.shape==(69,6,2)
    source=json.loads((SEL/'o3-provenance.json').read_text());assert digest(SEL/'o3-bank.npz')==source['cache_sha256']
    previous=json.loads((SEL/'selection-results.json').read_text());policy=json.loads((SEL/'pipeline-results.json').read_text())
    with np.load(SEL/'o3-bank.npz',allow_pickle=False) as bank:
        full={k[5:]:bank[k] for k in bank.files if k.startswith('full_')}
        old={k[9:]:bank[k] for k in bank.files if k.startswith('combined_') and k!='combined_o3_to_full'}
    search_keys=sorted(k for k in full if k.startswith('ifar_'));matched=[k for k in search_keys if k!='ifar_cwb'];assert len(matched)==4
    union=np.any([full[k]>1 for k in search_keys],axis=0);mf=np.any([full[k]>1 for k in matched],axis=0)
    early=(old['name']!=b'o3')&(np.any([old[k]>1 for k in search_keys],axis=0)|(old['optimal_snr_net']>10))
    eweights=weights(physical({k:v[early] for k,v in old.items()}))
    eweights*=old['mixture_weight'][early,None,None]/source['combined_attributes']['total_generated']
    oweights=weights(physical({k:v[union] for k,v in full.items()}))
    factor=previous['o3_time_fraction']/source['attributes']['total_generated']
    early_sum=eweights.sum(axis=0);all_sum=oweights.sum(axis=0)*factor;mf_sum=oweights[mf[union]].sum(axis=0)*factor
    selections=dict(matched_filter=early_sum+mf_sum,with_cwb=early_sum+all_sum,without_selection_diagnostic=np.ones_like(early_sum))
    bridges={}
    for name,sel in selections.items():
        if name=='without_selection_diagnostic':continue
        oldsel=np.array(policy['selection_components']['matched_filter' if name=='matched_filter' else 'union'])
        error=float(np.max(abs((sel[[0,3],0]+sel[[0,3],1])/2/oldsel[[0,3],1]-1)))
        assert error<1e-10;bridges[name]=error
    summaries={};curves=[];grid_details={}
    for name,sel in selections.items():
        g,coarse,p=posterior(means,sel,1001);checks=[]
        for points in [2001,4001,8001,16001]:
            fine_g,fine,fine_p=posterior(means,sel,points);error=float(np.max(abs(fine-coarse)))
            checks.append(dict(points=points,maximum_quantile_change=error))
            if error<2e-6:break
            coarse=fine
        else:raise ValueError('Quantile grid did not meet the fixed2e-6 absolute tolerance')
        summaries[name]=dict(quantiles=fine.tolist(),median_negative_probability=fine[:,1].tolist(),
            magnitude_model_median_difference_pp=(100*(fine[3:,1]-fine[:3,1])).tolist(),
            lambda16_minus0_median_difference_pp=[float(100*(fine[2,1]-fine[0,1])),float(100*(fine[5,1]-fine[3,1]))])
        grid_details[name]=checks
    # Save curves on one common grid for independent product verification.
    points=max(c[-1]['points'] for c in grid_details.values())
    for name,sel in selections.items():
        g,q,p=posterior(means,sel,points)
        summaries[name].update(quantiles=q.tolist(),median_negative_probability=q[:,1].tolist(),
            magnitude_model_median_difference_pp=(100*(q[3:,1]-q[:3,1])).tolist(),
            lambda16_minus0_median_difference_pp=[float(100*(q[2,1]-q[0,1])),float(100*(q[5,1]-q[3,1]))])
        curves.append(p)
    q=np.array(summaries['matched_filter']['quantiles'])
    concentration=[]
    for row in details:
        mu=np.array(row['mean']);second=np.array(row['second']);f=q[:,1]
        m=f*mu[:,0]+(1-f)*mu[:,1]
        m2=f*f*second[:,0,0]+2*f*(1-f)*second[:,0,1]+(1-f)**2*second[:,1,1]
        concentration.append(dict(index=row['index'],event=row['event'],rows=row['rows'],ess_at_medians=(row['rows']*m*m/m2).tolist(),
            limits='Contribution ESS only; excludes PE-chain correlations and actual-campaign design. No ESS-based parameter or event rejection.'))
    report=dict(schema='blackhole-reference-invariant-control/1',keys=KEYS,event_count=69,full_source_rows=sum(r['rows'] for r in details),
        functional='g=P(chi_eff<0), binary-level; not primary negative tilt or a formation-channel fraction',prior=dict(family='Uniform',minimum=0.,maximum=1.),
        approximation='2PN conserved effective spin and preserved isotropic angular measure; not exact full-GR invariance.',
        scenarios=summaries,grid_checks=grid_details,common_grid_points=points,isotropic_selection_bridges=bridges,event_concentration=concentration,
        cwb_minus_matched_median_pp=(100*(np.array(summaries['with_cwb']['quantiles'])[:,1]-q[:,1])).tolist(),
        no_selection_minus_matched_median_pp=(100*(np.array(summaries['without_selection_diagnostic']['quantiles'])[:,1]-q[:,1])).tolist(),
        input_hashes=dict(native_results=digest(native_path),native_validation=digest(CAT/'native-catalog-validation.json'),design=digest(BASE/'design.json'),model=digest(BASE/'invariant_model.py'),event_integrals=digest(cache_path)),
        limits='Restricted sign-shape families. Unknown injection execution reference, prior-excluded likelihood, O3 redshift support, early semianalytic approximations, finite-bank error and PE-chain correlation are not repaired. Does not replace original native-tilt audit or identify astrophysical channels.')
    (BASE/'control-results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    np.savez_compressed(BASE/'control-checkpoints.npz',g=g,means=means,selection_matched=selections['matched_filter'],selection_union=selections['with_cwb'],selection_none=selections['without_selection_diagnostic'],densities=np.array(curves),
        early_weights=eweights,o3_union_weights=oweights,o3_matched_mask=mf[union],o3_factor=factor,zgrid=model.ZGRID,zdv=model.ZDV,znorm=model.ZNORM[2])
    print(json.dumps(dict(events=69,rows=report['full_source_rows'],scenarios=summaries,cwb_minus_matched_median_pp=report['cwb_minus_matched_median_pp'],grid_checks=grid_details),indent=2))


if __name__=='__main__':main()
