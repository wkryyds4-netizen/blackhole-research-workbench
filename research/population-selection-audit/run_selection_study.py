"""Full-O3 representation and joint event/selection redshift-support controls."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
from scipy.integrate import cumulative_trapezoid, quad
from astropy.cosmology import Planck15

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
sys.path.insert(0,str(ROOT/'research/population-stage5'))
from run_tail_study import init_model,weights,coordinates,density_pair,KEYS,PARAMS,LAL


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def physical(raw):
    # Combined file stores float64 copies of original float32 physical values.
    # Promote before ratios/norms so identical source values use identical arithmetic.
    m=raw['mass1_source'].astype(float)
    d=dict(mass_1=m,mass_ratio=raw['mass2_source'].astype(float)/m,redshift=raw['redshift'].astype(float))
    for i in [1,2]:
        a=np.sqrt(sum(raw[f'spin{i}{axis}'].astype(float)**2 for axis in 'xyz'))
        d[f'a_{i}']=a;d[f'cos_tilt_{i}']=raw[f'spin{i}z'].astype(float)/a
    d['prior']=raw['sampling_pdf']*m*(2*np.pi*d['a_1']**2)*(2*np.pi*d['a_2']**2)
    return d


def posterior(means,selection,points=1001):
    r,fs=coordinates(points)
    logp=np.log(fs[None,:,:]*means[:,:,0,None]+(1-fs[None,:,:])*means[:,:,1,None]).sum(axis=0)
    logp-=len(means)*np.log(fs*selection[:,0,None]+(1-fs)*selection[:,1,None])
    p=np.exp(logp-logp.max(axis=1)[:,None]);p/=np.trapezoid(p,r,axis=1)[:,None]
    c=cumulative_trapezoid(p,r,axis=1,initial=0);c/=c[:,-1,None]
    q=np.array([np.interp([.05,.5,.95],v,r) for v in c])
    return q,p


def relative(a,b):
    v=np.abs(a/b-1)
    return dict(median=float(np.median(v)),p99=float(np.quantile(v,.99)),maximum=float(v.max()))


def proposal_audit(raw,attrs):
    z=raw['redshift'].astype(float);grid=np.linspace(0,attrs['max_z'],attrs['num_redshift_norm_grid_points'])
    dv=LAL.differential_comoving_volume(grid).value*(1+grid)**(attrs['pow_z']-1)
    norm=np.trapezoid(dv,grid)
    interp=np.interp(z,grid,dv)/norm
    direct=LAL.differential_comoving_volume(z).value*(1+z)**(attrs['pow_z']-1)/norm
    zcheck=relative(interp,raw['redshift_sampling_pdf']);assert zcheck['maximum']<1e-9
    m=raw['mass1_source'].astype(float);n=raw['mass2_source'].astype(float)
    slope=attrs['pow_mass1'];p1=m**slope/((attrs['max_mass1']**(slope+1)-attrs['min_mass1']**(slope+1))/(slope+1))
    p2=2*n/(m*m-attrs['min_mass2']**2)
    factors={'mass1':relative(p1,raw['mass1_source_sampling_pdf']),
             'joint_mass':relative(p1*p2,raw['mass1_source_mass2_source_sampling_pdf'])}
    for i in [1,2]:
        mag2=sum(raw[f'spin{i}{ax}'].astype(float)**2 for ax in 'xyz')
        p=1/(4*np.pi*attrs[f'max_spin{i}']*mag2)
        factors[f'spin{i}']=relative(p,raw[f'spin{i}x_spin{i}y_spin{i}z_sampling_pdf'])
    product=raw['mass1_source_mass2_source_sampling_pdf'].astype(float)*raw['redshift_sampling_pdf']*raw['spin1x_spin1y_spin1z_sampling_pdf']*raw['spin2x_spin2y_spin2z_sampling_pdf']
    productcheck=relative(product,raw['sampling_pdf']);assert productcheck['maximum']<1e-6
    return dict(redshift_grid_interpolated=zcheck,redshift_continuous_instead_of_stored_grid=relative(direct,raw['redshift_sampling_pdf']),
        formula_factors_on_rounded_coordinates=factors,stored_factor_product=productcheck,
        note='Redshift interpolation reproduces recorded density. Float32 masses/spins prevent exact inversion of all original factors; near the lower conditional-mass boundary relative errors amplify. Calculations retain the official joint sampling_pdf.')


def main():
    init_model();frozen=json.loads((ROOT/'results/model.json').read_text())
    stage7=json.loads((ROOT/'research/population-stage7/comparison-results.json').read_text())
    source=json.loads((BASE/'o3-provenance.json').read_text());assert digest(BASE/'o3-bank.npz')==source['cache_sha256']
    for path,sha in stage7['frozen_hashes'].items():assert digest(ROOT/path)==sha
    with np.load(BASE/'o3-bank.npz',allow_pickle=False) as bank:
        full={k[5:]:bank[k] for k in bank.files if k.startswith('full_')}
        old={k[9:]:bank[k] for k in bank.files if k.startswith('combined_') and k!='combined_o3_to_full'}
    proposal=proposal_audit(full,source['attributes'])
    local=dict(np.load(ROOT/'data/injections.npz',allow_pickle=False));total=float(local['total_generated'])
    found_old=np.any(np.stack([old[k]>1 for k in old if k.startswith('ifar_')]),axis=0)|((old['name']!=b'o3')&(old['optimal_snr_net']>10))
    found_full=np.any(np.stack([full[k]>1 for k in full if k.startswith('ifar_')]),axis=0)
    old_d=physical({k:v[found_old] for k,v in old.items()});full_d=physical({k:v[found_full] for k,v in full.items()})
    for k in PARAMS+['prior']:assert np.array_equal(old_d[k],local[k])
    old_w=weights(old_d);full_w=weights(full_d)
    nominal=np.array([frozen['scenarios'][k]['selection'][:2] for k in KEYS])
    baseline=old_w.sum(axis=0)/total
    baseline_error=float(np.max(abs(baseline/nominal-1)));assert baseline_error<1e-10
    names=old['name'][found_old];mw=old['mixture_weight'][found_old]
    attrs=source['attributes'];time_fraction=attrs['analysis_time_s']/source['combined_attributes']['analysis_time_s']
    old_total=source['combined_attributes']['total_generated']
    by_run={run:(old_w[names==run.encode()]*mw[names==run.encode(),None,None]).sum(axis=0)/old_total for run in ['o1','o2','o3']}
    old_exact=sum(by_run.values());full_o3=full_w.sum(axis=0)*time_fraction/attrs['total_generated']
    new=by_run['o1']+by_run['o2']+full_o3
    cap_old=old_d['redshift']<=1.9;cap_full=full_d['redshift']<=1.9;assert cap_full.all()
    cap_selection=(old_w[(names!=b'o3')&cap_old]*mw[(names!=b'o3')&cap_old,None,None]).sum(axis=0)/old_total+full_o3
    original=np.array([[[*frozen['scenarios'][key]['event_moments'][i][:2]] for key in KEYS] for i in range(69)])
    with np.load(ROOT/'research/population-stage7/comparison-checkpoints.npz',allow_pickle=False) as cp:
        full5=cp['all_event_means'][cp['names'].tolist().index('new_plus_38_nocosmo')]
    cap_tutorial=[];cap_five=[];tails=[]
    with np.load(ROOT/'data/posteriors.npz',allow_pickle=False) as archive:
        for pos,event in enumerate(frozen['events']):
            ix=event['index'];d={k:archive[f'{ix}_{k}'] for k in PARAMS+['prior']};w=weights(d)
            assert np.max(abs(w.mean(axis=0)/original[pos]-1))<1e-10
            cap=(w*(d['redshift']<=1.9)[:,None,None]).mean(axis=0);cap_tutorial.append(cap)
            if ix in [19,20,25,38,50]:
                path=ROOT/('research/population-stage5/representations.npz' if ix==38 else f'research/population-stage7/event-{ix}.npz')
                with np.load(path,allow_pickle=False) as bank:d={k:bank['nocosmo_'+k] for k in PARAMS}
                _,prior,_=density_pair(d);w=weights(dict(d,prior=prior))
                assert np.max(abs(w.mean(axis=0)/full5[pos]-1))<1e-10
            cap=(w*(d['redshift']<=1.9)[:,None,None]).mean(axis=0);cap_five.append(cap)
            tails.append(dict(index=ix,event=event.get('name',str(ix)),rows_above_1p9=int((d['redshift']>1.9).sum()),
                component_integral_fraction_above_1p9=((w*(d['redshift']>1.9)[:,None,None]).sum(axis=0)/w.sum(axis=0)).tolist()))
    cap_tutorial=np.array(cap_tutorial);cap_five=np.array(cap_five)
    cases={'tutorial_frozen':(original,baseline),'tutorial_exact_mixture':(original,old_exact),
           'tutorial_full_o3':(original,new),'five_frozen':(full5,baseline),
           'five_exact_mixture':(full5,old_exact),'five_full_o3':(full5,new),
           'tutorial_full_o3_zcap':(cap_tutorial,cap_selection),'five_full_o3_zcap':(cap_five,cap_selection)}
    redshift_mass={}
    for lam in [0,2,4]:
        fn=lambda z:Planck15.differential_comoving_volume(z).value*(1+z)**(lam-1)
        redshift_mass[str(lam)]=quad(fn,1.9,2.3)[0]/quad(fn,0,2.3)[0]
    # All six settings use lambda=2. Its normalization cancels between N events and xi^N.
    pcap=1-redshift_mass['2'];cases['five_zcap_normalized']=(cap_five/pcap,cap_selection/pcap)
    summaries={};curves=[];max_grid=0.;r,fs=coordinates()
    qref,_=posterior(original,baseline)
    assert np.max(abs(qref-np.array(stage7['scenarios']['original']['quantiles'])))<1e-11
    for name,(means,selection) in cases.items():
        q,p=posterior(means,selection);fine,_=posterior(means,selection,2001)
        max_grid=max(max_grid,float(np.max(abs(q-fine))))
        summaries[name]=dict(quantiles=q.tolist(),median_shift_from_frozen_pp=(100*(q[:,1]-qref[:,1])).tolist(),
            sigma03_minus12_gap_pp=[float(100*(q[0,1]-q[2,1])),float(100*(q[3,1]-q[5,1]))])
        curves.append(p)
    assert max_grid<1e-6
    error=float(np.max(abs(np.array(summaries['five_zcap_normalized']['quantiles'])-summaries['five_full_o3_zcap']['quantiles'])));assert error<1e-11
    q5=np.array(summaries['five_frozen']['quantiles']);qnew=np.array(summaries['five_full_o3']['quantiles']);qcap=np.array(summaries['five_full_o3_zcap']['quantiles'])
    selected_diag=[]
    for k,key in enumerate(KEYS):
        f=np.interp(q5[k,1],r,fs[k]);a=old_w[names==b'o3',k];b=full_w[:,k]
        w1=f*a[:,0]+(1-f)*a[:,1];w2=f*b[:,0]+(1-f)*b[:,1]
        selected_diag.append(dict(key=key,o3_combined_ess=float(w1.sum()**2/np.sum(w1*w1)),
            o3_full_ess=float(w2.sum()**2/np.sum(w2*w2)),
            old_o3_fraction=float((f*by_run['o3'][k,0]+(1-f)*by_run['o3'][k,1])/(f*old_exact[k,0]+(1-f)*old_exact[k,1])),
            full_over_exact_selection=float((f*new[k,0]+(1-f)*new[k,1])/(f*old_exact[k,0]+(1-f)*old_exact[k,1]))))
    result=dict(schema='blackhole-full-o3-selection-study/1',keys=KEYS,scenarios=summaries,proposal_audit=proposal,
        n_events=69,redshift_cut=1.9,population_redshift_mass_1p9_to_2p3=redshift_mass,
        o3_time_fraction=time_fraction,baseline_selection_relative_error=baseline_error,
        normalized_cap_quantile_error=error,max_grid_quantile_error=max_grid,
        five_event_full_o3_shift_pp=(100*(qnew[:,1]-q5[:,1])).tolist(),
        five_event_redshift_cap_shift_pp=(100*(qcap[:,1]-qnew[:,1])).tolist(),
        selection_diagnostics=selected_diag,selection_components={k:v.tolist() for k,v in dict(frozen=baseline,exact_mixture=old_exact,full_o3=new,zcap=cap_selection,**by_run).items()},
        event_tail_diagnostics=tails,frozen_hashes=stage7['frozen_hashes'],o3_bank_sha256=source['cache_sha256'],
        limits='Related finite injection representations; native spin references and original search rule retained. A joint z cap changes the target population. No measurement of detection efficiency above the O3 proposal support, no calibrated MC error or formation-channel fractions.')
    (BASE/'selection-results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    np.savez_compressed(BASE/'selection-checkpoints.npz',r=r,fs=fs,names=np.array(list(cases)),
        means=np.array([v[0] for v in cases.values()]),selection=np.array([v[1] for v in cases.values()]),densities=np.array(curves))
    print(json.dumps({k:result[k] for k in ['five_event_full_o3_shift_pp','five_event_redshift_cap_shift_pp','population_redshift_mass_1p9_to_2p3','max_grid_quantile_error','selection_diagnostics']},indent=2))


if __name__=='__main__':main()
