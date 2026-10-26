"""Replace only O3 sensitivity with the official four-population mixture.

Calls established scientific density functions; never runs model.py as a script
or rewrites existing fits, event integrals, or their caches.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import hashlib
import json
import sys
import h5py
import numpy as np
from scipy.integrate import quad
from scipy.special import betainc, ndtr

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
SEL = ROOT/'research/population-selection-audit'
INV = ROOT/'research/population-invariant'
sys.path.insert(0,str(SEL))
sys.path.insert(0,str(INV))
from run_selection_study import init_model, physical, weights as native_weights, posterior as native_posterior
from run_control import weights as invariant_weights, posterior as invariant_posterior
import model


def sha(path):
    with path.open('rb') as stream: return hashlib.file_digest(stream,'sha256').hexdigest()


def relative(a,b):
    assert np.array_equal(a == 0,b == 0)
    active = b != 0
    return float(np.max(abs(a[active]/b[active]-1))) if active.any() else 0.


def summarize(metric, means, selection, points=4001):
    if metric == 'g':
        grid,q,density = invariant_posterior(means,selection,points)
        _,fine,_ = invariant_posterior(means,selection,2*points-1)
        contrasts = 100*(q[3:,1]-q[:3,1])
    else:
        q,density = native_posterior(means,selection,points)
        fine,_ = native_posterior(means,selection,2*points-1)
        from run_selection_study import coordinates
        grid,_ = coordinates(points)
        contrasts = 100*(q[[0,3],1]-q[[2,5],1])
    error = float(np.max(abs(q-fine)))
    assert error < 2e-6, error
    return dict(quantiles=q.tolist(),model_median_contrasts_pp=contrasts.tolist(),
                max_grid_quantile_error=error), grid, density


def support_mass():
    # Union within our m1<=100,m2>=3 target: BBH below z=1.9;
    # IMBH can add only m1>=90,m2>=10 through z=2.5. BNS/NSBH add none.
    a,b=3.,100.
    gaussian_norm=5*np.sqrt(2*np.pi)*(ndtr(13)-ndtr(-6.4))
    mass=lambda m: .9*m**-2/(1/a-1/b)+.1*np.exp(-.5*((m-35)/5)**2)/gaussian_norm
    added_mass=quad(lambda m: mass(m)*(1-(10/m)**2)/(1-(3/m)**2),90,100,epsabs=1e-13)[0]
    tail=quad(lambda z: np.interp(z,model.ZGRID,model.ZDV)*(1+z),1.9,2.3,epsabs=1e-3)[0]/model.ZNORM[2]
    # This mass/shape-independent fraction applies to the symmetric total spin
    # model and each normalized sign component; isotropic sign symmetry holds.
    result=[]
    for aa,bb in [(1.5,8.5),(2.,5.)]:
        spin_supported=float(betainc(aa,bb,.998)**2)
        old=spin_supported*(1-tail)
        new=spin_supported*((1-tail)+tail*added_mass)
        result.append(dict(magnitude=[aa,bb],old_bbh_covered_target_mass=old,
                           mixture_union_covered_target_mass=new,
                           mixture_added_target_mass=new-old,
                           still_uncovered_target_mass=1-new))
    return dict(target_redshift_tail_mass=tail,conditional_mass_fraction_covered_by_imbh_at_high_z=added_mass,
                magnitude_families=result,
                interpretation='Proposal support only, before hopeless cuts and detection. A zero saved-row tail contribution is not a detection-probability bound.')


def main():
    init_model()
    provenance=json.loads((BASE/'provenance.json').read_text())
    path=BASE/'sources'/provenance['source']['file']
    assert sha(path)==provenance['source']['sha256']
    with np.load(INV/'control-checkpoints.npz',allow_pickle=False) as cp:
        gmeans=cp['means']; gearly=cp['early_weights'].sum(axis=0)
        gbaseline={'matched_filter':cp['selection_matched'],'with_cwb':cp['selection_union']}
        assert np.array_equal(model.ZGRID,cp['zgrid'])
        assert relative(model.ZDV,cp['zdv'])<1e-13
    with np.load(ROOT/'research/population-catalog/native-catalog-checkpoints.npz',allow_pickle=False) as cp:
        rmeans=cp['means']
    oldresult=json.loads((SEL/'selection-results.json').read_text())
    rearly=sum(np.array(oldresult['selection_components'][key]) for key in ['o1','o2'])
    policy=json.loads((SEL/'pipeline-results.json').read_text())
    rbaseline={'matched_filter':np.array(policy['selection_components']['matched_filter']),
               'with_cwb':np.array(policy['selection_components']['union'])}
    with h5py.File(path,'r') as h:
        group=h['injections']; fields=['mass1_source','mass2_source','redshift','sampling_pdf','mixture_weight']
        fields += [f'spin{i}{axis}' for i in [1,2] for axis in 'xyz']
        fields += [k for k in group if k.startswith('ifar_')]
        raw={k:group[k][:] for k in fields}
        attrs={k:float(h.attrs[k]) for k in ['analysis_time_s','total_generated','n_accepted','n_rejected']}
    oldsource=json.loads((SEL/'o3-provenance.json').read_text())
    assert sha(SEL/'o3-bank.npz')==oldsource['cache_sha256']
    assert attrs['analysis_time_s']==oldsource['attributes']['analysis_time_s']
    time_fraction=attrs['analysis_time_s']/oldsource['combined_attributes']['analysis_time_s']
    factor=time_fraction/attrs['total_generated']
    union=np.any([raw[k]>1 for k in raw if k.startswith('ifar_')],axis=0)
    matched=np.any([raw[k]>1 for k in raw if k.startswith('ifar_') and k!='ifar_cwb'],axis=0)
    found={k:v[union] for k,v in raw.items()}
    coordinates=physical(found)
    scales=found['mixture_weight'].astype(float)*factor
    m=found['mass1_source'].astype(float); n=found['mass2_source'].astype(float); z=found['redshift'].astype(float)
    old_domain=(m>=2)&(m<=100)&(n>=2)&(z<=1.9)&(coordinates['a_1']<=.998)&(coordinates['a_2']<=.998)
    masks={'matched_filter':matched[union],'with_cwb':np.ones(union.sum(),bool)}
    report=dict(schema='blackhole-official-o3-mixture-followup/1',event_count=69,
                timing='Design recorded before mixture file inspection/population integration; extension of previously observed results.',
                source_sha256=sha(path),o3_time_fraction=time_fraction,total_generated=attrs['total_generated'],
                fixed_early_term=True,normalization='sum found p_target/p_joint_mixture_draw * mixture_weight * (T_O3/T_all)/N_generated',
                metrics={},support=support_mass(),input_sha256={})
    saved={}
    for metric,means,early,baselines,weight_fn in [('g',gmeans,gearly,gbaseline,invariant_weights),('r',rmeans,rearly,rbaseline,native_weights)]:
        rowweights=weight_fn(coordinates)*scales[:,None,None]
        assert np.isfinite(rowweights).all() and np.all(rowweights>=0)
        metric_result=dict(functional='P(chi_eff<0)' if metric=='g' else 'P(cos_tilt1<0) at mixed native references',
                           scenarios={},selection_components={},new_support_contribution={},concentration={})
        saved[metric+'_means']=means; saved[metric+'_early']=early
        for policy_name,mask in masks.items():
            o3=rowweights[mask].sum(axis=0)
            restricted=rowweights[mask&old_domain].sum(axis=0)
            additional=rowweights[mask&~old_domain].sum(axis=0)
            assert relative(o3,restricted+additional)<1e-12
            selections={'pure_bbh_baseline':baselines[policy_name],
                        'mixture_bbh_domain_only':early+restricted,
                        'mixture_full_union':early+o3}
            policy_cases={}
            for case,selection in selections.items():
                result,grid,density=summarize(metric,means,selection)
                policy_cases[case]=result
                saved[metric+'_'+policy_name+'_'+case+'_selection']=selection
                saved[metric+'_'+policy_name+'_'+case+'_density']=density
            reference=np.array(policy_cases['pure_bbh_baseline']['quantiles'])
            baseline_density=saved[metric+'_'+policy_name+'_pure_bbh_baseline_density']
            for case,result in policy_cases.items():
                result['median_shift_from_pure_bbh_pp']=(100*(np.array(result['quantiles'])[:,1]-reference[:,1])).tolist()
                result['quantile_shift_from_pure_bbh_pp']=(100*(np.array(result['quantiles'])-reference)).tolist()
                result['posterior_total_variation_from_pure_bbh']=(.5*np.trapezoid(abs(saved[metric+'_'+policy_name+'_'+case+'_density']-baseline_density),grid,axis=1)).tolist()
            metric_result['scenarios'][policy_name]=policy_cases
            metric_result['selection_components'][policy_name]=dict(early=early.tolist(),o3_mixture=o3.tolist(),o3_old_domain=restricted.tolist(),o3_added_union=additional.tolist())
            metric_result['new_support_contribution'][policy_name]=dict(absolute=additional.tolist(),fraction_of_total_selection=(additional/(early+o3)).tolist())
            at_median=np.array(policy_cases['mixture_full_union']['quantiles'])[:,1]
            if metric=='r':
                from run_selection_study import coordinates as native_coordinates
                r,fs=native_coordinates(4001)
                mix=np.array([np.interp(at_median[i],r,fs[i]) for i in range(6)])
            else: mix=at_median
            combined=mix*rowweights[mask,:,0]+(1-mix)*rowweights[mask,:,1]
            metric_result['concentration'][policy_name]=dict(ess_at_median=(combined.sum(axis=0)**2/np.sum(combined**2,axis=0)).tolist(),
                largest_row_fraction_at_median=(combined.max(axis=0)/combined.sum(axis=0)).tolist(),
                scope='Weight concentration only, no PE/campaign independence or coverage claim.')
        report['metrics'][metric]=metric_result
        saved[metric+'_grid']=grid
    # Bridge the unchanged event terms to the already published aggregate fit.
    gold=json.loads((INV/'control-results.json').read_text())
    rold=json.loads((ROOT/'research/population-catalog/native-catalog-results.json').read_text())
    gerror=float(np.max(abs(np.array(report['metrics']['g']['scenarios']['matched_filter']['pure_bbh_baseline']['quantiles'])-gold['scenarios']['matched_filter']['quantiles'])))
    rerror=float(np.max(abs(np.array(report['metrics']['r']['scenarios']['matched_filter']['pure_bbh_baseline']['quantiles'])-rold['quantiles'])))
    assert max(gerror,rerror)<2e-6
    report['baseline_quantile_bridges']=dict(g=gerror,r=rerror)
    # Source-density diagnostic prompted by the completed comparison. Below
    # m1=90,m2>=3 only BBH has support. Do not assume the separately released
    # PDFs are exactly related by their total-trial mixture coefficient.
    with np.load(BASE/'mixture-membership.npz',allow_pickle=False) as mapping:
        join=mapping['pure_bbh_to_mixture']
    with np.load(SEL/'o3-bank.npz',allow_pickle=False) as bank:
        region=(bank['full_mass1_source']<90)&(bank['full_mass2_source']>=3)
        old_pdf=bank['full_sampling_pdf'][region]
    with h5py.File(path,'r') as h: new_pdf=h['injections/sampling_pdf'][:][join][region]
    ratio=new_pdf/old_pdf/(oldsource['attributes']['total_generated']/attrs['total_generated'])
    report['joint_proposal_bridge']=dict(timing='Post-comparison source diagnostic; no tolerance retrofitted as a reconstruction gate.',
        single_population_overlap_rows=int(region.sum()),joint_pdf_ratio_0_50_99_100=np.quantile(ratio,[0,.5,.99,1]).tolist(),
        maximum_relative_residual=float(np.max(abs(ratio-1))),
        interpretation='Physical/search rows are identical, but released mixture joint PDFs are not exactly a trial-fraction rescaling of the separately released BBH PDF. This study uses each official joint PDF as provided; it has not recovered the numerical density-construction history.')
    for p in [BASE/'design.md',BASE/'provenance.json',Path(__file__),ROOT/'model.py',INV/'run_control.py',INV/'invariant_model.py',INV/'control-checkpoints.npz',INV/'control-results.json',INV/'control-validation.json',ROOT/'research/population-catalog/native-catalog-checkpoints.npz',ROOT/'research/population-catalog/native-catalog-results.json',ROOT/'research/population-catalog/native-catalog-validation.json',SEL/'pipeline-results.json',SEL/'selection-results.json',SEL/'o3-provenance.json']:
        report['input_sha256'][p.relative_to(ROOT).as_posix()]=sha(p)
    report['limits']='Related fixed injection representations, unchanged all69 event likelihoods and restricted population models. Not calibrated actual-campaign uncertainty, complete target support, historical spin-reference recovery or formation-channel identification.'
    np.savez_compressed(BASE/'study-checkpoints.npz',**saved)
    report['checkpoint_sha256']=sha(BASE/'study-checkpoints.npz')
    (BASE/'results.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({metric:{policy:cases['mixture_full_union'] for policy,cases in value['scenarios'].items()} for metric,value in report['metrics'].items()},indent=2))
    print(json.dumps(report['support'],indent=2))


if __name__=='__main__': main()
