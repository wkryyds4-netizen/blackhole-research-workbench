"""Retained-row O1/O2 SNR control on one fixed cumulative O1/O2/O3 bank.

The historical-to-cumulative comparison is a bundled representation change;
observed-versus-optimal SNR here changes only masks on retained rows. The early
observed SNR has an empirical floor near 8.5, so the optimal mask cannot be
treated as a complete deterministic-threshold counterfactual.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from run_selection_study import init_model,physical,weights,posterior,KEYS
import model

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
LOGPDF='lnpdraw_mass1_source_mass2_source_redshift_spin1x_spin1y_spin1z_spin2x_spin2y_spin2z'


def main():
    init_model()
    provenance=json.loads((BASE/'cumulative-provenance.json').read_text())
    assert hashlib.sha256((BASE/'cumulative-bank.npz').read_bytes()).hexdigest()==provenance['cache_sha256']
    with np.load(BASE/'cumulative-bank.npz',allow_pickle=False) as bank:d={k:bank[k] for k in bank.files}
    early=d['early'];o3=~early
    # Preserve the release's mixture weights and global normalization.
    d['sampling_pdf']=np.exp(d[LOGPDF]);coords=physical(d)
    w=weights(coords)*d['weights'][:,None,None]/provenance['attributes']['total_generated']
    assert np.isfinite(w).all() and (w>=0).all()
    support=(coords['mass_1']>=3)&(coords['mass_1']<=100)&(d['mass2_source']>=3)
    assert np.all(w[~support]==0)
    assert np.all((coords['a_1']>0)&(coords['a_1']<=1)&(coords['a_2']>0)&(coords['a_2']<=1))
    matched=np.any([d[k]<1 for k in d if k.startswith('o3_') and k.endswith('_far') and k!='o3_cwb_far'],axis=0)
    assert not matched[early].any()
    observed=early&(d['semianalytic_observed_phase_maximized_snr_net']>10)
    optimal=early&(d['estimated_optimal_snr_net']>10)
    gained=observed&~optimal;lost=optimal&~observed
    with np.load(BASE/'pipeline-checkpoints.npz',allow_pickle=False) as cp:
        means=cp['means'];r=cp['r'];fs=cp['fs'];historical=cp['selection_matched_filter']
    partitions={name:w[sel].sum(axis=0) for name,sel in [('o3_matched',matched),('early_optimal',optimal),('early_observed',observed),('gained',gained),('lost',lost)]}
    selections={'historical_full_o3_matched':historical,
        'cumulative_optimal':partitions['o3_matched']+partitions['early_optimal'],
        'cumulative_observed':partitions['o3_matched']+partitions['early_observed']}
    assert np.max(abs(partitions['early_observed']-partitions['early_optimal']-partitions['gained']+partitions['lost']))<1e-15
    scenarios={};curves=[];quantiles=[];grid_errors=[]
    for name,selection in selections.items():
        q,p=posterior(means,selection);fine,_=posterior(means,selection,2001)
        grid_errors.append(float(np.max(abs(q-fine))));assert grid_errors[-1]<1e-6
        scenarios[name]=dict(quantiles=q.tolist(),width_gaps_pp=[float(100*(q[0,1]-q[2,1])),float(100*(q[3,1]-q[5,1]))])
        curves.append(p);quantiles.append(q)
    delta=100*(quantiles[2][:,1]-quantiles[1][:,1]);diagnostics=[]
    for k,key in enumerate(KEYS):
        f=np.interp(quantiles[2][k,1],r,fs[k]);v=f*w[:,k,0]+(1-f)*w[:,k,1]
        diag=dict(key=key)
        for label,sel in [('observed_total',observed|matched),('optimal_total',optimal|matched),('observed_early',observed),('o3_matched',matched),('gained',gained),('lost',lost)]:
            one=v[sel];diag[label]=dict(sum=float(one.sum()),ess=float(one.sum()**2/(one@one)),max_fraction=float(one.max()/one.sum()))
        diag['early_fraction']=float(v[observed].sum()/v[observed|matched].sum())
        diag['selection_ratio_observed_over_optimal']=float(v[observed|matched].sum()/v[optimal|matched].sum())
        xi_o=fs[k]*selections['cumulative_observed'][k,0]+(1-fs[k])*selections['cumulative_observed'][k,1]
        xi_p=fs[k]*selections['cumulative_optimal'][k,0]+(1-fs[k])*selections['cumulative_optimal'][k,1]
        ratio=xi_o/xi_p
        diag['selection_ratio_over_r_range']=[float(ratio.min()),float(ratio.max())]
        diag['paired_log_likelihood_change_range']=float(len(means)*np.ptp(np.log(ratio)))
        diagnostics.append(diag)
    # A common overall time/normalization factor must cancel in conditional shape.
    scaled,_=posterior(means,selections['cumulative_observed']*provenance['attributes']['total_analysis_time'])
    scale_error=float(np.max(abs(scaled-quantiles[2])));assert scale_error<1e-12
    result=dict(schema='blackhole-cumulative-paired-snr-control/1',keys=KEYS,scenarios=scenarios,
        counts={name:int(sel.sum()) for name,sel in [('all',np.ones(len(early),dtype=bool)),('target_mass_support',support),('early',early),('o3',o3),('early_observed',observed),('early_optimal',optimal),('observed_only',gained),('optimal_only',lost),('o3_matched',matched)]},
        paired_observed_minus_optimal_median_pp=delta.tolist(),
        paired_observed_minus_optimal_gap_pp=(np.array(scenarios['cumulative_observed']['width_gaps_pp'])-scenarios['cumulative_optimal']['width_gaps_pp']).tolist(),
        selection_components={k:v.tolist() for k,v in selections.items()},partitions={k:v.tolist() for k,v in partitions.items()},
        diagnostics=diagnostics,validation=dict(max_grid_quantile_error=max(grid_errors),global_normalization_quantile_error=scale_error),
        source_sha256=provenance['source_sha256'],cache_sha256=provenance['cache_sha256'],
        fixed='69 event means with five original-prior replacements; six models; common uniform r prior; same O3 matched-filter FAR<1/year rule.',
        optimal_scenario_scope='The cumulative_optimal key denotes optimal SNR>10 ONLY among retained rows. All early rows have observed SNR>8.500004; unrecorded rows may satisfy optimal>10. This is not the full optimal-threshold counterfactual.',
        early_observed_minimum=float(d['semianalytic_observed_phase_maximized_snr_net'][early].min()),
        limits='Retained-row mask control only. An empirical early observed-SNR floor near 8.5 invalidates interpreting cumulative_optimal as a complete deterministic-threshold estimate without a retention model. New-vs-historical difference combines bank, mixture, coverage and time changes. Native spin reference, proposal tails and PE support remain unresolved. No O4 data used.')
    (BASE/'cumulative-results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    np.savez_compressed(BASE/'cumulative-checkpoints.npz',r=r,fs=fs,means=means,names=np.array(list(selections)),
        selections=np.array(list(selections.values())),densities=np.array(curves),quantiles=np.array(quantiles),
        zgrid=model.ZGRID,zdv=model.ZDV,znorm=model.ZNORM[2])
    print(json.dumps({k:result[k] for k in ['counts','scenarios','paired_observed_minus_optimal_median_pp','paired_observed_minus_optimal_gap_pp','diagnostics','validation']},indent=2))


if __name__=='__main__':main()
