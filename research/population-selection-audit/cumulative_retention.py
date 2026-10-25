"""Explicitly hypothetical observed-SNR retention model for the early bank.

Assumes independent Gaussian quadratures, 2*N_IFO degrees of freedom, and
retention exactly at observed network SNR>8.5 (plus parameter-only cuts).
This describes the saved data well but is NOT recovered generation metadata.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import json
import hashlib
import numpy as np
from scipy.stats import ncx2,kstest
from run_selection_study import init_model,physical,weights,posterior
from cumulative_control import LOGPDF

BASE=Path(__file__).resolve().parent


def main():
    init_model();prov=json.loads((BASE/'cumulative-provenance.json').read_text())
    assert hashlib.sha256((BASE/'cumulative-bank.npz').read_bytes()).hexdigest()==prov['cache_sha256']
    with np.load(BASE/'cumulative-bank.npz',allow_pickle=False) as b:d={k:b[k] for k in b.files}
    with np.load(BASE/'cumulative-checkpoints.npz',allow_pickle=False) as cp:means=cp['means'];r=cp['r']
    previous=json.loads((BASE/'cumulative-results.json').read_text())
    early=d['early'];raw={k:v[early] for k,v in d.items() if v.shape==early.shape}
    raw['sampling_pdf']=np.exp(raw[LOGPDF]);w=weights(physical(raw))*raw['weights'][:,None,None]/prov['attributes']['total_generated']
    opt=raw['estimated_optimal_snr_net'];obs=raw['semianalytic_observed_phase_maximized_snr_net']
    dof=2*sum((raw['estimated_optimal_snr_'+ifo]>0).astype(int) for ifo in 'HLV')
    assert np.all(obs>8.5) and set(dof)=={4,6}
    quadrature=np.sqrt(sum(raw['estimated_optimal_snr_'+ifo]**2 for ifo in 'HLV'))
    assert np.max(abs(quadrature/opt-1))<1e-12
    survive=ncx2.sf(8.5**2,dof,opt**2)
    prob=ncx2.sf(10**2,dof,opt**2)/survive
    assert np.all((prob>=0)&(prob<=1))
    raw_pit=ncx2.cdf(obs**2,dof,opt**2)
    retained_pit=1-ncx2.sf(obs**2,dof,opt**2)/survive
    ks_all=kstest(raw_pit,'uniform');ks_ret=kstest(retained_pit,'uniform')
    diagnostics=[]
    for df in [4,6]:
        mask=dof==df;v=prob[mask];actual=obs[mask]>10
        ks=kstest(retained_pit[mask],'uniform')
        diagnostics.append(dict(dof=df,rows=int(mask.sum()),observed_minimum=float(obs[mask].min()),conditional_pit_ks_distance=float(ks.statistic),
            threshold_count_standardized_residual=float((actual.sum()-v.sum())/np.sqrt(np.sum(v*(1-v))))))
    # Inverse-retention identity: E[I(R>8.5) f(theta)/P(R>8.5|theta)] = E[f(theta)].
    # Use it only on opt>10; there, inclusion probability is high and weights stable.
    keep=opt>10
    o3=np.array(previous['partitions']['o3_matched'])
    corrected=o3+np.sum(w[keep]/survive[keep,None,None],axis=0)
    expected=o3+np.sum(w*prob[:,None,None],axis=0)
    selections={'retained_optimal_mask':np.array(previous['selection_components']['cumulative_optimal']),
        'hypothetical_retention_corrected_optimal':corrected,
        'released_observed_mask':np.array(previous['selection_components']['cumulative_observed']),
        'conditional_expected_observed':expected}
    cases={};curves=[];qs=[];maxgrid=0.
    for name,s in selections.items():
        q,p=posterior(means,s);fine,_=posterior(means,s,2001);maxgrid=max(maxgrid,float(np.max(abs(q-fine))))
        cases[name]=dict(quantiles=q.tolist(),width_gaps_pp=[float(100*(q[0,1]-q[2,1])),float(100*(q[3,1]-q[5,1]))]);qs.append(q);curves.append(p)
    assert maxgrid<1e-6
    # Resample only the remaining noise conditional on these saved parameters
    # and the assumed retention. Not a redraw of the injection campaign.
    rng=np.random.default_rng(20261026);repeats=256
    uncertain=(prob>0)&(prob<1);certain=prob==1
    flat=w[uncertain].reshape((-1,12));p=prob[uncertain]
    base=o3+w[certain].sum(axis=0)
    simulations=[]
    for start in range(0,repeats,32):
        masks=(rng.random((32,len(p)))<p).astype(float)
        simulations.extend(base[None,:,:]+(masks@flat).reshape((-1,6,2)))
    simulations=np.array(simulations)
    sq=np.array([posterior(means,s)[0] for s in simulations])
    shift=100*(sq[:,:,1]-qs[3][:,1])
    covariance=(flat.T*(p*(1-p)))@flat
    empirical=np.cov(simulations.reshape((repeats,12)),rowvar=False,ddof=1)
    report=dict(schema='blackhole-hypothetical-early-retention/1',generation_configuration_recovered=False,
        hypothesis='Observed-SNR>8.5 retention with noncentral chi-square R^2 ~ ncx2(df=2*N_IFO,nc=optimal_snr^2), plus parameter-only cuts; independent noise quadratures. Empirical model check, not historical configuration proof.',
        source_sha256=prov['source_sha256'],minimum_saved_observed_snr=float(obs.min()),
        diagnostics=dict(unconditional_pit_ks_distance=float(ks_all.statistic),conditional_pit_ks_distance=float(ks_ret.statistic),
            by_network=diagnostics,optimal_selected_retention_probability_range=[float(survive[keep].min()),float(survive[keep].max())],
            all_retention_probability_range=[float(survive.min()),float(survive.max())],
            note='The 8.5 cutoff hypothesis was suggested by the observed minimum. KS distances are descriptive; do not treat this as a prespecified significance test.'),
        scenarios=cases,corrected_optimal_minus_retained_optimal_median_pp=(100*(qs[1][:,1]-qs[0][:,1])).tolist(),
        expected_observed_minus_corrected_optimal_median_pp=(100*(qs[3][:,1]-qs[1][:,1])).tolist(),
        released_observed_minus_expected_median_pp=(100*(qs[2][:,1]-qs[3][:,1])).tolist(),
        conditional_noise_simulation=dict(seed=20261026,repeats=repeats,median_shift_5_50_95_pp=np.quantile(shift,[.05,.5,.95],axis=0).T.tolist(),
            median_shift_sd_pp=shift.std(axis=0,ddof=1).tolist(),empirical_to_analytic_selection_variance=(np.diag(empirical)/np.diag(covariance)).tolist(),
            scope='Fixed retained source parameters and fixed O3. Only truncated early Gaussian noise is resimulated. Not full selection MC or astrophysical coverage.'),
        max_grid_quantile_error=maxgrid,
        limits='Correction is conditional on an inferred, unconfirmed retention law and Gaussian model. Unobserved source regions, optimal-SNR cuts, campaign time mixing, spin reference and PE supports are not repaired. Do not substitute this for the published baseline.')
    (BASE/'cumulative-retention-results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    np.savez_compressed(BASE/'cumulative-retention-checkpoints.npz',r=r,names=np.array(list(selections)),selections=np.array(list(selections.values())),densities=np.array(curves),quantiles=np.array(qs),simulated_selections=simulations,simulated_quantiles=sq,analytic_covariance=covariance)
    print(json.dumps({k:report[k] for k in ['diagnostics','corrected_optimal_minus_retained_optimal_median_pp','expected_observed_minus_corrected_optimal_median_pp','released_observed_minus_expected_median_pp','conditional_noise_simulation','max_grid_quantile_error']},indent=2))


if __name__=='__main__':main()
