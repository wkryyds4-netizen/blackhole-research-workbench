"""Explain observed paired-policy dispersion; exploratory finite-bank diagnostic.

Rank rows by a delta-method influence on the *policy contrast*, not by their
measured thinning outcomes. Keep every row in the scientific estimator.
"""
from pathlib import Path
import json
import hashlib
import numpy as np
from scipy.integrate import cumulative_trapezoid
from run_selection_study import init_model,weights,physical,posterior,coordinates

BASE=Path(__file__).resolve().parent


def main():
    source=json.loads((BASE/'o3-provenance.json').read_text())
    policy=json.loads((BASE/'pipeline-results.json').read_text())
    assert hashlib.sha256((BASE/'o3-bank.npz').read_bytes()).hexdigest()==source['cache_sha256']
    with np.load(BASE/'o3-bank.npz',allow_pickle=False) as bank:
        raw={k[5:]:bank[k] for k in bank.files if k.startswith('full_')};actual=bank['combined_o3_to_full']
    with np.load(BASE/'pipeline-checkpoints.npz',allow_pickle=False) as bank:means=bank['means']
    with np.load(BASE/'paired-policy-checkpoints.npz',allow_pickle=False) as bank:
        xi=bank['selections'];gaps=bank['gaps'];cov=bank['analytic_covariance']
    full=np.array([policy['selection_components'][k] for k in ['union','matched_filter']])
    init_model();w=weights(physical(raw));m=len(w);n=len(actual)
    mf=np.maximum.reduce([raw[k] for k in ['ifar_gstlal','ifar_mbta','ifar_pycbc_bbh','ifar_pycbc_hyperbank']])>1
    union=mf|(raw['ifar_cwb']>1)
    masked=np.stack([w*union[:,None,None],w*mf[:,None,None]],axis=1)
    factor=source['attributes']['analysis_time_s']/source['combined_attributes']['analysis_time_s']/source['attributes']['total_generated']
    grads=[]
    for eps in [5e-7,1e-7]:
        gradient=np.zeros((2,6,2))
        for p in range(2):
            for k in [0,2,3,5]:
                for j in [0,1]:
                    plus=full[p].copy();minus=full[p].copy();h=eps*full[p,k,j]
                    plus[k,j]+=h;minus[k,j]-=h
                    a,_=posterior(means,plus);b,_=posterior(means,minus)
                    gradient[p,k,j]=(a[k,1]-b[k,1])/(2*h)
        grads.append(gradient)
    # Differentiate the same trapezoid CDF and linearly interpolated median exactly.
    r,fs=coordinates();grad=np.zeros((2,6,2))
    for p in range(2):
        _,density=posterior(means,full[p])
        for k in [0,2,3,5]:
            cdf=cumulative_trapezoid(density[k],r,initial=0);cdf/=cdf[-1]
            right=np.searchsorted(cdf,.5);left=right-1
            fraction=(.5-cdf[left])/(cdf[right]-cdf[left])
            efficiency=fs[k]*full[p,k,0]+(1-fs[k])*full[p,k,1]
            for j,basis in enumerate([fs[k],1-fs[k]]):
                score=-len(means)*basis/efficiency
                integral=cumulative_trapezoid(density[k]*score,r,initial=0)
                dc=integral-cdf*integral[-1]
                grad[p,k,j]=-(r[right]-r[left])/(cdf[right]-cdf[left])*((1-fraction)*dc[left]+fraction*dc[right])
    den=np.maximum(abs(grad),1e-30)
    derivative_difference=float(np.max(abs(grads[0]-grad)/den))
    derivative_fine=float(np.max(abs(grads[1]-grad)/den));assert derivative_fine<1e-5
    families=[];top_rows=[]
    for label,indices in [('low_spin',(0,2)),('baseline_spin',(3,5))]:
        a,b=indices;contrast=np.zeros_like(grad)
        contrast[1,a]=100*grad[1,a];contrast[1,b]=-100*grad[1,b]
        contrast[0,a]=-100*grad[0,a];contrast[0,b]=100*grad[0,b]
        row_score=np.einsum('npkj,pkj->n',masked,contrast);center=row_score-row_score.mean()
        order=np.argsort(abs(center))[::-1];top=int(order[0]);top_rows.append(top)
        shares=center*center/np.sum(center*center)
        linear=np.einsum('rpkj,pkj->r',xi-full,contrast)
        variance=float(contrast.ravel()@cov@contrast.ravel())
        direct_variance=float((factor*m)**2*(1-n/m)/n*np.sum(center*center)/(m-1))
        assert abs(variance/direct_variance-1)<1e-9
        fam=len(families);observed=gaps[:,1,fam]-gaps[:,0,fam]
        fullgap=policy['matched_filter_width_gaps_pp'][fam]-policy['union_width_gaps_pp'][fam]
        ranked=[]
        for row in order[:5]:
            ranked.append(dict(source_row=int(row),variance_share=float(shares[row]),cwb_only=bool(union[row] and not mf[row]),
                original_subset_contains=bool(row in actual),
                parameters={key:float(raw[key][row]) for key in ['mass1_source','mass2_source','redshift','spin1x','spin1y','spin1z','spin2x','spin2y','spin2z','optimal_snr_net']},
                ifar={key:float(raw[key][row]) for key in raw if key.startswith('ifar_')},
                predicted_conditional_inclusion_difference_pp=float(factor*m/n*m/(m-1)*center[row])))
        families.append(dict(label=label,ranked_rows=ranked,top5_linear_variance_share=float(shares[order[:5]].sum()),
            linear_predicted_std_pp=float(np.sqrt(variance)),empirical_policy_contrast_std_pp=float(np.std(observed,ddof=1)),
            linear_observed_correlation=float(np.corrcoef(linear,observed-fullgap)[0,1]),
            linear_prediction_rmse_pp=float(np.sqrt(np.mean((linear-(observed-fullgap))**2))),
            covariance_projection_relative_error=float(abs(variance/direct_variance-1))))
    # Reproduce exactly the previous selection indices, now record top-row membership.
    rng=np.random.default_rng(20261025);presence=[]
    for _ in range(512):
        take=rng.choice(m,n,replace=False);presence.append([bool(row in take) for row in top_rows])
    presence=np.array(presence)
    for k,fam in enumerate(families):
        delta=gaps[:,1,k]-gaps[:,0,k];yes=presence[:,k]
        fam['top_row_inclusion_groups']=dict(included_replicates=int(yes.sum()),excluded_replicates=int((~yes).sum()),
            included_contrast_5_50_95_pp=np.quantile(delta[yes],[.05,.5,.95]).tolist(),excluded_contrast_5_50_95_pp=np.quantile(delta[~yes],[.05,.5,.95]).tolist(),
            included_minus_excluded_mean_pp=float(delta[yes].mean()-delta[~yes].mean()))
    result=dict(schema='blackhole-policy-influence/1',families=families,
        derivative_coarse_relative_error=derivative_difference,derivative_fine_relative_error=derivative_fine,
        source_sha256=source['cache_sha256'],
        interpretation='Exploratory delta-method diagnosis on a finite fixed bank; row variance shares are linearized, not exact causal attribution of nonlinear medians. No rows deleted, reweighted away, or labeled physically invalid.')
    (BASE/'policy-influence-results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    np.savez_compressed(BASE/'policy-influence-checkpoints.npz',top_rows=np.array(top_rows),presence=presence)
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
