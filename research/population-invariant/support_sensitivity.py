"""Conditional likelihood-shape envelopes and explicit missing-component stress.

No sampled row is removed, and no missing likelihood is inferred from a zero
bank count. The event-completion family is an integral-level counterfactual.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.integrate import cumulative_trapezoid, quad
from scipy.optimize import brentq
from scipy.special import betaincc, expit, logit
from astropy.cosmology import Planck15

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
SEL=ROOT/'research/population-selection-audit'
LAMBDAS=[0,4,16]
MAGNITUDES=[(1.5,8.5),(2.,5.)]


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def summarize(grid,density):
    p=density/np.trapezoid(density,grid,axis=-1)[...,None]
    c=cumulative_trapezoid(p,grid,axis=-1,initial=0);c/=c[...,-1,None]
    q=np.array([np.interp([.05,.5,.95],cdf,grid) for cdf in c.reshape(-1,len(grid))])
    return p,c,q.reshape(*p.shape[:-1],3)


def main():
    report=json.loads((BASE/'control-results.json').read_text())
    validation=json.loads((BASE/'control-validation.json').read_text())
    assert validation['result_sha256']==sha(BASE/'control-results.json')
    with np.load(BASE/'control-checkpoints.npz',allow_pickle=False) as cp:
        g=cp['g'];means=cp['means'];saved=cp['densities'][:2]
        selection=np.stack([cp['selection_matched'],cp['selection_union']])
        zgrid=cp['zgrid'];zdv=cp['zdv'];znorm=float(cp['znorm'])
    integrals=json.loads((BASE/'event-integrals.json').read_text())['events']
    frozen=json.loads((ROOT/'results/model.json').read_text())
    events=frozen['events'];pos=next(i for i,e in enumerate(events) if e['index']==54)
    event=next(e for e in integrals if e['index']==54)
    assert np.all(means[pos,:,0]==0) and np.all(means[pos,:,1]>0)
    p,cdf,q=summarize(g,saved)
    baseline=np.array([report['scenarios'][name]['quantiles'] for name in ['matched_filter','with_cwb']])
    assert np.max(abs(q-baseline))<1e-12
    budgets=[0.,.05,.1,.2,.5,1.,2.]
    envelopes=[];thresholds=[]
    for policy in range(2):
        for k,lam in enumerate(LAMBDAS):
            lower=cdf[policy,k];upper=cdf[policy,k+3]
            median_lower=q[policy,k,1];median_upper=q[policy,k+3,1]
            def gap(B):
                return np.interp(expit(-B),upper,g)-np.interp(expit(B),lower,g)
            critical=brentq(gap,0,30,xtol=1e-12)
            thresholds.append(dict(policy=policy,lambda_shape=lam,full_median_difference_pp=float(100*(median_upper-median_lower)),
                critical_log_correction_oscillation=float(critical),
                equal_per_event_relative_omission_bound_at_threshold=float(np.expm1(critical/69)),
                uniform_relative_selection_completion_bound_at_threshold=float(np.expm1(critical/69)),
                interpretation='Sufficient outer-envelope ordering certificate is positive below this B. Failure of the certificate is not an observed physical reversal.'))
        for B in budgets:
            lo=np.array([np.interp(expit(-B),row,g) for row in cdf[policy]])
            hi=np.array([np.interp(expit(B),row,g) for row in cdf[policy]])
            envelopes.append(dict(policy=policy,B=B,median_lower=lo.tolist(),median_upper=hi.tolist(),
                amplitude_median_gap_outer_lower_pp=(100*(lo[3:]-hi[:3])).tolist(),
                amplitude_median_gap_outer_upper_pp=(100*(hi[3:]-lo[:3])).tolist(),
                total_variation_upper_bound=float(np.tanh(B/4))))
    # Evaluate new event factors directly; division by (1-g) would lose the
    # finite eta>0 endpoint contribution at g=1.
    completions=[];curves=[]
    eta_values=[0.,.0001,.001,.01,.1,1.,10.,100.,None]
    for eta in eta_values:
        mu=means.copy()
        if eta is None:mu[pos,:,0]=mu[pos,:,1];mu[pos,:,1]=0
        else:mu[pos,:,0]=eta*mu[pos,:,1]
        mix=g*mu[:,:,0,None]+(1-g)*mu[:,:,1,None]
        with np.errstate(divide='ignore'):lognum=np.log(mix).sum(axis=0)
        densities=[];quantiles=[]
        for policy in range(2):
            eff=g*selection[policy,:,0,None]+(1-g)*selection[policy,:,1,None]
            lp=lognum-69*np.log(eff)
            density=np.exp(lp-lp.max(axis=1)[:,None]);density,_,quant=summarize(g,density)
            densities.append(density);quantiles.append(quant)
        quantiles=np.array(quantiles);curves.append(densities)
        completions.append(dict(eta=eta,eta_label='infinity_limit' if eta is None else str(eta),
            quantiles=quantiles.tolist(),median_shift_pp=(100*(quantiles[...,1]-q[...,1])).tolist(),
            amplitude_median_difference_pp=(100*(quantiles[:,3:,1]-quantiles[:,:3,1])).tolist()))
    assert np.max(abs(np.array(completions[0]['quantiles'])-baseline))<2e-12
    source=json.loads((SEL/'o3-provenance.json').read_text())
    earlier=json.loads((SEL/'selection-results.json').read_text())
    attrs=source['attributes'];fraction=earlier['o3_time_fraction']
    fz=lambda z:Planck15.differential_comoving_volume(z).value*(1+z)
    continuous_tail=quad(fz,attrs['max_z'],2.3,epsabs=.0001,epsrel=1e-11)[0]/quad(fz,0,2.3,epsabs=.0001,epsrel=1e-11)[0]
    # Match the *implemented* linearly interpolated dV/dz and its stored
    # trapezoid normalizer. Integrating at interval midpoints is exact for
    # each quadratic (1+z)*interp(dV/dz). The tiny normalization offset
    # cancels in inference, but must also scale a selection-contribution bound.
    def coded_mass(lo,hi):
        nodes=np.r_[lo,zgrid[(zgrid>lo)&(zgrid<hi)],hi]
        mid=(nodes[:-1]+nodes[1:])/2
        f=lambda z:np.interp(z,zgrid,zdv)*(1+z)
        return float(np.sum(np.diff(nodes)/6*(f(nodes[:-1])+4*f(mid)+f(nodes[1:])))/znorm)
    total_mass=coded_mass(0,2.3);coded_tail=coded_mass(attrs['max_z'],2.3)
    ztail=coded_tail/total_mass
    domain=[]
    for j,(a,b) in enumerate(MAGNITUDES):
        s1=float(betaincc(a,b,attrs['max_spin1']));s2=float(betaincc(a,b,attrs['max_spin2']))
        spin_tail=s1+s2-s1*s2
        mass=ztail+(1-ztail)*spin_tail;coded_domain=coded_tail+(total_mass-coded_tail)*spin_tail
        upper=fraction*coded_domain
        pe99=float(betaincc(a,b,.99));pe_union=2*pe99-pe99**2
        entries=[]
        for policy in range(2):
            for k in range(3*j,3*j+3):
                eff=g*selection[policy,k,0]+(1-g)*selection[policy,k,1]
                atmedian=float(np.interp(q[policy,k,1],g,eff))
                entries.append(dict(policy=policy,model=k,represented_efficiency_at_median=atmedian,
                    omitted_domain_contribution_upper_over_represented_at_median=float(upper/atmedian),
                    conservative_log_shape_bound_from_this_domain_only=float(69*np.log1p(upper/eff.min()))))
        domain.append(dict(beta=[a,b],population_redshift_tail_mass=float(ztail),
            population_spin_union_above_injection_caps=float(spin_tail),known_omitted_union_population_mass=float(mass),
            known_omitted_union_mass_in_coded_selection_measure=float(coded_domain),
            o3_time_weighted_selection_contribution_upper=float(upper),
            population_spin_union_above_point99=float(pe_union),model_details=entries))
    result=dict(schema='blackhole-support-sensitivity/1',events=69,keys=report['keys'],budgets=envelopes,
        producer_sha256=sha(Path(__file__)),
        ordering_certificate_thresholds=thresholds,completion_event=dict(index=54,event=event['event'],rows=event['rows'],
            zero_negative_components=6,scope='An observed finite-bank zero, not a likelihood zero. eta is missing negative integral divided by the retained positive integral, not missing posterior probability.'),
        event_completions=completions,known_o3_omitted_domain=domain,o3_time_fraction=fraction,
        redshift_quadrature_audit=dict(continuous_normalized_tail_mass=float(continuous_tail),tabulated_normalized_tail_mass=float(ztail),
            implemented_density_total_mass=total_mass,implemented_density_tail_mass=coded_tail,
            bound_uses='Coded selection measure, including its stored-normalizer offset; normalized population mass is reported separately.'),
        formulas=dict(cdf_lower='F/(F+exp(B)*(1-F))',cdf_upper='exp(B)*F/(1-F+exp(B)*F)',
            quantile_lower='Q(expit(logit(q)-B))',quantile_upper='Q(expit(logit(q)+B))',
            log_shape_budget='sum_i log(1+epsilon_i) + 69 log(1+eta_selection)',
            assumption='Uniform relative bounds on missing event integrals and missing selection contributions, conditional on exact represented integrals. Undefined where a retained integral is zero and a completion is positive.'),
        input_sha256={name:sha(BASE/name) for name in ['support-design.json','control-results.json','control-validation.json','control-checkpoints.npz','event-integrals.json']},
        source_sha256={name:sha(SEL/name) for name in ['o3-provenance.json','selection-results.json']},
        limitations='B, epsilon and eta are assumptions, not estimated error bars. Adversarial envelopes optimize over a larger class than a single common physical likelihood/detection law. O3-domain bound covers only the stated omitted domain; other support and sampling errors remain. No channel fractions or unconditional robustness certification.')
    (BASE/'support-results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    np.savez_compressed(BASE/'support-checkpoints.npz',g=g,baseline=p,cdf=cdf,completion_densities=np.array(curves))
    print(json.dumps(dict(thresholds=thresholds,event_completion=event['event'],
        largest_finite_completion_shift_pp=float(np.max(abs(np.array([r['median_shift_pp'] for r in completions[:-1]])))),
        largest_infinite_limit_shift_pp=float(np.max(abs(np.array(completions[-1]['median_shift_pp'])))),domain=domain),indent=2))


if __name__=='__main__':main()
