"""Independent product, Simpson-CDF and direct-root support diagnostics."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.integrate import cumulative_simpson, quad
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.special import beta

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    path=BASE/'support-results.json';result=json.loads(path.read_text())
    assert result['producer_sha256']==sha(BASE/'support_sensitivity.py')
    for name,digest in result['input_sha256'].items():assert sha(BASE/name)==digest
    sel=ROOT/'research/population-selection-audit'
    for name,digest in result['source_sha256'].items():assert sha(sel/name)==digest
    with np.load(BASE/'control-checkpoints.npz',allow_pickle=False) as cp:
        g=cp['g'];means=cp['means'];selections=[cp['selection_matched'],cp['selection_union']]
        zgrid=cp['zgrid'];zdv=cp['zdv'];znorm=float(cp['znorm'])
    with np.load(BASE/'support-checkpoints.npz',allow_pickle=False) as cp:
        saved=cp['completion_densities'];nominal=cp['baseline']
    events=json.loads((ROOT/'results/model.json').read_text())['events']
    pos=next(i for i,e in enumerate(events) if e['index']==54)
    assert np.flatnonzero(np.any(means[:,:,0]==0,axis=1)).tolist()==[pos]
    maximum_density=0.;maximum_quantile=0.;cdf_functions={};quantile_functions={}
    for completion_index,completion in enumerate(result['event_completions']):
        eta=completion['eta'];other=np.delete(means,pos,axis=0)
        factors=g*other[:,:,0,None]+(1-g)*other[:,:,1,None]
        mid=len(g)//2
        num=np.prod((factors/factors[:,:,mid,None]).astype(np.longdouble),axis=0)
        replacement=g if eta is None else ((1-g)+eta*g)
        num*=replacement[None,:]/replacement[mid]
        for policy,selection in enumerate(selections):
            eff=g*selection[:,0,None]+(1-g)*selection[:,1,None]
            density=num*(eff[:,mid,None]/eff).astype(np.longdouble)**69
            density=np.asarray(density/density.max(axis=1)[:,None],dtype=float)
            density/=np.trapezoid(density,g,axis=1)[:,None]
            maximum_density=max(maximum_density,float(np.max(abs(density-saved[completion_index,policy]))))
            for k in range(6):
                c=cumulative_simpson(density[k],x=g,initial=0);c/=c[-1]
                cfun=PchipInterpolator(g,c);qfun=lambda prob,cfun=cfun:brentq(lambda x:float(cfun(x))-prob,0,1,xtol=1e-13)
                quant=np.array([qfun(prob) for prob in [.05,.5,.95]])
                maximum_quantile=max(maximum_quantile,float(np.max(abs(quant-completion['quantiles'][policy][k]))))
                if completion_index==0:
                    cdf_functions[policy,k]=cfun
                    quantile_functions[policy,k]=qfun
    assert maximum_density<1e-10
    assert maximum_quantile<2e-6
    # Opposing median-envelope endpoints touch at a shared x whose two CDFs
    # sum to one. This solves for x first, independently of the producer's
    # root in B with interpolated inverse CDFs.
    threshold_errors=[]
    for row in result['ordering_certificate_thresholds']:
        p=row['policy'];k=[0,4,16].index(row['lambda_shape'])
        F=cdf_functions[p,k];G=cdf_functions[p,k+3]
        x=brentq(lambda x:float(F(x)+G(x))-1,0,1,xtol=1e-13)
        prob=float(F(x));B=np.log(prob/(1-prob))
        threshold_errors.append(abs(B-row['critical_log_correction_oscillation']))
    assert max(threshold_errors)<2e-5
    # Explicit split-domain reweighting attains both CDF bounds on each
    # nominal density. Insert split points before integration; never apply a
    # quadrature rule across a discontinuous weight.
    bound_errors=[];tv_errors=[]
    for B in [0.,.05,.1,.2,.5,1.,2.]:
        R=np.exp(B)
        for policy in range(2):
            for k in range(6):
                density=nominal[policy,k]
                for probability in [.1,.5,.9]:
                    x=quantile_functions[policy,k](probability)
                    lo=np.r_[g[g<x],x];hi=np.r_[x,g[g>x]]
                    a=np.trapezoid(np.interp(lo,g,density),lo)
                    b=np.trapezoid(np.interp(hi,g,density),hi)
                    lower=a/(a+R*b);upper=R*a/(R*a+b)
                    bound_errors.extend([abs(lower-probability/(probability+R*(1-probability))),
                        abs(upper-R*probability/(1-probability+R*probability))])
        t=1/(1+np.sqrt(R));new=R*t/(R*t+1-t)
        tv_errors.append(abs((new-t)-np.tanh(B/4)))
    assert max(bound_errors)<2e-6
    assert max(tv_errors)<1e-14
    # Independent cosmology: integrate the FLRW distance equation directly
    # instead of importing Astropy or the model's tabulated volume factors.
    DH=299792.458/67.74
    # Planck15 has massive neutrinos; reproduce its expansion function from
    # its public immutable parameters rather than dropping them silently.
    # Astropy is used only for E(z); the two distance/volume integrals and
    # normalization are independently evaluated here.
    from astropy.cosmology import Planck15
    def volume(z):
        dc=DH*quad(lambda t:1/Planck15.efunc(t),0,z,epsabs=1e-12)[0]
        return DH*dc**2/Planck15.efunc(z)*(1+z)
    tail=quad(volume,1.9,2.3,epsabs=.0001,epsrel=1e-10)[0]/quad(volume,0,2.3,epsabs=.0001,epsrel=1e-10)[0]
    audit=result['redshift_quadrature_audit']
    volume_error=abs(tail-audit['continuous_normalized_tail_mass'])
    def polynomial_mass(lo,hi):
        edges=np.r_[lo,zgrid[(zgrid>lo)&(zgrid<hi)],hi]
        left=edges[:-1];h=np.diff(edges)
        values=np.interp(edges,zgrid,zdv);slope=np.diff(values)/h
        return float(np.sum(values[:-1]*(1+left)*h+(values[:-1]+slope*(1+left))*h*h/2+slope*h**3/3)/znorm)
    total=polynomial_mass(0,2.3);coded_tail=polynomial_mass(1.9,2.3)
    mass_errors=[abs(total-audit['implemented_density_total_mass']),abs(coded_tail-audit['implemented_density_tail_mass'])]
    tail_errors=[];spin_errors=[]
    for row in result['known_o3_omitted_domain']:
        a,b=row['beta'];tail_errors.append(abs(coded_tail/total-row['population_redshift_tail_mass']))
        for cap,target in [(.998,row['population_spin_union_above_injection_caps']),(.99,row['population_spin_union_above_point99'])]:
            one=quad(lambda x:x**(a-1)*(1-x)**(b-1)/beta(a,b),cap,1,epsabs=1e-40,epsrel=1e-10)[0]
            spin_errors.append(abs((2*one-one*one)/target-1))
    assert volume_error<1e-9 and max(tail_errors+mass_errors)<1e-13 and max(spin_errors)<1e-8
    validation=dict(schema='blackhole-support-independent-validation/1',
        independent_completions_checked=9,policies=2,models=6,maximum_density_absolute_error=maximum_density,
        maximum_simpson_root_quantile_error=maximum_quantile,quantile_gate=2e-6,
        maximum_direct_shared_location_B_error=float(max(threshold_errors)),threshold_gate=2e-5,
        maximum_split_quadrature_cdf_envelope_error=float(max(bound_errors)),cdf_gate=2e-6,
        maximum_two_point_extremal_tv_error=float(max(tv_errors)),
        maximum_independent_volume_tail_error=float(volume_error),
        maximum_tabulated_mass_antiderivative_error=float(max(tail_errors+mass_errors)),
        maximum_spin_survival_integral_relative_error=float(max(spin_errors)),
        result_sha256=sha(path),validator_sha256=sha(Path(__file__)),
        scope='Numerical/algebraic validation only. Input event means were previously source-validated. These checks do not establish any physical missing-likelihood budget or realizability of all integral-level completions by one waveform likelihood.')
    (BASE/'support-validation.json').write_text(json.dumps(validation,indent=2),encoding='utf-8')
    print(json.dumps(validation,indent=2))


if __name__=='__main__':main()
