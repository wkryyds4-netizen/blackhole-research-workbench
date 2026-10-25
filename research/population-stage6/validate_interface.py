"""Synthetic checks and stage-five regression for weighted event representations."""
from pathlib import Path
import sys
import hashlib
import json
import tempfile
import numpy as np

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
sys.path.insert(0, str(ROOT / 'research'))
sys.path.insert(0, str(ROOT / 'research/population-stage5'))
from population_event import WeightedEvent
from run_tail_study import density_pair, summarize, init_model, weights, PARAMS


def metadata(source, event='synthetic', coordinates=None):
    return dict(event=event, analysis_label='C01:IMRPhenomXPHM' if coordinates else 'analytic fixture',
                source=source, row_basis='zero-based original file row',
                measure=dict(coordinates=coordinates or ['x'],
                             units={k: ('solar_mass' if k == 'mass_1' else 'dimensionless')
                                    for k in (coordinates or ['x'])}),
                reference=dict(spin_frequency_hz=10 if coordinates else None,
                               frame='native PE reference; no spin evolution'),
                cosmology='Planck15_lal' if coordinates else 'not applicable',
                prior=dict(formula='DL^2 dDL/dz * m1_source*(1+z)^2/4',
                           normalization='hyperparameter-independent constants omitted'),
                support_status='event-specific global support not fully reconstructed',
                identity_status='exact source rows' if coordinates else 'synthetic')


def main():
    source=dict(file='synthetic', url='urn:synthetic:stage6', sha256='0'*64)
    m=metadata(source)
    m['prior']=dict(formula='[0.2,0.3,0.5]', normalization='synthetic density')
    ev=WeightedEvent({'x':np.array([.1,.4,.8])},np.array([4,8,12]),
                     np.log([.2,.3,.5]),np.log([1.,2.,3.]),m)
    target=np.array([0.,.6,.4])
    with np.errstate(divide='ignore'):
        lt=np.log(target)
    expected=np.sum(np.array([1.,2.,3.])*target/np.array([.2,.3,.5]))/6
    assert np.isclose(np.exp(ev.log_integral(lt)),expected,rtol=1e-14)
    new=ev.reexpress_prior(np.log([.4,.4,.2]),dict(formula='synthetic new'), 'fixture support known')
    # Re-expression changes only the event-wide normalizing constant.
    targets=np.column_stack([target,[.7,.1,.2],[.2,.2,.6]])
    with np.errstate(divide='ignore'):
        log_targets=np.log(targets)
    delta=new.log_integral(log_targets)-ev.log_integral(log_targets)
    assert np.ptp(delta)<1e-13
    diag=ev.diagnostics(lt)
    fractions=np.array([0.,4.,2.4])/6.4
    assert np.isclose(diag['contribution_ess'],1/np.sum(fractions**2))
    assert diag['largest_source_row']==8
    shifted=WeightedEvent(ev.samples,ev.row_ids,ev.log_pe_prior,ev.log_sample_weight+10000,m)
    assert np.allclose(shifted.log_integral(log_targets),ev.log_integral(log_targets),rtol=0,atol=2e-12)
    assert np.isneginf(ev.log_integral(np.full(3,-np.inf)))
    rejected=[]
    cases={
        'zero new prior':lambda: ev.reexpress_prior(np.array([0.,-np.inf,0.]),{},'unknown'),
        'duplicate source rows':lambda: WeightedEvent(ev.samples,np.array([1,1,2]),ev.log_pe_prior,ev.log_sample_weight,m),
        'all-zero representation weights':lambda: WeightedEvent(ev.samples,ev.row_ids,ev.log_pe_prior,np.full(3,-np.inf),m),
        'NaN target':lambda: ev.log_integral(np.array([np.nan,0.,0.])),
        'zero-integral diagnostics':lambda: ev.diagnostics(np.full(3,-np.inf)),
        'coordinate mismatch':lambda: WeightedEvent({'y':np.ones(3)},ev.row_ids,ev.log_pe_prior,ev.log_sample_weight,m),
    }
    for name,func in cases.items():
        try:func()
        except ValueError:rejected.append(name)
        else:raise AssertionError(name)
    with tempfile.TemporaryDirectory(dir=BASE) as folder:
        path=Path(folder)/'roundtrip.npz';new.save(path);loaded=WeightedEvent.load(path)
        assert loaded.metadata==new.metadata and np.array_equal(loaded.row_ids,new.row_ids)
        assert np.array_equal(loaded.log_sample_weight,new.log_sample_weight)
        assert np.array_equal(loaded.log_integral(log_targets),new.log_integral(log_targets))

    init_model()
    download=json.loads((ROOT/'research/population-stage5/download.json').read_text())
    raw=ROOT/'research/population-stage5/source-data'/download['file']
    assert hashlib.sha256(raw.read_bytes()).hexdigest()==download['sha256']
    stage5=json.loads((ROOT/'research/population-stage5/tail-results.json').read_text())
    frozen=json.loads((ROOT/'results/model.json').read_text())
    with np.load(ROOT/'research/population-stage5/representations.npz',allow_pickle=False) as bank:
        d={k:bank[f'nocosmo_{k}'] for k in PARAMS}
    pc,pn,_=density_pair(d)
    ev=WeightedEvent(d,np.arange(len(pn)),np.log(pn),np.zeros(len(pn)),
                     metadata({k:download[k] for k in ['file','url','sha256']},'GW190805_211137',PARAMS))
    converted=ev.reexpress_prior(np.log(pc),dict(formula='dVc/dz/(1+z)*m1*(1+z)^2/4',
                                                normalization='constants omitted'),'global support unresolved; same retained rows')
    # Unit prior obtains population density independently of the interface.
    population=weights(dict(d,prior=np.ones(len(pn))))
    with np.errstate(divide='ignore'):
        lp=np.log(population)
    a=ev.log_integral(lp);b=converted.log_integral(lp)
    direct=np.mean(population/pn[:,None,None],axis=0)
    error=float(np.max(abs(np.exp(a)/direct-1)));assert error<1e-11
    shape_error=float(np.ptp(b-a));assert shape_error<1e-11
    q,_=summarize(frozen,np.exp(a));qc,_=summarize(frozen,np.exp(b))
    previous=np.array(stage5['scenarios']['nocosmo_lal']['quantiles'])
    qerr=float(np.max(abs(q-previous)));assert qerr<1e-11
    assert np.max(abs(q-qc))<1e-11
    from population_functionals import aligned_cdf
    f=(.5-stage5['reference_r'][0])/(.5-float(aligned_cdf(0,.3)))
    mix=f*population[:,0,0]+(1-f)*population[:,0,1]
    with np.errstate(divide='ignore'):
        dm=ev.diagnostics(np.log(mix));dc=converted.diagnostics(np.log(mix))
    assert np.isclose(dm['contribution_ess'],stage5['weight_diagnostics']['nocosmo_lal'][0]['ess'],rtol=1e-11)
    assert np.isclose(dm['contribution_ess'],dc['contribution_ess'],rtol=1e-11)
    ev.save(BASE/'gw190805-weighted-event.npz')
    result=dict(schema='blackhole-weighted-event-validation/1',status='pass',
                synthetic_direct_sum=True,nonuniform_weights=True,log_weight_shift_10000=True,
                zero_population_density=True,rejected_invalid_cases=rejected,serialization_roundtrip=True,
                gw190805_direct_integral_max_relative_error=error,
                reexpression_log_constant_range=shape_error,
                stage5_max_quantile_error=qerr,original_representation=dm,reexpressed_representation=dc,
                event_metadata=ev.metadata,source_sha256=download['sha256'],
                implementation_sha256=hashlib.sha256((ROOT/'research/population_event.py').read_bytes()).hexdigest(),
                limit='Regression and algebraic checks; no missing-support recovery, new samples or convergence proof.')
    (BASE/'interface-validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='event_metadata'},indent=2))


if __name__=='__main__':main()
