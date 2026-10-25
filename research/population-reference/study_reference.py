"""Audit available infinity tilts and test transport without inventing missing phases."""
from pathlib import Path
import hashlib
import json
import platform
import warnings
import numpy as np
import scipy
import precession
from scipy.stats import ks_2samp
from transport import to_infinity

BASE=Path(__file__).resolve().parent


def main():
    assert precession.__version__=='2.1.2'
    source=json.loads((BASE/'spin-provenance.json').read_text())
    assert hashlib.sha256((BASE/'spin-bank.npz').read_bytes()).hexdigest()==source['cache_sha256']
    summaries=[];comparison=[];warning_log=[];failures=[]
    # Counterexample fixed before the calculation; not selected for a large shift.
    proof=[]
    for phi in np.linspace(0,2*np.pi,25):
        a=to_infinity(1.1,2.,float(phi),.65,.7,.6,80.,20.,roundtrip=True)
        b=to_infinity(1.1,2.,float(phi),.65,.7,.6,80.,20.,rtol=1e-11,atol=1e-13,solver='DOP853')
        loose=to_infinity(1.1,2.,float(phi),.65,.7,.6,80.,20.,rtol=1e-9,atol=1e-11,solver='DOP853')
        a.update(phi12=float(phi),tight_cosine_max_difference=max(abs(loose[f'cos_tilt_{j}_infinity']-b[f'cos_tilt_{j}_infinity']) for j in [1,2]),
                 independent_solver_max_difference=max(abs(a[f'cos_tilt_{j}_infinity']-b[f'cos_tilt_{j}_infinity']) for j in [1,2]))
        proof.append(a)
    with np.load(BASE/'spin-bank.npz',allow_pickle=False) as bank:
        for event in source['events']:
            index=event['index'];d={k[len(str(index))+1:]:bank[k] for k in bank.files if k.startswith(str(index)+'_')}
            n=event['rows'];record=dict(event=event['event'],index=index,rows=n,native_reference_hz=event['native_spin_reference_hz'],spins=[])
            for j in [1,2]:
                native=d[f'cos_tilt_{j}'];infinite=d[f'cos_tilt_{j}_infinity_only_prec_avg'];valid=np.isfinite(native)&np.isfinite(infinite)
                # Report completeness; no undisclosed dropping of failed conversions.
                record['spins'].append(dict(spin=j,finite_rows=int(valid.sum()),
                    native_negative_fraction=float(np.mean(native[valid]<0)),infinity_negative_fraction=float(np.mean(infinite[valid]<0)),
                    sign_change_fraction=float(np.mean((native[valid]<0)!=(infinite[valid]<0))),
                    absolute_cosine_change_quantiles=np.quantile(abs(native[valid]-infinite[valid]),[.5,.9,.99,1]).tolist(),
                    empirical_cdf_max_difference=float(ks_2samp(native[valid],infinite[valid],method='asymp').statistic)))
            summaries.append(record)
            for row in np.linspace(0,n-1,16,dtype=int):
                kwargs=dict(theta1=float(d['tilt_1'][row]),theta2=float(d['tilt_2'][row]),phi12=float(d['phi_12'][row]),
                    q=float(d['mass_ratio'][row]),a1=float(d['a_1'][row]),a2=float(d['a_2'][row]),
                    total_detector_mass=float(d['mass_1'][row]+d['mass_2'][row]),reference_hz=event['native_spin_reference_hz'])
                try:
                    with warnings.catch_warnings(record=True) as caught:
                        warnings.simplefilter('always')
                        value=to_infinity(**kwargs,roundtrip=True)
                        tight=to_infinity(**kwargs,rtol=1e-11,atol=1e-13,solver='DOP853',roundtrip=True)
                        loose=to_infinity(**kwargs,rtol=1e-9,atol=1e-11,solver='DOP853')
                    meaningful=sorted({str(w.message) for w in caught if not issubclass(w.category,DeprecationWarning)})
                    if meaningful:warning_log.append(dict(index=index,source_row=int(row),warnings=meaningful))
                    value.update(index=index,source_row=int(row),inputs=kwargs,
                        released_cosines=[float(d[f'cos_tilt_{j}_infinity_only_prec_avg'][row]) for j in [1,2]],
                        implementation_difference=[value[f'cos_tilt_{j}_infinity']-float(d[f'cos_tilt_{j}_infinity_only_prec_avg'][row]) for j in [1,2]],
                        tight_cosine_max_difference=max(abs(loose[f'cos_tilt_{j}_infinity']-tight[f'cos_tilt_{j}_infinity']) for j in [1,2]),
                        independent_solver_max_difference=max(abs(value[f'cos_tilt_{j}_infinity']-tight[f'cos_tilt_{j}_infinity']) for j in [1,2]),
                        independent_roundtrip_kappa_error=tight['roundtrip_kappa_abs_error'])
                    comparison.append(value)
                except Exception as exc:
                    failures.append(dict(index=index,source_row=int(row),error=type(exc).__name__,message=str(exc)))
            print(f'{event["event"]}: {sum(x["index"]==index for x in comparison)}/16 transport rows completed',flush=True)
    validations=dict(proof_max_roundtrip_kappa_error=max(x['roundtrip_kappa_abs_error'] for x in proof),
        proof_max_tight_cosine_difference=max(x['tight_cosine_max_difference'] for x in proof),
        proof_primary_cosine_range=[min(x['cos_tilt_1_infinity'] for x in proof),max(x['cos_tilt_1_infinity'] for x in proof)],
        proof_secondary_cosine_range=[min(x['cos_tilt_2_infinity'] for x in proof),max(x['cos_tilt_2_infinity'] for x in proof)],
        sample_max_roundtrip_kappa_error=max(x['roundtrip_kappa_abs_error'] for x in comparison),
        sample_max_tight_cosine_difference=max(x['tight_cosine_max_difference'] for x in comparison),
        sample_max_independent_solver_difference=max(x['independent_solver_max_difference'] for x in comparison),
        sample_max_independent_roundtrip_kappa_error=max(x['independent_roundtrip_kappa_error'] for x in comparison),
        sample_max_chieff_drift=max(x['chieff_transport_error'] for x in comparison),
        released_implementation_abs_difference_quantiles=np.quantile(abs(np.array([x['implementation_difference'] for x in comparison])),[.5,.9,.99,1]).tolist())
    result=dict(schema='blackhole-spin-reference-study/1',environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,precession=precession.__version__),
        spin_bank_sha256=source['cache_sha256'],events=summaries,
        counterexample_fixed_inputs=dict(theta1=1.1,theta2=2.,q=.65,a1=.7,a2=.6,total_detector_mass=80.,reference_hz=20.),
        counterexample=proof,transport_comparison=comparison,failures=failures,warnings=warning_log,validation=validations,
        common_reference_population_ready=False,
        upstream_tolerance_override='precession 2.1.2 integrator_precav tests rol/aol, not rtol/atol, and resets tolerances to 1e-10. Independent DOP853 integration of its documented RHS performs the actual tolerance sweep; upstream source unmodified.',
        limits='Event posterior summaries only; no population weighting or formation-channel inference. Infinity approximations differ between library and historical release; no assumption of equality. Missing tutorial phi12 and unresolved injection native reference prevent deterministic common-reference population transport.')
    (BASE/'reference-results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(dict(events=summaries,validation=validations,failures=failures,warnings=warning_log),indent=2),flush=True)


if __name__=='__main__':main()
