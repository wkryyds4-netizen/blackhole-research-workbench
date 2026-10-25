"""Held-out validation after fixing the initial-L convention from diagnosis.

Before running: 64 random rows per event, seed 20261025, exclude discovery rows;
plus distinct rows at maximum q, spin magnitude, detector mass, and in-plane spin.
No residual-based row selection, exclusions, or fitting of coefficients.
"""
from pathlib import Path
import hashlib
import json
import warnings
import numpy as np
from diagnose_initialization import initial_state
from transport import to_infinity

BASE=Path(__file__).resolve().parent


def main():
    provenance=json.loads((BASE/'spin-provenance.json').read_text())
    discovery=json.loads((BASE/'reference-results.json').read_text())
    assert hashlib.sha256((BASE/'spin-bank.npz').read_bytes()).hexdigest()==provenance['cache_sha256']
    rng=np.random.default_rng(20261025);records=[];failures=[];warning_log=[];selection=[]
    with np.load(BASE/'spin-bank.npz',allow_pickle=False) as bank:
        for event in provenance['events']:
            index=event['index'];prefix=str(index)+'_'
            d={k[len(prefix):]:bank[k] for k in bank.files if k.startswith(prefix)}
            used=[r['source_row'] for r in discovery['transport_comparison'] if r['index']==index]
            available=np.setdiff1d(np.arange(event['rows']),used)
            random_rows=rng.choice(available,64,replace=False)
            criteria=dict(q=d['mass_ratio'],a1=d['a_1'],a2=d['a_2'],
                detector_mass=d['mass_1']+d['mass_2'],
                inplane1=d['a_1']*np.sin(d['tilt_1']),inplane2=d['a_2']*np.sin(d['tilt_2']))
            stress={name:int(available[np.argmax(v[available])]) for name,v in criteria.items()}
            chosen=np.unique(np.concatenate([random_rows,list(stress.values())]))
            selection.append(dict(index=index,random_rows=random_rows.tolist(),stress=stress,distinct_rows=len(chosen)))
            for row in chosen:
                inputs=dict(theta1=float(d['tilt_1'][row]),theta2=float(d['tilt_2'][row]),phi12=float(d['phi_12'][row]),
                    q=float(d['mass_ratio'][row]),a1=float(d['a_1'][row]),a2=float(d['a_2'][row]),
                    total_detector_mass=float(d['mass_1'][row]+d['mass_2'][row]),reference_hz=event['native_spin_reference_hz'])
                try:
                    state=initial_state(inputs,'L_2PN_nonspin')
                    with warnings.catch_warnings(record=True) as caught:
                        warnings.simplefilter('always')
                        tight=to_infinity(**state,solver='DOP853',rtol=1e-11,atol=1e-13,roundtrip=True)
                        loose=to_infinity(**state,solver='DOP853',rtol=1e-9,atol=1e-11)
                    meaningful=sorted({str(w.message) for w in caught if not issubclass(w.category,DeprecationWarning)})
                    if meaningful:warning_log.append(dict(index=index,source_row=int(row),warnings=meaningful))
                    target=[float(d[f'cos_tilt_{j}_infinity_only_prec_avg'][row]) for j in [1,2]]
                    records.append(dict(index=index,source_row=int(row),inputs=inputs,transport=tight,released_cosines=target,
                        difference=[tight[f'cos_tilt_{j}_infinity']-target[j-1] for j in [1,2]],
                        tolerance_difference=max(abs(tight[f'cos_tilt_{j}_infinity']-loose[f'cos_tilt_{j}_infinity']) for j in [1,2])))
                except Exception as exc:failures.append(dict(index=index,source_row=int(row),error=type(exc).__name__,message=str(exc)))
            print(f'{event["event"]}: {sum(r["index"]==index for r in records)}/{len(chosen)} held-out rows',flush=True)
    summary=dict(requested=sum(s['distinct_rows'] for s in selection),completed=len(records),
        abs_difference_quantiles=np.quantile(abs(np.array([r['difference'] for r in records])),[.5,.9,.99,1]).tolist(),
        max_tolerance_difference=max(r['tolerance_difference'] for r in records),
        max_roundtrip_kappa_error=max(r['transport']['roundtrip_kappa_abs_error'] for r in records),
        max_chieff_drift=max(r['transport']['chieff_transport_error'] for r in records))
    result=dict(schema='blackhole-spin-reference-validation/1',selection=selection,seed=20261025,summary=summary,
        records=records,failures=failures,warnings=warning_log,
        interpretation='Independent implementation agrees with release under matched initial-L convention on this held-out set only; no population weights, injection transport, or hybrid physical-accuracy validation.',
        source_sha256=provenance['cache_sha256'])
    (BASE/'reference-validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(dict(summary=summary,failures=failures,warnings=warning_log),indent=2),flush=True)


if __name__=='__main__':main()
