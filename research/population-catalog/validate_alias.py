"""Held-out frequency-convention validation; prepare and integrate separately."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import warnings
import numpy as np

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]


def prepare():
    import h5py
    sys.path.insert(0,str(ROOT/'research/population-stage7'))
    from extract_events import join,SOURCE
    diagnosis=json.loads((BASE/'alias-source-diagnosis.json').read_text());excluded={r['canonical_source_row'] for r in diagnosis['diagnostic_rows']}
    rec=next(r for r in json.loads((BASE/'download-results.json').read_text()) if r['index']==14);path=BASE/'source-data'/rec['source']['file']
    assert hashlib.sha256(path.read_bytes()).hexdigest()==rec['sha256']
    with h5py.File(path,'r') as h:
        a=h[diagnosis['canonical']+'/posterior_samples'][:];b=h[diagnosis['alias']+'/posterior_samples'][:]
        mapping,info=join([a[k] for k in SOURCE],[b[k] for k in SOURCE]);assert info['matches']==len(a) and info['ambiguous']==0
    rng=np.random.default_rng(20261025);allowed=np.array(sorted(set(range(len(a)))-excluded))
    chosen=list(rng.choice(allowed,64,replace=False));reasons={int(i):'random' for i in chosen}
    extremes=[('q',a['mass_ratio']),('a1',a['a_1']),('a2',a['a_2']),('Mdet',a['mass_1']+a['mass_2']),('a1_perp',a['a_1']*np.sqrt(1-a['cos_tilt_1']**2)),('a2_perp',a['a_2']*np.sqrt(1-a['cos_tilt_2']**2))]
    for name,v in extremes:
        i=next(int(i) for i in np.argsort(v)[::-1] if int(i) not in excluded and int(i) not in chosen);chosen.append(i);reasons[i]='maximum '+name
    assert len(set(chosen))==70 and not set(chosen)&excluded
    rows=[]
    for i in chosen:
        rows.append(dict(source_row=int(i),reason=reasons[int(i)],inputs=dict(theta1=float(np.arccos(a['cos_tilt_1'][i])),theta2=float(np.arccos(a['cos_tilt_2'][i])),phi12=float(a['phi_12'][i]),q=float(a['mass_ratio'][i]),a1=float(a['a_1'][i]),a2=float(a['a_2'][i]),total_detector_mass=float(a['mass_1'][i]+a['mass_2'][i]),reference_hz=16.),
            canonical_cosines=[float(a[f'cos_tilt_{j}_infinity_only_prec_avg'][i]) for j in [1,2]],alias_cosines=[float(b[f'cos_tilt_{j}_infinity_only_prec_avg'][mapping[i]]) for j in [1,2]]))
    payload=dict(schema='blackhole-alias-heldout-input/1',seed=20261025,diagnosis_sha256=hashlib.sha256((BASE/'alias-source-diagnosis.json').read_bytes()).hexdigest(),rows=rows,
        selection='64 random rows excluding diagnostics, then six distinct maxima of native inputs, selected independently of derived target columns.')
    (BASE/'alias-heldout-input.json').write_text(json.dumps(payload,indent=2),encoding='utf-8');print('Prepared 70 held-out rows.',flush=True)


def integrate():
    sys.path.insert(0,str(ROOT/'research/population-reference'))
    from transport import to_infinity
    from diagnose_initialization import initial_state
    source=json.loads((BASE/'alias-heldout-input.json').read_text());rows=[];refinements=[];warnings_log=[]
    def run(record,freq,tolerance):
        inputs=dict(record['inputs'],reference_hz=freq)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always');answer=to_infinity(**initial_state(inputs,'L_2PN_nonspin'),solver='DOP853',rtol=tolerance,atol=tolerance/100,roundtrip=True)
        messages=sorted({str(w.message) for w in caught if not issubclass(w.category,DeprecationWarning)})
        if messages:warnings_log.append(dict(row=record['source_row'],frequency=freq,rtol=tolerance,messages=messages))
        return answer
    for n,record in enumerate(source['rows']):
        for freq in [16.,20.]:
            results=[run(record,freq,t) for t in [1e-11,1e-12]];target=record['canonical_cosines'] if freq==16 else record['alias_cosines']
            def cos(a):return np.array([a[f'cos_tilt_{j}_infinity'] for j in [1,2]])
            delta=float(np.max(abs(cos(results[-1])-cos(results[-2]))))
            if delta>1e-6 or max(a['roundtrip_kappa_abs_error'] for a in results)>1e-8 or np.max(abs(cos(results[-1])-target))>1e-6:
                results += [run(record,freq,t) for t in [1e-13,3e-14]]
                refinements.append(dict(row=record['source_row'],frequency=freq,all_results=results))
            rows.append(dict(source_row=record['source_row'],reason=record['reason'],frequency=freq,
                difference=(cos(results[-1])-target).tolist(),final_tolerance_difference=float(np.max(abs(cos(results[-1])-cos(results[-2])))),
                roundtrip_kappa_error=results[-1]['roundtrip_kappa_abs_error']))
        if (n+1)%10==0:print(f'{n+1}/70 held-out source rows',flush=True)
    summaries=[]
    for freq in [16.,20.]:
        group=[r for r in rows if r['frequency']==freq]
        summaries.append(dict(frequency=freq,rows=len(group),max_target_cosine_difference=float(np.max(abs(np.array([r['difference'] for r in group])))),
            max_final_tolerance_difference=max(r['final_tolerance_difference'] for r in group),max_roundtrip_kappa_error=max(r['roundtrip_kappa_error'] for r in group)))
    passed=all(v['max_target_cosine_difference']<1e-6 and v['max_final_tolerance_difference']<1e-6 and v['max_roundtrip_kappa_error']<1e-8 for v in summaries)
    report=dict(schema='blackhole-alias-heldout-validation/1',input_sha256=hashlib.sha256((BASE/'alias-heldout-input.json').read_bytes()).hexdigest(),summaries=summaries,rows=rows,
        refinements=refinements,warnings=warnings_log,numerical_checks_passed=passed,
        limits='Validates numerical compatibility of the two released labels with 16 and 20 Hz conventions on held-out/stress rows. Does not recover the historical call configuration or justify substituting either label into a population without a consistent selection function.')
    (BASE/'alias-heldout-validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(dict(summaries=summaries,refinements=len(refinements),warnings=warnings_log,passed=passed),indent=2))
    if not passed:raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--integrate',action='store_true');args=parser.parse_args()
    if args.integrate:integrate()
    else:prepare()
