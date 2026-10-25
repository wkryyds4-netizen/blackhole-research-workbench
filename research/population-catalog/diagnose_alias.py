"""Native-posterior equivalence and derived-field differences between PE labels.

Extraction uses h5py. Evolution uses the separately pinned precession runtime.
Reference-frequency alternatives are diagnostic controls, not asserted history.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import warnings
import numpy as np

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
CANONICAL='C01:IMRPhenomXPHM'
ALIAS='C01:Mixed'
SOURCE=['mass_1_source','mass_ratio','redshift','a_1','a_2','cos_tilt_1','cos_tilt_2']
NATIVE=SOURCE+['mass_1','mass_2','mass_2_source','luminosity_distance','chirp_mass','phi_12','phi_jl','theta_jn','psi','phase','ra','dec','log_likelihood','chi_eff']+[f'spin_{i}{ax}' for i in [1,2] for ax in 'xyz']


def extract():
    import h5py
    sys.path.insert(0,str(ROOT/'research/population-stage7'))
    from extract_events import join
    rec=next(r for r in json.loads((BASE/'download-results.json').read_text()) if r['index']==14)
    path=BASE/'source-data'/rec['source']['file'];assert hashlib.sha256(path.read_bytes()).hexdigest()==rec['sha256']
    with h5py.File(path,'r') as h:
        a=h[CANONICAL+'/posterior_samples'][:];b=h[ALIAS+'/posterior_samples'][:]
        mapping,info=join([a[k] for k in SOURCE],[b[k] for k in SOURCE])
        assert info['matches']==len(a)==len(b) and info['ambiguous']==0 and len(np.unique(mapping))==len(a)
        equal={k:bool(np.array_equal(a[k],b[k][mapping],equal_nan=True)) for k in NATIVE};assert all(equal.values())
        different={}
        for k in sorted(set(a.dtype.names)&set(b.dtype.names)):
            aa=a[k];bb=b[k][mapping]
            if np.array_equal(aa,bb,equal_nan=True):continue
            finite=np.isfinite(aa)&np.isfinite(bb);delta=abs(aa[finite]-bb[finite])
            different[k]=dict(max_abs=float(delta.max()),median_abs=float(np.median(delta)),finite_rows=int(finite.sum()))
        fref=float(h[CANONICAL+'/meta_data/meta_data/f_ref'][0]);assert fref==16
        all_delta=np.maximum(abs(a['cos_tilt_1_infinity_only_prec_avg']-b['cos_tilt_1_infinity_only_prec_avg'][mapping]),abs(a['cos_tilt_2_infinity_only_prec_avg']-b['cos_tilt_2_infinity_only_prec_avg'][mapping]))
        diagnostic=np.unique(np.r_[np.linspace(0,len(a)-1,16,dtype=int),np.argsort(all_delta)[-4:]])
        records=[]
        for i in diagnostic:
            records.append(dict(canonical_source_row=int(i),alias_source_row=int(mapping[i]),
                inputs=dict(theta1=float(np.arccos(a['cos_tilt_1'][i])),theta2=float(np.arccos(a['cos_tilt_2'][i])),phi12=float(a['phi_12'][i]),q=float(a['mass_ratio'][i]),a1=float(a['a_1'][i]),a2=float(a['a_2'][i]),total_detector_mass=float(a['mass_1'][i]+a['mass_2'][i]),reference_hz=fref),
                canonical_infinity=[float(a[f'cos_tilt_{j}_infinity_only_prec_avg'][i]) for j in [1,2]],
                alias_infinity=[float(b[f'cos_tilt_{j}_infinity_only_prec_avg'][mapping[i]]) for j in [1,2]]))
    report=dict(schema='blackhole-label-alias-diagnosis/1',index=14,event=rec['event'],source_sha256=rec['sha256'],canonical=CANONICAL,alias=ALIAS,
        full_native_join=info,native_columns_exact=equal,native_posterior_equivalent=True,
        different_derived_fields=different,diagnostic_rows=records,
        note='All native input columns tested agree after an exact bijective row permutation. Derived columns differ. Labels are not asserted identical in all fields; evolution controls will diagnose possible conventions.')
    (BASE/'alias-source-diagnosis.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='diagnostic_rows'},indent=2))


def evolve():
    sys.path.insert(0,str(ROOT/'research/population-reference'))
    from transport import to_infinity
    from diagnose_initialization import initial_state
    source=json.loads((BASE/'alias-source-diagnosis.json').read_text());rows=[];failures=[];warning_log=[]
    for record in source['diagnostic_rows']:
        for freq in [16.,20.]:
            try:
                inputs=dict(record['inputs'],reference_hz=freq)
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always')
                    result=to_infinity(**initial_state(inputs,'L_2PN_nonspin'),solver='DOP853',rtol=1e-11,atol=1e-13,roundtrip=True)
                meaningful=sorted({str(w.message) for w in caught if not issubclass(w.category,DeprecationWarning)})
                if meaningful:warning_log.append(dict(row=record['canonical_source_row'],frequency=freq,warnings=meaningful))
                values=np.array([result[f'cos_tilt_{j}_infinity'] for j in [1,2]])
                rows.append(dict(source_row=record['canonical_source_row'],frequency=freq,
                    difference_from_canonical=(values-record['canonical_infinity']).tolist(),difference_from_alias=(values-record['alias_infinity']).tolist(),result=result))
            except Exception as error:failures.append(dict(row=record['canonical_source_row'],frequency=freq,error=str(error)))
    summaries=[]
    for freq in [16.,20.]:
        group=[r for r in rows if r['frequency']==freq]
        if group:summaries.append(dict(frequency=freq,rows=len(group),max_difference_from_canonical=float(np.max(abs(np.array([r['difference_from_canonical'] for r in group])))),
            max_difference_from_alias=float(np.max(abs(np.array([r['difference_from_alias'] for r in group])))),
            max_roundtrip_kappa_error=max(r['result']['roundtrip_kappa_abs_error'] for r in group)))
    report=dict(schema='blackhole-label-alias-evolution-diagnosis/1',source_diagnosis_sha256=hashlib.sha256((BASE/'alias-source-diagnosis.json').read_bytes()).hexdigest(),
        summaries=summaries,rows=rows,failures=failures,warnings=warning_log,
        scope='Exploratory reference-frequency controls on 16 uniform-spaced plus four largest-disagreement source rows. Numerical agreement is not recovered historical execution metadata.')
    (BASE/'alias-evolution-diagnosis.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(dict(summaries=summaries,failures=failures,warnings=warning_log),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--evolve',action='store_true');args=parser.parse_args()
    if args.evolve:evolve()
    else:extract()
