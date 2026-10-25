"""Retain, then diagnose every held-out row failing numerical tolerances."""
from pathlib import Path
import json
import warnings
from diagnose_initialization import initial_state
from transport import to_infinity

BASE=Path(__file__).resolve().parent


def main():
    prior=json.loads((BASE/'reference-validation.json').read_text())
    rows=[r for r in prior['records'] if r['tolerance_difference']>1e-6 or r['transport']['roundtrip_kappa_abs_error']>1e-7]
    records=[];failures=[]
    for r in rows:
        for rtol,atol in [(1e-12,1e-14),(1e-13,1e-15),(3e-14,3e-16)]:
            try:
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always')
                    out=to_infinity(**initial_state(r['inputs'],'L_2PN_nonspin'),solver='DOP853',rtol=rtol,atol=atol,roundtrip=True)
                records.append(dict(index=r['index'],source_row=r['source_row'],q=r['inputs']['q'],rtol=rtol,atol=atol,result=out,
                    released_difference=[out[f'cos_tilt_{j}_infinity']-r['released_cosines'][j-1] for j in [1,2]],
                    warnings=sorted({str(w.message) for w in caught if not issubclass(w.category,DeprecationWarning)})))
            except Exception as exc:failures.append(dict(index=r['index'],source_row=r['source_row'],rtol=rtol,error=str(exc)))
        print(r['index'],r['source_row'],'done',flush=True)
    result=dict(schema='blackhole-near-equal-mass-probe/1',trigger='tolerance_difference>1e-6 or roundtrip kappa error>1e-7',records=records,failures=failures)
    (BASE/'equal-mass-results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    for r in records:print(r['index'],r['q'],r['rtol'],r['released_difference'],r['result']['roundtrip_kappa_abs_error'],flush=True)


if __name__=='__main__':main()
