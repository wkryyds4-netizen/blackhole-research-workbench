"""Keep primary validation and adaptive follow-up distinct in the final audit."""
from pathlib import Path
import hashlib
import json
import numpy as np

BASE=Path(__file__).resolve().parent


def main():
    primary=json.loads((BASE/'reference-validation.json').read_text())
    probe=json.loads((BASE/'equal-mass-results.json').read_text())
    assert not primary['failures'] and not probe['failures']
    refinements={}
    for r in probe['records']:refinements.setdefault((r['index'],r['source_row']),[]).append(r)
    rows=[];final_warnings=[]
    for r in primary['records']:
        group=refinements.get((r['index'],r['source_row']))
        if group:
            b,c=sorted(group,key=lambda x:x['rtol'],reverse=True)[-2:]
            final_warnings.extend(b['warnings']+c['warnings'])
            convergence=max(abs(b['result'][f'cos_tilt_{i}_infinity']-c['result'][f'cos_tilt_{i}_infinity']) for i in [1,2])
            rows.append(dict(index=r['index'],source_row=r['source_row'],refined=True,difference=c['released_difference'],
                tolerance_difference=convergence,roundtrip_kappa_error=c['result']['roundtrip_kappa_abs_error']))
        else:rows.append(dict(index=r['index'],source_row=r['source_row'],refined=False,difference=r['difference'],
            tolerance_difference=r['tolerance_difference'],roundtrip_kappa_error=r['transport']['roundtrip_kappa_abs_error']))
    metrics=dict(rows=len(rows),refined_rows=len(refinements),
        release_abs_difference_quantiles=np.quantile(abs(np.array([r['difference'] for r in rows])),[.5,.9,.99,1]).tolist(),
        max_tolerance_difference=max(r['tolerance_difference'] for r in rows),
        max_roundtrip_kappa_error=max(r['roundtrip_kappa_error'] for r in rows))
    # Engineering checks on this fixed bank, not a physical-error certificate.
    checks=dict(all_requested_completed=len(rows)==primary['summary']['requested'],
        cosine_release_difference_below_1e_minus6=metrics['release_abs_difference_quantiles'][-1]<1e-6,
        cosine_tolerance_difference_below_1e_minus6=metrics['max_tolerance_difference']<1e-6,
        roundtrip_kappa_difference_below_1e_minus8=metrics['max_roundtrip_kappa_error']<1e-8,
        final_tolerance_pairs_without_warnings=not primary['warnings'] and not final_warnings)
    result=dict(schema='blackhole-spin-reference-audit/1',metrics=metrics,checks=checks,rows=rows,
        primary_validation=primary['summary'],
        intermediate_warnings=[dict(index=r['index'],source_row=r['source_row'],rtol=r['rtol'],warnings=r['warnings']) for r in probe['records'] if r['warnings']],
        sources={name:hashlib.sha256((BASE/name).read_bytes()).hexdigest() for name in ['reference-validation.json','equal-mass-results.json','initialization-results.json','historical-metadata.json']},
        common_reference_population_ready=False,
        scope='Post-diagnosis numerical equivalence on 350 held-out rows after 5 adaptive refinements. Primary failures are retained separately. Not a hybrid evolution accuracy test, a calibrated coverage statement, or full-catalog/injection transport.')
    (BASE/'reference-audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(dict(metrics=metrics,checks=checks),indent=2))
    assert all(checks.values()),'Reference audit did not pass; inspect retained diagnostics.'


if __name__=='__main__':main()
