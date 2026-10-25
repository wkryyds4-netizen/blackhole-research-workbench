"""Retain a failed release-agreement gate and test independent solver stability."""
from pathlib import Path
import json
import sys
import warnings
import numpy as np
import precession
from scipy.integrate import solve_ivp

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
sys.path.insert(0,str(ROOT/'research/population-reference'))
from diagnose_initialization import initial_state


def main():
    source=json.loads((BASE/'alias-heldout-input.json').read_text());validation=json.loads((BASE/'alias-heldout-validation.json').read_text())
    worst=max(validation['rows'],key=lambda v:max(map(abs,v['difference'])))['source_row']
    near_equal=next(v['source_row'] for v in source['rows'] if v['reason']=='maximum q')
    rows=[];summaries=[]
    for index in sorted({worst,near_equal}):
        record=next(r for r in source['rows'] if r['source_row']==index)
        for freq in [16.,20.]:
            inputs=initial_state(dict(record['inputs'],reference_hz=freq),'L_2PN_nonspin')
            q=inputs['q'];a1=inputs['a1'];a2=inputs['a2'];r0=inputs['initial_separation']
            _,kappa,ce=precession.angles_to_conserved(theta1=inputs['theta1'],theta2=inputs['theta2'],deltaphi=inputs['phi12'],q=q,chi1=a1,chi2=a2,r=r0)
            ce=float(ce.item());kappa=float(kappa.item());u0=float(precession.eval_u(r0,q).item())
            def rhs(u,k):return [float(np.asarray(precession.rhs_precav(k,u,ce,q,a1,a2)).item())]
            group=[]
            for method in ['DOP853','RK45']:
                for tolerance in [1e-13,3e-14]:
                    with warnings.catch_warnings(record=True) as caught:
                        warnings.simplefilter('always');answer=solve_ivp(rhs,(u0,0),[kappa],method=method,rtol=tolerance,atol=tolerance/100)
                    assert answer.success
                    dc=precession.eval_deltachiinf(float(answer.y[0,-1]),ce,q,a1,a2)
                    cosines=[float(precession.eval_costheta1(dc,ce,q,a1).item()),float(precession.eval_costheta2(dc,ce,q,a2).item())]
                    target=record['canonical_cosines'] if freq==16 else record['alias_cosines']
                    row=dict(source_row=index,frequency=freq,method=method,rtol=tolerance,cosines=cosines,target_difference=(np.array(cosines)-target).tolist(),rhs_evaluations=answer.nfev,
                        warnings=sorted({str(w.message) for w in caught if not issubclass(w.category,DeprecationWarning)}));rows.append(row);group.append(row)
            spread=float(np.max(np.ptp(np.array([v['cosines'] for v in group]),axis=0)))
            summaries.append(dict(source_row=index,frequency=freq,two_solver_two_tolerance_spread=spread,max_release_residual=max(max(map(abs,v['target_difference'])) for v in group)))
    report=dict(schema='blackhole-alias-residual-refinement/1',initial_release_agreement_gate_passed=validation['numerical_checks_passed'],summaries=summaries,rows=rows,
        interpretation='Two ODE solvers test whether the small retained release discrepancy is integration instability. A converged residual still exceeding 1e-6 does not pass the original gate; do not overwrite that failure or call it exact release reproduction.')
    (BASE/'alias-residual-refinement.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(dict(initial_release_agreement_gate_passed=report['initial_release_agreement_gate_passed'],summaries=summaries),indent=2))


if __name__=='__main__':main()
