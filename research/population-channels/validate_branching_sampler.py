"""Exact conjugate/no-information controls before observed branching fits."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.special import gammaln
from scipy.stats import beta as beta_distribution
from branching_posterior import sample

BASE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    design=json.loads((BASE/'branching-sampler-design.json').read_text());gates=design['gates'];cases=[]
    for name,counts in [('no_information',None),('known_69_assignments',[35,0,0,0,34])]:
        if counts is None:means=np.ones((69,5));expected=np.ones(5);evidence=0.
        else:
            means=np.repeat(np.eye(5),counts,axis=0);expected=1+np.array(counts)
            evidence=gammaln(5)+gammaln(expected).sum()-gammaln(expected.sum())
        trueq=np.array([beta_distribution.ppf([.05,.5,.95],a,expected.sum()-a) for a in expected])
        rows=[]
        for seed in design['replicate_seeds']:
            r=sample(means,np.ones(5),design['pilot_power'],design['final_power'],seed)
            error=float(np.max(abs(np.array(r['intrinsic_quantiles'])-trueq)))
            log_error=abs(r['log_evidence']-evidence)
            assert error<gates['analytic_calibration_maximum_quantile_error'],(name,'quantiles',error)
            assert log_error<gates['analytic_calibration_log_evidence_error'],(name,'evidence',log_error)
            assert r['concentration_ess']>gates['minimum_concentration_ess_per_replicate'],(name,'ess',r['concentration_ess'])
            rows.append(dict(**r,maximum_exact_quantile_error=error,exact_log_evidence_error=log_error))
            print(json.dumps(dict(case=name,seed=seed,error=error,log_evidence_error=log_error,ess=r['concentration_ess'])),flush=True)
        span=float(np.ptp(np.array([r['intrinsic_quantiles'] for r in rows]),axis=0).max())
        assert span<gates['replicate_maximum_quantile_span']
        cases.append(dict(case=name,exact_dirichlet=expected.tolist(),exact_quantiles=trueq.tolist(),replicate_quantile_span=span,replicates=rows))
    report=dict(schema='blackhole-channel-branching-sampler-calibration/1',status='pass',cases=cases,
                input_sha256={p.name:sha(p) for p in [Path(__file__),BASE/'branching_posterior.py',BASE/'branching-sampler-design.json']},
                limits='Exact special-case computational calibration, not real-event, PE-chain, injection-campaign or astrophysical coverage calibration.')
    (BASE/'branching-sampler-validation.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    try:main()
    except Exception as error:
        p=BASE/'branching-sampler-failures.json';h=json.loads(p.read_text()) if p.exists() else []
        h.append(dict(tool_sha256=sha(Path(__file__)),module_sha256=sha(BASE/'branching_posterior.py'),error=repr(error)))
        p.write_text(json.dumps(h,indent=2)+'\n',encoding='utf-8');raise
