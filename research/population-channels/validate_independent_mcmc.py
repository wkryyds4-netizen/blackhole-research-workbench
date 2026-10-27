"""Analytic calibration before a second algorithm sees observed channel fits."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import hashlib
import json
import argparse
import numpy as np
from scipy.stats import beta as beta_distribution
from independent_mcmc import sample

BASE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def check(r,gates):
    assert r['maximum_rhat']<gates['maximum_rank_folded_split_rhat'],('rhat',r['maximum_rhat'])
    assert r['minimum_bulk_ess']>gates['minimum_bulk_ess'],('bulk_ess',r['minimum_bulk_ess'])
    assert r['minimum_tail_ess']>gates['minimum_tail_ess'],('tail_ess',r['minimum_tail_ess'])


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--refined',action='store_true');args=parser.parse_args()
    design=json.loads((BASE/'mcmc-design.json').read_text());gates=design['gates'];cases=[]
    if args.refined:
        refinement=json.loads((BASE/'mcmc-refinement.json').read_text())
        for key in ['production_steps','save_every']:design[key]=refinement[key]
    for i,(name,counts) in enumerate([('no_information',None),('known_69_assignments',[35,0,0,0,34])]):
        means=np.ones((69,5)) if counts is None else np.repeat(np.eye(5),counts,axis=0)
        alpha=np.ones(5) if counts is None else 1+np.array(counts)
        exact=np.array([beta_distribution.ppf([.05,.5,.95],a,alpha.sum()-a) for a in alpha])
        result=sample(means,np.ones(5),design,design['seed']+i)
        error=float(np.max(abs(np.array(result['intrinsic_quantiles'])-exact)))
        attempt=BASE/'independent-mcmc-attempts.json';history=json.loads(attempt.read_text()) if attempt.exists() else []
        history.append(dict(case=name,design=design,result=result,maximum_exact_quantile_error=error))
        attempt.write_text(json.dumps(history,indent=2)+'\n',encoding='utf-8')
        check(result,gates)
        assert error<gates['maximum_exact_calibration_quantile_error'],('quantiles',error)
        cases.append(dict(case=name,**result,maximum_exact_quantile_error=error,exact_quantiles=exact.tolist()))
        print(json.dumps({k:cases[-1][k] for k in ['case','maximum_exact_quantile_error','maximum_rhat','minimum_bulk_ess','minimum_tail_ess']}),flush=True)
    report=dict(status='pass',schema='blackhole-independent-mcmc-calibration/1',cases=cases,
                design=design,
                input_sha256={p.name:sha(p) for p in [BASE/'mcmc-design.json',BASE/'independent_mcmc.py',Path(__file__),*([BASE/'mcmc-refinement.json'] if args.refined else [])]},
                scope='Computational calibration only. No event, campaign or physical coverage claim.')
    (BASE/'independent-mcmc-calibration.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    try:main()
    except Exception as error:
        p=BASE/'independent-mcmc-failures.json';h=json.loads(p.read_text()) if p.exists() else []
        h.append(dict(tool_sha256=sha(Path(__file__)),module_sha256=sha(BASE/'independent_mcmc.py'),error=repr(error)))
        p.write_text(json.dumps(h,indent=2)+'\n',encoding='utf-8');raise
