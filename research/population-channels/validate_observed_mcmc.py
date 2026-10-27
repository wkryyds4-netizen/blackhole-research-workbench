"""Second-algorithm checks for every complete conditional branching fit."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import hashlib
import json
import numpy as np
from independent_mcmc import sample
from validate_independent_mcmc import check

BASE=Path(__file__).resolve().parent
def sha(p):
    with p.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def save(p,d):
    temporary=p.with_suffix('.tmp');temporary.write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8');temporary.replace(p)


def main():
    fits=json.loads((BASE/'branching-results.json').read_text());assert len(fits['cases'])==24
    calibration=json.loads((BASE/'independent-mcmc-calibration.json').read_text());assert calibration['status']=='pass'
    for name,digest in calibration['input_sha256'].items():assert sha(BASE/name)==digest
    design=calibration['design'];gates=design['gates']
    channel_design=json.loads((BASE/'likelihood-design.json').read_text())
    normalizers=json.loads((BASE/'channel-log-normalization.json').read_text())
    norm=np.array([next(r['fraction'] for r in normalizers['channels'] if r['channel']==c) for c in channel_design['channels']])
    events=json.loads((BASE/'channel-log-event-integrals.json').read_text())['events']
    selection=json.loads((BASE/'channel-log-selection-integrals.json').read_text())
    hashes={p.name:sha(p) for p in [Path(__file__),BASE/'independent_mcmc.py',BASE/'validate_independent_mcmc.py',BASE/'independent-mcmc-calibration.json',BASE/'branching-results.json',BASE/'channel-log-integrals-complete.json']}
    protocol=hashlib.sha256(json.dumps(hashes,sort_keys=True).encode()).hexdigest()
    path=BASE/'observed-mcmc-validation.json'
    report=json.loads(path.read_text()) if path.exists() else dict(protocol=protocol,status='partial',cases={},input_sha256=hashes)
    assert report['protocol']==protocol
    for ci,(key,fit) in enumerate(fits['cases'].items()):
        if key in report['cases']:
            assert report['cases'][key]['status']=='pass','A preserved failed case requires an explicit diagnosed refinement, not silent skipping.'
            continue
        fi=channel_design['hidden_spin_conditionals'].index(fit['family']);eps=fit['epsilon']
        # Independently construct normalized/regularized event and selection
        # columns rather than importing the first algorithm's transform.
        raw=np.array([r['mean'][fi] for r in events])
        rs=np.array(selection['early']['total'])[fi]+np.array(selection['banks'][fit['bank']]['policies'][fit['policy']]['total'])[fi]
        means=(1-eps)*raw[:,:5]/norm[None,:]+eps*raw[:,5,None]
        efficiency=(1-eps)*rs[:5]/norm+eps*rs[5]
        result=sample(means,efficiency,design,design['seed']+100+ci)
        # Preserve diagnostics even when a fixed gate fails.
        difference=float(np.max(abs(np.array(result['intrinsic_quantiles'])-fit['intrinsic_quantiles'])))
        report['cases'][key]=dict(**result,maximum_qmc_quantile_difference=difference,status='not_yet_gated')
        save(path,report)
        check(result,gates)
        assert difference<gates['maximum_observed_qmc_quantile_difference'],(key,'quantiles',difference)
        report['cases'][key]['status']='pass';save(path,report)
        print(json.dumps(dict(case=key,maximum_qmc_difference=difference,rhat=result['maximum_rhat'],bulk_ess=result['minimum_bulk_ess'],tail_ess=result['minimum_tail_ess'])),flush=True)
    assert len(report['cases'])==24 and all(r['status']=='pass' for r in report['cases'].values())
    report.update(status='pass',scope='Independent Metropolis vs importance-Sobol posterior calculations for24 complete saved-integral cases. Calibration and convergence gates do not establish physical, PE-chain or injection-campaign coverage.')
    save(path,report)


if __name__=='__main__':
    try:main()
    except Exception as error:
        path=BASE/'observed-mcmc-failures.json';h=json.loads(path.read_text()) if path.exists() else []
        h.append(dict(tool_sha256=sha(Path(__file__)),error=repr(error)))
        save(path,h);raise
