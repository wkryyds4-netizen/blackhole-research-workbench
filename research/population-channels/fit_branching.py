"""Fit all declared conditional five-channel cases after full69 integration.

These are exploratory restricted-domain fits, not calibrated physical rates.
Row-weight verification and analytic sampler calibration are mandatory gates.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.special import logsumexp
from branching_posterior import sample,weighted_quantiles,log_likelihood

BASE=Path(__file__).resolve().parent
def sha(p):
    with p.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def transform(normalizers,epsilon):
    result=np.zeros((6,5));result[:5]=np.diag((1-epsilon)/normalizers);result[5]=epsilon
    return result


def write(path,value):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8');tmp.replace(path)


def main():
    complete=json.loads((BASE/'channel-log-integrals-complete.json').read_text())
    assert complete['events']==69 and complete['rows']==11180578
    for name,digest in complete['input_sha256'].items():assert sha(BASE/name)==digest
    product_check=json.loads((BASE/'channel-product-validation.json').read_text());assert product_check['status']=='pass'
    for name,digest in product_check['input_sha256'].items():assert sha(BASE/name)==digest
    design=json.loads((BASE/'likelihood-design.json').read_text())
    sampler=json.loads((BASE/'branching-sampler-design.json').read_text());gates=sampler['gates']
    calibration=json.loads((BASE/'branching-sampler-validation.json').read_text());assert calibration['status']=='pass'
    for name,digest in calibration['input_sha256'].items():assert sha(BASE/name)==digest
    norm=json.loads((BASE/'channel-log-normalization.json').read_text())
    events=json.loads((BASE/'channel-log-event-integrals.json').read_text())['events']
    selected=json.loads((BASE/'channel-log-selection-integrals.json').read_text())
    channels=design['channels']
    normalizers=np.array([next(r['fraction'] for r in norm['channels'] if r['channel']==c) for c in channels])
    halves=np.array([next(r['half_fractions'] for r in norm['channels'] if r['channel']==c) for c in channels]).T
    n=np.array([r['rows'] for r in events]);raw_means=np.array([r['mean'] for r in events]);raw_outer=np.array([r['outer'] for r in events])
    assert len(events)==69 and n.sum()==11180578
    hashes={p.name:sha(p) for p in [Path(__file__),BASE/'branching_posterior.py',BASE/'branching-sampler-design.json',BASE/'branching-sampler-validation.json',BASE/'likelihood-design.json',BASE/'channel-log-integrals-complete.json',BASE/'channel-weight-validation.json',BASE/'channel-product-validation.json']}
    protocol=hashlib.sha256(json.dumps(hashes,sort_keys=True).encode()).hexdigest()
    path=BASE/'branching-results.json'
    report=json.loads(path.read_text()) if path.exists() else dict(protocol=protocol,status='partial_exploratory_fits',cases={},input_sha256=hashes)
    assert report['protocol']==protocol,'Preserve/version posterior outputs if the protocol changes.'
    for family_index,family in enumerate(design['hidden_spin_conditionals']):
        raw=raw_means[:,family_index,:]
        for bank,record in selected['banks'].items():
            for policy,selection in record['policies'].items():
                raw_efficiency=np.array(selected['early']['total'])[family_index]+np.array(selection['total'])[family_index]
                raw_selection_outer=np.array(selected['early']['outer'])[family_index]+np.array(selection['outer'])[family_index]
                for epsilon in design['regularization_epsilon']:
                    key=f'{family}:{bank}:{policy}:epsilon={epsilon:g}'
                    if key in report['cases']:continue
                    T=transform(normalizers,epsilon);means=raw@T;efficiency=raw_efficiency@T
                    assert np.all(efficiency>0)
                    reps=[];normalization_tests=[];diagnostics={}
                    for ri,seed in enumerate(sampler['replicate_seeds']):
                        result=sample(means,efficiency,sampler['pilot_power'],sampler['final_power'],seed,keep=ri==0)
                        assert result['concentration_ess']>gates['minimum_concentration_ess_per_replicate'],(key,'ess',result['concentration_ess'])
                        if ri==0:
                            beta=result.pop('samples');weights=result.pop('weights');result.pop('log_sampling_density')
                            base_log=log_likelihood(beta,means,efficiency)
                            for half in halves:
                                HT=transform(half,epsilon);hm=raw@HT;hx=raw_efficiency@HT
                                logw=np.log(weights)+log_likelihood(beta,hm,hx)-base_log
                                hw=np.exp(logw-logsumexp(logw))
                                normalization_tests.append(dict(intrinsic_quantiles=weighted_quantiles(beta,hw).tolist(),concentration_ess=float(1/np.sum(hw*hw))))
                            center=np.array(result['intrinsic_mean'])
                            mixed=means@center
                            second=np.einsum('c,ncd,d->n',T@center,raw_outer[:,family_index],T@center)/n
                            event_ess=n*mixed*mixed/second
                            selection_second=(T@center)@raw_selection_outer@(T@center)
                            diagnostics=dict(minimum_event_weight_concentration=float(event_ess.min()),
                                             minimum_event_index=int(events[int(np.argmin(event_ess))]['index']),
                                             selection_weight_concentration=float((efficiency@center)**2/selection_second),
                                             scope='Weight concentration at posterior mean; not PE-chain/campaign effective size or coverage calibration.')
                        reps.append(result)
                    quantiles=np.array([r['intrinsic_quantiles'] for r in reps]);detected=np.array([r['detected_quantiles'] for r in reps])
                    span=float(np.ptp(quantiles,axis=0).max())
                    assert span<gates['replicate_maximum_quantile_span'],(key,'replicate_span',span)
                    half_shift=float(np.max(abs(np.array([r['intrinsic_quantiles'] for r in normalization_tests])-quantiles[0])))
                    report['cases'][key]=dict(family=family,bank=bank,policy=policy,epsilon=epsilon,
                        intrinsic_quantiles=np.mean(quantiles,axis=0).tolist(),detected_quantiles=np.mean(detected,axis=0).tolist(),
                        replicate_quantile_span=span,replicates=reps,normalization_half_tests=normalization_tests,
                        maximum_normalizer_half_quantile_shift=half_shift,concentration=diagnostics,
                        efficiency=efficiency.tolist(),zero_event_components=[dict(index=events[i]['index'],channel=channels[j]) for i,j in np.argwhere(means==0)])
                    write(path,report)
                    print(json.dumps(dict(case=key,median=np.mean(quantiles,axis=0)[:,1].tolist(),replicate_span=span)),flush=True)
    assert len(report['cases'])==24
    report.update(status='exploratory_conditional_fits_computed_full_product_validation_pending',events=69,rows=int(n.sum()),channels=channels,
                  context=design['context'],functional=design['branching_parameter'],limits=design['limits'])
    write(path,report)


if __name__=='__main__':
    try:main()
    except Exception as error:
        path=BASE/'branching-fit-failures.json';h=json.loads(path.read_text()) if path.exists() else []
        h.append(dict(tool_sha256=sha(Path(__file__)),error=repr(error)))
        write(path,h);raise
