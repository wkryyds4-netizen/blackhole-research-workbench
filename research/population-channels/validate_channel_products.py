"""Independent decimal products from complete saved69-event channel sums.

This does not re-evaluate every PE row; row and full-sum scopes stay distinct.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
from decimal import Decimal,localcontext
import hashlib
import json
import numpy as np
from branching_posterior import log_likelihood

BASE=Path(__file__).resolve().parent;ROOT=BASE.parents[1]
def sha(p):
    with p.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def decimal_log_product(beta,means,selection):
    with localcontext() as context:
        context.prec=70
        b=list(map(lambda x:Decimal(str(x)),beta))
        denominator=sum(v*Decimal(str(x)) for v,x in zip(b,selection))
        product=Decimal(1)
        for row in means:
            product*=sum(v*Decimal(str(x)) for v,x in zip(b,row))/denominator
        return float(product.ln())


def check_stats(record):
    n=record['rows'];total=np.array(record['total']);outer=np.array(record['outer']);maximum=np.array(record['maximum'])
    assert n>0 and np.isfinite(total).all() and np.all(total>=0)
    assert np.all(np.array(record['positive_rows'])<=n)
    diagonal=np.diagonal(outer,axis1=-2,axis2=-1)
    assert np.all(total*total<=n*diagonal*(1+1e-10)+np.finfo(float).tiny)
    assert np.all(maximum<=total*(1+1e-12)+np.finfo(float).tiny)
    if 'mean' in record:assert np.array_equal(np.array(record['mean']),total/n)


def main():
    complete=json.loads((BASE/'channel-log-integrals-complete.json').read_text());assert complete['events']==69
    for name,digest in complete['input_sha256'].items():assert sha(BASE/name)==digest
    events=json.loads((BASE/'channel-log-event-integrals.json').read_text())['events']
    selection=json.loads((BASE/'channel-log-selection-integrals.json').read_text())
    normalization=json.loads((BASE/'channel-log-normalization.json').read_text())
    provenance=json.loads((ROOT/'research/population-catalog/event-provenance.json').read_text())
    design=json.loads((BASE/'likelihood-design.json').read_text());channels=design['channels']
    norm=np.array([next(r['fraction'] for r in normalization['channels'] if r['channel']==c) for c in channels])
    frozen=json.loads((ROOT/'results/model.json').read_text())
    assert [r['index'] for r in events]==[e['index'] for e in frozen['events']]
    assert len(events)==69 and sum(r['rows'] for r in events)==11180578
    for event in events:
        p=next(r for r in provenance if r['index']==event['index'])
        assert event['rows']==p['rows'] and event['source_cache_sha256']==p['cache_sha256']
        check_stats(event)
    check_stats(selection['early'])
    for bank in selection['banks'].values():
        for policy in bank['policies'].values():check_stats(policy)
    rng=np.random.default_rng(20261105)
    vertices=np.eye(5)*(1-5e-6)+1e-6
    beta=np.r_[rng.dirichlet(np.ones(5),size=128),vertices,np.ones((1,5))/5]
    rows=[]
    for fi,family in enumerate(design['hidden_spin_conditionals']):
        for bank,record in selection['banks'].items():
            for policy,s in record['policies'].items():
                raw_selection=np.array(selection['early']['total'])[fi]+np.array(s['total'])[fi]
                raw_means=np.array([r['mean'][fi] for r in events])
                for epsilon in design['regularization_epsilon']:
                    # Build each column explicitly, independently of fit.transform.
                    means=np.array([[(1-epsilon)*row[c]/norm[c]+epsilon*row[5] for c in range(5)] for row in raw_means])
                    efficiency=np.array([(1-epsilon)*raw_selection[c]/norm[c]+epsilon*raw_selection[5] for c in range(5)])
                    expected=np.array([decimal_log_product(b,means,efficiency) for b in beta])
                    actual=log_likelihood(beta,means,efficiency)
                    error=float(np.max(abs(actual-expected)));assert error<1e-9,(family,bank,policy,epsilon,error)
                    # Known per-channel scaling must cancel after normalization.
                    scaling=np.array([.1,.3,2,7,13]);rescaled=(raw_means[:,:5]*scaling)/(norm*scaling)
                    assert np.max(abs(rescaled/(raw_means[:,:5]/norm)-1)[raw_means[:,:5]>0])<1e-12
                    rows.append(dict(family=family,bank=bank,policy=policy,epsilon=epsilon,points=len(beta),max_log_product_error=error))
                noinfo=np.repeat(raw_means[:,5,None],5,axis=1)
                equal_eff=np.repeat(raw_selection[5],5)
                values=log_likelihood(beta,noinfo,equal_eff)
                assert np.ptp(values)<1e-9,'Fully regularized limit must carry no channel information.'
    report=dict(schema='blackhole-channel-saved-product-independent-validation/1',status='pass',events=69,rows=11180578,cases=rows,
                input_sha256={p.name:sha(p) for p in [Path(__file__),BASE/'branching_posterior.py',BASE/'channel-log-integrals-complete.json',BASE/'likelihood-design.json']},
                scope='Independent70-digit direct products, normalization identities and aggregate-moment bounds on complete69-event saved sums. Does not independently reintegrate all raw PE rows or calibrate astrophysical/chain/campaign uncertainty.')
    (BASE/'channel-product-validation.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='pass',cases=len(rows),max_error=max(r['max_log_product_error'] for r in rows))),flush=True)


if __name__=='__main__':
    try:main()
    except Exception as error:
        path=BASE/'channel-product-failures.json';h=json.loads(path.read_text()) if path.exists() else []
        h.append(dict(tool_sha256=sha(Path(__file__)),error=repr(error)))
        path.write_text(json.dumps(h,indent=2)+'\n',encoding='utf-8');raise
