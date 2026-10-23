"""Independent synthetic API audit for population_mc_bootstrap.py.

Small arrays make direct G-by-G covariance construction practical, in contrast
to the optimized implementation's two-component Jacobian reduction.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import population_mc_bootstrap as impl

ROOT=Path(__file__).resolve().parents[1]


def scalar_curves(events,selection,fs):
    b,e,k,_=events.shape
    out=np.empty((b,k,fs.shape[1]))
    for bi in range(b):
        for ki in range(k):
            for gi,f in enumerate(fs[ki]):
                out[bi,ki,gi]=sum(np.log(f*events[bi,ei,ki,0]+(1-f)*events[bi,ei,ki,1]) for ei in range(e))-e*np.log(f*selection[bi,ki,0]+(1-f)*selection[bi,ki,1])
    return out


def run():
    rng=np.random.default_rng(120261024)
    k=6
    weights=[rng.uniform(.7,1.4,(n,k,2)) for n in [23,31,43]]
    injection=rng.uniform(.5,1.7,(57,k,2))
    total=160
    r,fs=impl.r_grid(41)
    original_events=np.stack([w.mean(axis=0) for w in weights])[None]
    original_selection=(injection.sum(axis=0)/total)[None]
    logs=impl.log_curves(original_events,original_selection,fs)
    error_curves=float(np.max(abs(logs-scalar_curves(original_events,original_selection,fs))))
    assert error_curves<1e-13
    nominal,p=impl.curve_stats(r,logs)
    pquad=p[0]*impl.quadrature_weights(r)
    result=impl.delta_shape_diagnostics(weights,injection,total,fs,pquad)
    exact_shape={'events':np.zeros(k),'injections':np.zeros(k)}
    exact_point={'events':np.zeros(k),'injections':np.zeros(k)}
    errors=[]
    for w,n,mode,factor in [(w,len(w),'events',1) for w in weights]+[(injection,total,'injections',-len(weights))]:
        for ki in range(k):
            # Explicit zeros + centered empirical vectors, independent formula.
            full=np.vstack((w[:,ki],np.zeros((n-len(w),2))))
            means=full.mean(axis=0)
            centered=full-means
            expected=centered.T@centered/n**2
            errors.append(np.max(abs(expected-impl.mean_covariance(w[:,ki],n))))
            # Materialize every grid point's per-sample log influence.
            influence=np.array([factor*(f*full[:,0]+(1-f)*full[:,1])/(f*means[0]+(1-f)*means[1]) for f in fs[ki]]).T
            influence-=influence.mean(axis=0)
            c=influence.T@influence/n**2
            pair=np.diag(c)[:,None]+np.diag(c)[None,:]-2*c
            exact_shape[mode][ki]+=pquad[ki]@pair@pquad[ki]
            exact_point[mode][ki]+=pquad[ki]@np.diag(c)
    shape_error=max(abs(result[impl.KEYS[ki]][mode]['posterior_mean_pair_difference_variance']-exact_shape[mode][ki]) for ki in range(k) for mode in exact_shape)
    point_error=max(abs(result[impl.KEYS[ki]][mode]['posterior_mean_pointwise_variance']-exact_point[mode][ki]) for ki in range(k) for mode in exact_point)
    assert shape_error<1e-13 and point_error<1e-13 and max(errors)<1e-14
    # Generated counts must preserve fixed total and allow zeros among omitted trials.
    count=impl.empirical_counts(np.random.default_rng(52),7,30,10000)
    assert count.shape==(10000,7) and np.all(count>=0) and np.all(count.sum(axis=1)<=30)
    complete=impl.empirical_counts(np.random.default_rng(52),30,30,100)
    assert np.all(complete.sum(axis=1)==30)
    counts_cov=np.cov(count,rowvar=False)
    expected_counts=np.eye(7)-np.full((7,7),1/30)
    assert np.max(abs(counts_cov-expected_counts))<.06
    # Deterministic shared draws; verify direct covariance of output curves.
    ev,se=impl.generate_means(weights,injection,total,128,441)
    ev2,se2=impl.generate_means(weights,injection,total,128,441)
    assert np.array_equal(ev,ev2) and np.array_equal(se,se2)
    qs,shapes=impl.evaluate_replicates(r,fs,ev,se,original_events,original_selection,pquad)
    direct_errors=[]
    all_logs={}
    for mode in impl.MODES:
        eval_ev=ev if mode!='injections' else np.broadcast_to(original_events,ev.shape)
        eval_se=se if mode!='events' else np.broadcast_to(original_selection,se.shape)
        curve=impl.log_curves(eval_ev,eval_se,fs)
        all_logs[mode]=curve
        for ki in range(k):
            c=np.cov(curve[:,ki],rowvar=False,ddof=1)
            direct=pquad[ki]@(np.diag(c)[:,None]+np.diag(c)[None,:]-2*c)@pquad[ki]
            direct_errors.append(abs(direct-shapes[mode][ki]))
        qdirect,_=impl.curve_stats(r,curve)
        assert np.max(abs(qdirect-qs[mode]))<1e-13
    assert max(direct_errors)<1e-11
    additive=float(np.max(abs(all_logs['both']-(all_logs['events']+all_logs['injections']-logs))))
    assert additive<1e-13
    # Identical conditional curves (same f mapping) with per-model scale factors:
    # each event and selection scaling cancels in logL; shared draws essential.
    scales=np.arange(1,k+1)[None,:,None]
    clone_weights=[np.broadcast_to(w[:,:1],w.shape)*scales for w in weights]
    clone_inj=np.broadcast_to(injection[:,:1],injection.shape)*scales
    ce,cs=impl.generate_means(clone_weights,clone_inj,total,64,442)
    samefs=np.broadcast_to(fs[0],fs.shape)
    clones=impl.log_curves(ce,cs,samefs)
    clone_error=float(np.max(abs(clones-clones[:,:1])))
    clone_q,_=impl.curve_stats(r,clones)
    assert clone_error<1e-12 and np.max(abs(clone_q-clone_q[:,:1]))<1e-13
    # Linear log density admits exact exponential CDF on an interval.
    dense=np.linspace(r[0],r[-1],8001)
    slopes=np.array([-45.,-1.,0.,1.,5.,25.])
    exponential=slopes[:,None]*(dense-dense[0])
    expq,_=impl.curve_stats(dense,exponential[None])
    expectedq=[]
    for s in slopes:
        expectedq.append(dense[0]+(np.log1p(np.array([.05,.5,.95])*np.expm1(s*np.ptp(dense)))/s if s!=0 else np.array([.05,.5,.95])*np.ptp(dense)))
    exponential_error=float(np.max(abs(expq[0]-expectedq)))
    assert exponential_error<3e-9
    report={'status':'pass','date':'2026-10-24','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'implementation_sha256':hashlib.sha256(Path(impl.__file__).read_bytes()).hexdigest(),
            'errors':{'scalar_log_curves':error_curves,'exact_mean_covariance':float(max(errors)),
                      'direct_pair_covariance_shape':float(shape_error),'pointwise_variance':float(point_error),
                      'empirical_pair_covariance_shape':float(max(direct_errors)),'mode_log_additivity':additive,
                      'shared_proportional_model_cancellation':clone_error,'exponential_cdf_quantiles':exponential_error},
            'checks':['retained and full counts','deterministic generation','exact covariance including zeros','direct grid-pair shape covariance',
                      'empirical output shape covariance','shared-model cancellation','analytic posterior quantiles'],
            'scope':'Synthetic API correctness. No verification of physical injection trial interpretation, original sampling design, or missed-tail bias.'}
    dest=ROOT/'research/population-stage2/independent-implementation-validation.json'
    dest.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':run()
