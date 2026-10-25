"""Incrementally propagate every available verified native posterior.

This preserves the fixed event set, six models, common r prior, and matched
selection. It does not silently reinterpret native tilts as infinity tilts.
"""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
sys.path.insert(0,str(ROOT/'research/population-selection-audit'))
from run_selection_study import init_model,weights,density_pair,posterior,KEYS,PARAMS,coordinates
import model


def main():
    init_model();frozen=json.loads((ROOT/'results/model.json').read_text());provenance=json.loads((BASE/'event-provenance.json').read_text())
    config_audit=json.loads((BASE/'configuration-prior-row-audit.json').read_text())
    supported={r['index'] for r in config_audit['events'] if r['source_supported_default_family'] and not any(r['bounds_violations'].values())}
    protocol=hashlib.sha256(b'fixed-native-six-model-original-distance-squared-family/1'+(ROOT/'model.py').read_bytes()+(ROOT/'research/population-stage5/run_tail_study.py').read_bytes()).hexdigest()
    cache=BASE/'native-event-integrals.json';saved=json.loads(cache.read_text()) if cache.exists() else dict(protocol=protocol,events=[])
    assert saved['protocol']==protocol,'Numerical protocol changed; version the cached integrals before reuse.'
    with np.load(ROOT/'research/population-selection-audit/pipeline-checkpoints.npz',allow_pickle=False) as cp:baseline=cp['means'];selection=cp['selection_matched_filter']
    means=baseline.copy();used=[];pending=[]
    for info in sorted(provenance,key=lambda r:r['index']):
        index=info['index'];analytic=all(info.get('prior_audit',{}).get('shape_checks',{}).values()) and bool(info.get('prior_audit'))
        if analytic:
            assert not any(info['prior_audit']['retained_row_bound_violations'].values()),(index,'Released prior bound violation')
        if not info.get('native_source_identity_resolved') or not (analytic or index in supported):
            pending.append(dict(index=index,reason='source identity or prior family not yet established'));continue
        path=BASE/f'event-{index}.npz';assert hashlib.sha256(path.read_bytes()).hexdigest()==info['cache_sha256']
        prior_mode='analytic_declarations' if analytic else 'released_configuration_and_historical_defaults'
        old=next((r for r in saved['events'] if r['index']==index and r['cache_sha256']==info['cache_sha256'] and r['prior_mode']==prior_mode),None)
        if old is None:
            with np.load(path,allow_pickle=False) as bank:d={k:bank[k] for k in PARAMS+['luminosity_distance']}
            _,prior,dl=density_pair(d);w=weights(dict(d,prior=prior));new=w.mean(axis=0)
            assert np.isfinite(w).all() and np.all(w>=0) and np.all(new>0)
            error=float(np.max(abs(dl/d['luminosity_distance']-1)));assert error<1e-3
            old=dict(index=index,event=info['event'],cache_sha256=info['cache_sha256'],prior_mode=prior_mode,rows=len(prior),mean=new.tolist(),
                distance_max_relative_error=error,min_component_ess=float(np.min(w.sum(axis=0)**2/np.sum(w*w,axis=0))))
            saved['events']=[r for r in saved['events'] if r['index']!=index]+[old]
            temporary=cache.with_suffix('.tmp');temporary.write_text(json.dumps(saved,indent=2),encoding='utf-8');temporary.replace(cache)
        pos=next(i for i,e in enumerate(frozen['events']) if e['index']==index);means[pos]=old['mean'];used.append(old)
    q0,_=posterior(baseline,selection);q,p=posterior(means,selection);fine,_=posterior(means,selection,2001)
    error=float(np.max(abs(q-fine)));assert error<1e-6
    indices=sorted(r['index'] for r in used);complete=sorted(set(indices)|{19,20,25,38,50})
    result=dict(schema='blackhole-native-catalog-expansion/1',keys=KEYS,newly_recomputed_source_indices=indices,full_source_indices_in_likelihood=complete,
        completed_full_source_events=len(complete),fixed_total_events=69,pending_provenance= pending,
        remaining_full_source_indices=sorted({e['index'] for e in frozen['events']}-set(complete)),
        quantiles=q.tolist(),baseline_quantiles=q0.tolist(),median_shift_pp=(100*(q[:,1]-q0[:,1])).tolist(),
        width_gaps_pp=[float(100*(q[0,1]-q[2,1])),float(100*(q[3,1]-q[5,1]))],max_grid_quantile_error=error,
        integral_protocol=protocol,prior_modes={mode:[r['index'] for r in used if r['prior_mode']==mode] for mode in ['analytic_declarations','released_configuration_and_historical_defaults']},
        source_rows=sum(r['rows'] for r in used),
        limits='Fixed69 native-reference conditional inference. Some event families use released configuration plus historical defaults. Common spin reference, source-support likelihood tails, injection support and campaign MC are not resolved. Not channel fractions.')
    (BASE/'native-catalog-results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    r,fs=coordinates()
    np.savez_compressed(BASE/'native-catalog-checkpoints.npz',means=means,baseline=baseline,selection=selection,quantiles=q,densities=p,r=r,fs=fs,zgrid=model.ZGRID,zdv=model.ZDV,znorm=model.ZNORM[2])
    print(json.dumps({k:result[k] for k in ['completed_full_source_events','full_source_indices_in_likelihood','prior_modes','median_shift_pp','width_gaps_pp','max_grid_quantile_error']},indent=2))


if __name__=='__main__':main()
