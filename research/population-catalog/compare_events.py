"""Conditional full-posterior expansion on a fixed selection and common r prior."""
from pathlib import Path
import argparse
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
    parser=argparse.ArgumentParser();parser.add_argument('--indices',nargs='+',type=int,default=[0,14]);args=parser.parse_args();assert set(args.indices)<={0,14}
    init_model();frozen=json.loads((ROOT/'results/model.json').read_text());provenance=json.loads((BASE/'event-provenance.json').read_text())
    with np.load(ROOT/'research/population-selection-audit/pipeline-checkpoints.npz',allow_pickle=False) as cp:
        means=cp['means'];selection=cp['selection_matched_filter']
    r,fs=coordinates();qref,_=posterior(means,selection);replacements={};details=[]
    with np.load(ROOT/'data/posteriors.npz',allow_pickle=False) as archive:
        for index in args.indices:
            info=next(v for v in provenance if v['index']==index);assert info.get('native_source_identity_resolved',info['uniquely_verified_tutorial_source'])
            assert all(info['prior_audit']['shape_checks'].values()) and not any(info['prior_audit']['retained_row_bound_violations'].values())
            path=BASE/f'event-{index}.npz';assert hashlib.sha256(path.read_bytes()).hexdigest()==info['cache_sha256']
            with np.load(path,allow_pickle=False) as bank:
                data={k:bank[k] for k in PARAMS+['luminosity_distance']}
            pc,pn,dl=density_pair(data);w=weights(dict(data,prior=pn))
            after=w.mean(axis=0);identity=(weights(dict(data,prior=pc))*(pc/pn)[:,None,None]).mean(axis=0)
            identity_error=float(np.max(abs(identity/after-1)));assert identity_error<1e-11
            dl_error=float(np.max(abs(dl/data['luminosity_distance']-1)));assert dl_error<1e-3
            pos=next(i for i,e in enumerate(frozen['events']) if e['index']==index)
            local={k:archive[f'{index}_{k}'] for k in PARAMS+['prior']};before=weights(local).mean(axis=0)
            frozen_moments=np.array([frozen['scenarios'][k]['event_moments'][pos][:2] for k in KEYS])
            moment_error=float(np.max(abs(before/frozen_moments-1)));assert moment_error<1e-11
            assert np.max(abs(before/means[pos]-1))<1e-11
            replacements[index]=(pos,after);diagnostics=[]
            for k,key in enumerate(KEYS):
                f=np.interp(qref[k,1],r,fs[k]);one=f*w[:,k,0]+(1-f)*w[:,k,1]
                diagnostics.append(dict(key=key,ess_at_previous_median=float(one.sum()**2/(one@one)),largest_weight_fraction=float(one.max()/one.sum())))
            details.append(dict(index=index,event=info['event'],rows=len(w),moment_relative_error=moment_error,
                density_identity_relative_error=identity_error,distance_relative_error=dl_error,weight_diagnostics=diagnostics,
                original_mean=before.tolist(),new_mean=after.tolist()))
    cases={'five_event_baseline':means.copy()}
    for index,(pos,after) in replacements.items():
        change=means.copy();change[pos]=after;cases[f'plus_{index}']=change
    joint=means.copy()
    for pos,after in replacements.values():joint[pos]=after
    cases['all_new_joint']=joint
    scenarios={};curves=[];qs=[];maxgrid=0.
    for name,em in cases.items():
        q,p=posterior(em,selection);fine,_=posterior(em,selection,2001);maxgrid=max(maxgrid,float(np.max(abs(q-fine))))
        scenarios[name]=dict(quantiles=q.tolist(),median_shift_pp=(100*(q[:,1]-qref[:,1])).tolist(),
            width_gaps_pp=[float(100*(q[0,1]-q[2,1])),float(100*(q[3,1]-q[5,1]))]);curves.append(p);qs.append(q)
    assert maxgrid<1e-6
    result=dict(schema='blackhole-catalog-event-expansion/1',indices=args.indices,keys=KEYS,scenarios=scenarios,event_diagnostics=details,max_grid_quantile_error=maxgrid,
        inputs=dict(frozen_model_sha256=hashlib.sha256((ROOT/'results/model.json').read_bytes()).hexdigest(),
            pipeline_checkpoint_sha256=hashlib.sha256((ROOT/'research/population-selection-audit/pipeline-checkpoints.npz').read_bytes()).hexdigest(),
            provenance_sha256=hashlib.sha256((BASE/'event-provenance.json').read_bytes()).hexdigest()),
        limits='Fixed 69 events and historical O1/O2 plus full O3 matched-filter selection. Added original-prior full posteriors remain in native reference coordinates. No joint reference transport, full catalog restoration, or astrophysical channel fractions.')
    prefix='catalog' if args.indices==[0,14] else 'partial-'+'-'.join(map(str,args.indices))
    (BASE/f'{prefix}-comparison.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    np.savez_compressed(BASE/f'{prefix}-comparison.npz',r=r,fs=fs,selection=selection,names=np.array(list(cases)),means=np.array(list(cases.values())),quantiles=np.array(qs),densities=np.array(curves),zgrid=model.ZGRID,zdv=model.ZDV,znorm=model.ZNORM[2])
    print(json.dumps(dict(scenarios=scenarios,max_grid_quantile_error=maxgrid),indent=2))


if __name__=='__main__':main()
