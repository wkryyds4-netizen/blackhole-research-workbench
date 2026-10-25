"""O3 search-union control using LVK GWTC-3 Appendix A's matched-filter rule.

The 69 candidate identifiers match the paper's BBH list (stage three).
The paper states that observed cWB detections were corroborated by the other
pipelines. This is documentary support; individual search logs are not recovered.
Only the O3 detection union is changed. No new data-driven event deletion.
"""
from pathlib import Path
import json
import hashlib
import numpy as np
from run_selection_study import init_model,physical,weights,posterior,KEYS

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]


def main():
    init_model()
    provenance=json.loads((BASE/'o3-provenance.json').read_text())
    previous=json.loads((BASE/'selection-results.json').read_text())
    membership=json.loads((ROOT/'research/population-stage3/event-membership.json').read_text())
    assert membership['exact_identifier_set_match'] and membership['local_selected_count']==69
    assert hashlib.sha256((BASE/'o3-bank.npz').read_bytes()).hexdigest()==provenance['cache_sha256']
    with np.load(BASE/'o3-bank.npz',allow_pickle=False) as bank:
        full={k[5:]:bank[k] for k in bank.files if k.startswith('full_')}
        old={k[9:]:bank[k] for k in bank.files if k.startswith('combined_') and k!='combined_o3_to_full'}
    names=sorted(k for k in full if k.startswith('ifar_'))
    matched=[k for k in names if k!='ifar_cwb']
    assert len(names)==5 and len(matched)==4
    mask=np.any([full[k]>1 for k in names],axis=0)
    cwb=full['ifar_cwb']>1;mf=np.any([full[k]>1 for k in matched],axis=0)
    assert np.array_equal(mask,mf|cwb)
    raw={k:v[mask] for k,v in full.items()};w=weights(physical(raw))
    specific=(cwb&~mf)[mask];assert np.all(mf<=mask)
    factor=previous['o3_time_fraction']/provenance['attributes']['total_generated']
    o3_all=w.sum(axis=0)*factor;o3_mf=w[~specific].sum(axis=0)*factor;removed=w[specific].sum(axis=0)*factor
    early=np.array(previous['selection_components']['o1'])+np.array(previous['selection_components']['o2'])
    all_selection=early+o3_all;mf_selection=early+o3_mf
    error=float(np.max(abs(all_selection/np.array(previous['selection_components']['full_o3'])-1)))
    partition=float(np.max(abs((o3_mf+removed)/o3_all-1)))
    assert error<1e-11 and partition<1e-11
    with np.load(BASE/'selection-checkpoints.npz',allow_pickle=False) as cp:
        means=cp['means'][cp['names'].tolist().index('five_full_o3')]
        r=cp['r'];fs=cp['fs']
    qa,pa=posterior(means,all_selection);qm,pm=posterior(means,mf_selection)
    qfine,_=posterior(means,mf_selection,2001)
    assert np.max(abs(qa-np.array(previous['scenarios']['five_full_o3']['quantiles'])))<1e-11
    assert np.max(abs(qm-qfine))<1e-6
    diagnostics=[]
    for k,key in enumerate(KEYS):
        full_curve=fs[k]*all_selection[k,0]+(1-fs[k])*all_selection[k,1]
        drop_curve=fs[k]*removed[k,0]+(1-fs[k])*removed[k,1]
        fraction=drop_curve/full_curve
        diagnostics.append(dict(key=key,removed_fraction_of_total_selection_range=[float(fraction.min()),float(fraction.max())],
            removed_fraction_at_union_median=float(np.interp(qa[k,1],r,fraction)),
            median_shift_pp=float(100*(qm[k,1]-qa[k,1]))))
    old_o3=old['name']==b'o3';old_mf=np.any([old[k]>1 for k in matched],axis=0);old_cwb=old['ifar_cwb']>1
    result=dict(schema='blackhole-o3-pipeline-policy/1',keys=KEYS,
        policy_source='https://arxiv.org/html/2111.03634v3#A1',
        policy_description='Appendix A omits cWB from BBH sensitivity because observed binary detections were corroborated by other pipelines. IFAR>1 year union of GstLAL, MBTA, PyCBC-BBH and PyCBC-hyperbank.',
        event_support='All 69 proposed identifiers match the paper BBH set; policy compatibility relies on the paper statement, not a new independent audit of each search FAR.',
        counts=dict(full_union=int(mask.sum()),full_matched_filter=int(mf.sum()),full_cwb=int(cwb.sum()),full_cwb_only=int(specific.sum()),
            combined_o3_union=int(((old_mf|old_cwb)&old_o3).sum()),combined_o3_matched_filter=int((old_mf&old_o3).sum()),combined_o3_cwb_only=int((old_cwb&~old_mf&old_o3).sum())),
        diagnostics=diagnostics,union_quantiles=qa.tolist(),matched_filter_quantiles=qm.tolist(),
        union_width_gaps_pp=[float(100*(qa[0,1]-qa[2,1])),float(100*(qa[3,1]-qa[5,1]))],
        matched_filter_width_gaps_pp=[float(100*(qm[0,1]-qm[2,1])),float(100*(qm[3,1]-qm[5,1]))],
        selection_components=dict(union=all_selection.tolist(),matched_filter=mf_selection.tolist(),removed=removed.tolist()),
        validation=dict(original_selection_relative_error=error,partition_relative_error=partition,max_grid_quantile_error=float(np.max(abs(qm-qfine)))),
        inputs=dict(o3_bank=provenance['cache_sha256'],selection_results=hashlib.sha256((BASE/'selection-results.json').read_bytes()).hexdigest(),membership=hashlib.sha256((ROOT/'research/population-stage3/event-membership.json').read_bytes()).hexdigest()),
        limits='Conditional fixed-69-event six-model comparison. Retains unresolved native spin references, prior support, O1/O2 semianalytic sensitivity, O3 proposal/tail and finite-MC limitations. Not a full LVK reproduction or formation-channel inference.')
    (BASE/'pipeline-results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    np.savez_compressed(BASE/'pipeline-checkpoints.npz',r=r,fs=fs,means=means,union=pa,matched_filter=pm,selection_union=all_selection,selection_matched_filter=mf_selection)
    print(json.dumps({k:result[k] for k in ['counts','diagnostics','union_width_gaps_pp','matched_filter_width_gaps_pp','validation']},indent=2))


if __name__=='__main__':main()
