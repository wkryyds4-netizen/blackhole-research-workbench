"""Identify a complete single-waveform component inside index72 Mixed.

This resolves the waveform-origin analysis of all tutorial points while
retaining ambiguity about containing labels and preprocessing commands.
Mixed itself is never declared equivalent to the selected Phenom analysis.
"""
from pathlib import Path
import hashlib
import json
import sys
import h5py
import numpy as np
from diagnose_alias import NATIVE,SOURCE

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
sys.path.insert(0,str(ROOT/'research/population-stage7'))
from extract_events import join


def main():
    record=next(r for r in json.loads((BASE/'download-results.json').read_text()) if r['index']==72)
    path=ROOT/record.get('local_file','research/population-catalog/source-data/'+record['source']['file'])
    with path.open('rb') as stream:assert hashlib.file_digest(stream,'sha256').hexdigest()==record['sha256']
    canonical='C01:IMRPhenomXPHM';other='C01:SEOBNRv4PHM';mixed='C01:Mixed'
    with h5py.File(path,'r') as h:
        a=h[canonical+'/posterior_samples'].fields(NATIVE)[:]
        b=h[other+'/posterior_samples'].fields(NATIVE)[:]
        m=h[mixed+'/posterior_samples'].fields(NATIVE)[:]
    pm,pinfo=join([a[k] for k in SOURCE],[m[k] for k in SOURCE])
    assert pinfo['matches']==len(a) and pinfo['ambiguous']==0 and len(np.unique(pm))==len(a)
    pequal={k:bool(np.array_equal(a[k],m[k][pm],equal_nan=True)) for k in NATIVE};assert all(pequal.values())
    remainder=np.ones(len(m),dtype=bool);remainder[pm]=False
    eb,einfo=join([m[k][remainder] for k in SOURCE],[b[k] for k in SOURCE])
    assert einfo['matches']==remainder.sum() and einfo['ambiguous']==0 and len(np.unique(eb))==len(eb)
    eequal={k:bool(np.array_equal(m[k][remainder],b[k][eb],equal_nan=True)) for k in NATIVE};assert all(eequal.values())
    # Show the two generating analyses have no overlapping seven-coordinate points.
    _,overlap=join([a[k] for k in SOURCE],[b[k] for k in SOURCE]);assert overlap['matches']==0
    with np.load(ROOT/'data/posteriors.npz',allow_pickle=False) as archive:
        tutorial=[archive[f'72_{k}'] for k in ['mass_1']+SOURCE[1:]]
    ta,ainfo=join(tutorial,[a[k] for k in SOURCE]);tb,binfo=join(tutorial,[b[k] for k in SOURCE]);tm,minfo=join(tutorial,[m[k] for k in SOURCE])
    assert ainfo['matches']==len(tutorial[0]) and ainfo['ambiguous']==0 and binfo['matches']==0
    assert np.array_equal(pm[ta],tm) and minfo['ambiguous']==0
    report=dict(schema='blackhole-complete-waveform-component-resolution/1',index=72,event=record['event'],source_sha256=record['sha256'],canonical=canonical,
        canonical_rows=len(a),mixed_rows=len(m),other_analysis_rows=len(b),
        canonical_to_mixed=pinfo,mixed_remainder_to_other=einfo,canonical_native_columns_exact=pequal,other_native_columns_exact=eequal,
        disjoint_analysis_input_points=overlap['matches']==0,complete_mixed_partition=True,
        tutorial_to_canonical=ainfo,tutorial_to_other=binfo,tutorial_to_mixed=minfo,tutorial_origin_waveform_unique=True,
        containing_label_unique=False,mixed_full_posterior_equivalent=False,
        resolution='All tutorial points belong uniquely to the Phenom waveform-origin analysis. Its entire native posterior is an exact component of Mixed; the disjoint remainder maps to SEOB. Select the full Phenom analysis for the fixed single-waveform control.',
        limits='Containing-label identity and exact historical preprocessing remain ambiguous. This is not equivalence of the full Mixed and Phenom posteriors or an LVK waveform-marginalized population reproduction.')
    (BASE/'mixture-component-resolution.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['event','canonical_rows','mixed_rows','other_analysis_rows','complete_mixed_partition','tutorial_origin_waveform_unique','containing_label_unique']},indent=2))


if __name__=='__main__':main()
