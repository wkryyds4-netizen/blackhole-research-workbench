"""Resolve index48 using the independently matched cosmological release.

The nocosmo Mixed label is not equivalent to Phenom. Every tutorial row is
instead traced through the unique cosmo Phenom analysis and then injectively
back to the original-prior Phenom analysis with unchanged native inputs.
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
from extract_events import join,content_hash,read


def digest(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    old=json.loads((ROOT/'research/population-stage3/pilot-downloads.json').read_text())
    co=next(r for r in old if r['event']=='GW191103_012549')
    no=next(r for r in json.loads((BASE/'download-results.json').read_text()) if r['index']==48)
    copath=ROOT/'research/population-stage3/source-data'/co['file']
    nopath=ROOT/no.get('local_file','research/population-catalog/source-data/'+no['source']['file'])
    assert digest(copath)==co['sha256'] and digest(nopath)==no['sha256']
    with np.load(ROOT/'data/posteriors.npz',allow_pickle=False) as archive:
        left=[archive[f'48_{k}'] for k in ['mass_1']+SOURCE[1:]]
    canonical='C01:IMRPhenomXPHM';co_matches=[]
    with h5py.File(copath,'r') as hc,h5py.File(nopath,'r') as hn:
        for label,group in hc.items():
            if not isinstance(group,h5py.Group) or 'posterior_samples' not in group:continue
            ds=group['posterior_samples']
            if not isinstance(ds,h5py.Dataset) or not set(SOURCE)<=set(ds.dtype.names or []):continue
            data=ds.fields(SOURCE)[:];_,info=join(left,[data[k] for k in SOURCE]);co_matches.append(dict(label=label,**info))
        full=[r['label'] for r in co_matches if r['matches']==len(left[0]) and r['ambiguous']==0]
        assert full==[canonical],co_matches
        a=hc[canonical+'/posterior_samples'].fields(NATIVE)[:]
        b=hn[canonical+'/posterior_samples'].fields(NATIVE)[:]
        mapping,info=join([a[k] for k in SOURCE],[b[k] for k in SOURCE])
        assert info['matches']==len(a) and info['ambiguous']==0 and len(np.unique(mapping))==len(a)
        equality={k:bool(np.array_equal(a[k],b[k][mapping],equal_nan=True)) for k in NATIVE};assert all(equality.values()),equality
        tutorial,identity=join(left,[a[k] for k in SOURCE]);assert identity['matches']==len(left[0]) and identity['ambiguous']==0
        original=mapping[tutorial];assert len(np.unique(original))==len(left[0])
        direct,direct_info=join(left,[b[k] for k in SOURCE]);assert np.array_equal(direct,original)
        config_equal=content_hash(hc[canonical+'/config_file'])==content_hash(hn[canonical+'/config_file'])
        prior_equal=content_hash(hc[canonical+'/priors/analytic'])==content_hash(hn[canonical+'/priors/analytic'])
        assert config_equal and prior_equal
        metadata={k:read(hn[canonical],k) for k in ['meta_data/meta_data/f_ref','meta_data/meta_data/cosmology','approximant']}
        assert all(read(hc[canonical],k)==v for k,v in metadata.items())
        assert metadata['approximant']==['IMRPhenomXPHM']
    report=dict(schema='blackhole-cosmo-native-lineage-resolution/1',index=48,event=no['event'],canonical=canonical,
        cosmo_source=dict(file=co['file'],sha256=co['sha256']),nocosmo_source=dict(file=no['source']['file'],sha256=no['sha256']),
        cosmo_label_matches=co_matches,unique_cosmo_tutorial_label=canonical,cosmo_to_original=info,
        native_columns_exact=equality,tutorial_to_cosmo=identity,tutorial_to_original=direct_info,
        three_way_tutorial_rows_exact=True,config_equal=config_equal,original_analytic_prior_declarations_equal=prior_equal,metadata=metadata,
        nocosmo_container_label_unique=False,mixed_full_posterior_equivalent=False,
        resolution='Unique matched cosmo Phenom analysis, injective all-native-column/config/prior identity into full original-prior Phenom analysis. Multiple nocosmo containing labels remain.',
        limits='Recovers the single-waveform analysis lineage. Does not assert that the nocosmo Mixed posterior is equivalent or uniquely identify every intermediate preprocessing command.')
    (BASE/'cosmo-lineage-resolution.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['event','unique_cosmo_tutorial_label','cosmo_to_original','three_way_tutorial_rows_exact','config_equal','mixed_full_posterior_equivalent']},indent=2))


if __name__=='__main__':main()
