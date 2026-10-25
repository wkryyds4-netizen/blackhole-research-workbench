"""Exact tutorial -> cosmo -> original-prior row provenance for four new events."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import numpy as np
import h5py

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
PARAMS=['mass_1','mass_ratio','redshift','a_1','a_2','cos_tilt_1','cos_tilt_2']
SOURCE=['mass_1_source']+PARAMS[1:]
EXTRA=['luminosity_distance','mass_1','mass_2','chirp_mass','log_likelihood']
PRIORS=['chirp_mass','mass_ratio','mass_1','mass_2','luminosity_distance','a_1','a_2','tilt_1','tilt_2']


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def read(group,path):
    if path not in group or not isinstance(group[path],h5py.Dataset):return None
    return [x.decode() if isinstance(x,bytes) else x.item() if isinstance(x,np.generic) else x
            for x in np.atleast_1d(group[path][()])]


def content_hash(group):
    result={}
    def visit(name,item):
        if isinstance(item,h5py.Dataset):result[name]=repr(item[()].tolist())
    group.visititems(visit)
    return hashlib.sha256(json.dumps(result,sort_keys=True).encode()).hexdigest()


def join(left,right):
    """Exact seven-dimensional matching, returning -1 for absent rows."""
    l=np.rec.fromarrays(left);r=np.rec.fromarrays(right)
    unique,first,counts=np.unique(r,return_index=True,return_counts=True)
    pos=np.searchsorted(unique,l);valid=pos<len(unique)
    safe=np.minimum(pos,len(unique)-1);valid &= unique[safe]==l
    rows=np.where(valid,first[safe],-1)
    return rows,dict(matches=int(valid.sum()),ambiguous=int(np.sum(valid&(counts[safe]>1))),
                     left_rows=len(l),right_rows=len(r),unique_right_rows=len(unique))


def prior_audit(group,samples):
    strings={key:read(group,'priors/analytic/'+key) for key in PRIORS}
    expected={'chirp_mass':'UniformInComponentsChirpMass(', 'mass_ratio':'UniformInComponentsMassRatio(',
              'luminosity_distance':'PowerLaw(alpha=2,','a_1':'Uniform(','a_2':'Uniform(',
              'tilt_1':'Sine(','tilt_2':'Sine('}
    shapes={key:bool(strings[key] and strings[key][0].startswith(value)) for key,value in expected.items()}
    bounds={};violations={}
    for key,values in strings.items():
        if not values:continue
        limits={name:float(value) for name,value in re.findall(r'(minimum|maximum)=([\d.eE+\-]+)',values[0])}
        if len(limits)!=2:continue
        bounds[key]=limits
        vals=samples[key] if key in samples.dtype.names else np.arccos(samples['cos_'+key])
        violations[key]=int(np.sum((vals<limits['minimum']-1e-10)|(vals>limits['maximum']+1e-10)))
    shapes['isotropic_full_tilt_support']=all(
        key in bounds and bounds[key]['minimum']==0 and abs(bounds[key]['maximum']-np.pi)<1e-12
        for key in ['tilt_1','tilt_2'])
    return dict(analytic_strings=strings,shape_checks=shapes,bounds=bounds,
                retained_row_bound_violations=violations,
                support_scope='Analytic file declarations and retained rows only; not proof of absent likelihood tails or all historic constraints.')


def extract(index,downloads,archive):
    files=[x for x in downloads if x['index']==index]
    assert len(files)==2,'download both representations first'
    co=next(x for x in files if x['file'].endswith('_cosmo.h5'))
    no=next(x for x in files if x['file'].endswith('_nocosmo.h5'))
    for x in files:assert digest(BASE/'source-data'/x['file'])==x['sha256']
    local=[archive[f'{index}_{k}'] for k in PARAMS]
    label_results=[];maps={}
    with h5py.File(BASE/'source-data'/co['file'],'r') as hc,h5py.File(BASE/'source-data'/no['file'],'r') as hn:
        for label in hc:
            if not isinstance(hc[label],h5py.Group) or 'posterior_samples' not in hc[label]:continue
            samples=hc[label]['posterior_samples'][:]
            if not set(SOURCE)<=set(samples.dtype.names or []):continue
            mapping,info=join(local,[samples[k] for k in SOURCE])
            label_results.append(dict(label=label,**info));maps[label]=mapping
        full=[r for r in label_results if r['matches']==len(local[0]) and r['ambiguous']==0]
        assert len(full)==1, f'ambiguous or unmatched tutorial label: {label_results}'
        label=full[0]['label'];assert label=='C01:IMRPhenomXPHM',f'Prior recovery required for {label}'
        c=hc[label]['posterior_samples'][:];n=hn[label]['posterior_samples'][:]
        cn,info=join([c[k] for k in SOURCE],[n[k] for k in SOURCE])
        assert info['matches']==len(c) and info['ambiguous']==0
        shared={key:bool(np.array_equal(c[key],n[key][cn])) for key in SOURCE+EXTRA}
        assert all(shared.values()),shared
        lc=maps[label];ln=cn[lc]
        assert len(np.unique(lc))==len(lc) and len(np.unique(cn))==len(cn), 'Repeated source rows require weight-aware diagnostics'
        assert all(np.array_equal(local[j],n[key][ln]) for j,key in enumerate(SOURCE))
        meta={key:read(hn[label],key) for key in ['approximant','meta_data/meta_data/f_ref','meta_data/meta_data/cosmology']}
        assert all(read(hc[label],key)==value for key,value in meta.items())
        config_ref=read(hn[label],'config_file/config/reference-frequency')
        assert config_ref and float(config_ref[0])==float(meta['meta_data/meta_data/f_ref'][0])
        config_same=content_hash(hn[label]['config_file'])==content_hash(hc[label]['config_file'])
        prior_same=content_hash(hn[label]['priors/analytic'])==content_hash(hc[label]['priors/analytic'])
        assert config_same and prior_same
        audit=prior_audit(hn[label],n)
        assert all(audit['shape_checks'].values()),audit['shape_checks']
        assert not any(audit['retained_row_bound_violations'].values()),audit['retained_row_bound_violations']
        assert meta['meta_data/meta_data/cosmology']==['Planck15_lal']
        bank=dict(tutorial_to_cosmo=lc,tutorial_to_nocosmo=ln,cosmo_to_nocosmo=cn)
        for prefix,samples in [('cosmo',c),('nocosmo',n)]:
            for key in SOURCE+EXTRA:
                output='mass_1' if key=='mass_1_source' else 'mass_1_detector' if key=='mass_1' else key
                bank[f'{prefix}_{output}']=samples[key]
        np.savez_compressed(BASE/f'event-{index}.npz',**bank)
        result=dict(index=index,event=co['event'],analysis_label=label,labels=label_results,
            tutorial_rows=len(lc),cosmo_rows=len(c),nocosmo_rows=len(n),cosmo_to_nocosmo=info,
            unique_tutorial_source_rows=len(np.unique(ln)),unique_cosmo_source_rows=len(np.unique(cn)),
            config_reference_frequency=config_ref,
            exact_shared_fields=shared,all_tutorial_rows_exact_in_both=True,config_equal=config_same,
            prior_strings_equal=prior_same,metadata=meta,prior_audit=audit,sources=dict(cosmo=co,nocosmo=no),
            cache_sha256=digest(BASE/f'event-{index}.npz'))
    print(json.dumps({k:result[k] for k in ['index','event','analysis_label','tutorial_rows','cosmo_rows','nocosmo_rows','metadata']},indent=2),flush=True)
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--indices',type=int,nargs='+',default=[19,20,25,50]);args=parser.parse_args()
    downloads=json.loads((BASE/'downloads.json').read_text())
    expected=json.loads((ROOT/'research/population-stage3/validation.json').read_text())['frozen_hashes']['data/posteriors.npz']
    assert digest(ROOT/'data/posteriors.npz')==expected
    path=BASE/'event-provenance.json';results=json.loads(path.read_text()) if path.exists() else []
    with np.load(ROOT/'data/posteriors.npz',allow_pickle=False) as archive:
        for index in args.indices:
            result=extract(index,downloads,archive)
            results=[x for x in results if x['index']!=index]+[result]
            path.write_text(json.dumps(results,indent=2),encoding='utf-8')


if __name__=='__main__':main()
