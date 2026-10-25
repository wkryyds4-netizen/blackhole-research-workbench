"""Trace cosmo rows to the same official nocosmo analysis and inspect RIFT gaps."""
from pathlib import Path
import hashlib
import json
import numpy as np
import h5py

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
PARAMS=['mass_1_source','mass_ratio','redshift','a_1','a_2','cos_tilt_1','cos_tilt_2']
EXTRA=['luminosity_distance','chirp_mass','mass_1','mass_2','log_likelihood']

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def read(g,key):
    if key not in g or not isinstance(g[key],h5py.Dataset):return None
    return [x.decode() if isinstance(x,bytes) else float(x) for x in np.atleast_1d(g[key][()])]

def config_content(g):
    result={}
    def visit(name,obj):
        if isinstance(obj,h5py.Dataset):result[name]=repr(obj[()].tolist())
    g.visititems(visit)
    return result

def main():
    dl=json.loads((BASE/'download.json').read_text())
    no=BASE/'source-data'/dl['file'];assert digest(no)==dl['sha256']
    co=next((ROOT/'research/population-stage3/source-data').glob('*190805*_cosmo.h5'))
    dl_co=next(x for x in json.loads((ROOT/'research/population-stage3/pilot-downloads.json').read_text()) if x['event']=='GW190805_211137')
    assert digest(co)==dl_co['sha256']
    label='C01:IMRPhenomXPHM';bank={};rift=[]
    with h5py.File(no,'r') as hn,h5py.File(co,'r') as hc:
        n=hn[label]['posterior_samples'][:];c=hc[label]['posterior_samples'][:]
        right=np.rec.fromarrays([n[k] for k in PARAMS[:1]+['redshift']])
        left=np.rec.fromarrays([c[k] for k in PARAMS[:1]+['redshift']])
        unique,first,counts=np.unique(right,return_index=True,return_counts=True)
        pos=np.searchsorted(unique,left)
        assert np.all(pos<len(unique)) and np.array_equal(unique[pos],left)
        indices=first[pos]
        checks={k:bool(np.array_equal(c[k],n[k][indices])) for k in PARAMS+EXTRA}
        assert all(checks.values())
        # Independent full seven-coordinate membership, irrespective of duplicated source rows.
        source_rows=set(zip(*(n[k].tolist() for k in PARAMS)))
        assert all(tuple(row) in source_rows for row in zip(*(c[k].tolist() for k in PARAMS)))
        config_n=config_content(hn[label]['config_file']);config_c=config_content(hc[label]['config_file'])
        assert config_n==config_c
        prior_n=config_content(hn[label]['priors/analytic']);prior_c=config_content(hc[label]['priors/analytic'])
        assert prior_n==prior_c
        meta={}
        for key in ['approximant','meta_data/meta_data/f_ref','meta_data/meta_data/cosmology']:
            assert read(hn[label],key)==read(hc[label],key)
            meta[key]=read(hn[label],key)
        for short,d in [('nocosmo',n),('cosmo',c)]:
            for k in PARAMS+EXTRA:
                out='mass_1' if k=='mass_1_source' else ('mass_1_detector' if k=='mass_1' else k)
                bank[f'{short}_{out}']=d[k]
        bank['cosmo_to_nocosmo_rows']=indices
        info=dict(schema='blackhole-population-prior-representation-match/1',
            label=label,nocosmo_n=len(n),cosmo_n=len(c),nocosmo_unique_mass_redshift_pairs=len(unique),
            ambiguous_source_pairs_for_cosmo=int(np.sum(counts[pos]>1)),
            exact_matches=checks,unique_mapped_rows=len(np.unique(indices)),
            config_equal=True,analytic_prior_strings_equal=True,metadata=meta,
            source_hashes=dict(cosmo=dl_co['sha256'],nocosmo=dl['sha256']),
            row2691_nocosmo_index=int(indices[2691]),
            scope='Same released analysis and exact shared samples; files are dependent, not independent observations.')
    for path in (ROOT/'research/population-stage3/source-data').glob('*.h5'):
        with h5py.File(path,'r') as h:
            g=h['C01:SEOBNRv4PHM'];samples=g['posterior_samples'][:]
            bound=read(g,'config_file/engine/chirpmass-max')
            rift.append(dict(file=path.name,sha256=digest(path),
                analytic_prior_entries=list(g['priors/analytic']),prior_sample_entries=list(g['priors/samples']),
                version=read(g,'version'),engine_chirpmass_max=bound,
                observed_detector_chirpmass_max=float(samples['chirp_mass'].max()),
                n_above_configured_max=int(np.sum(samples['chirp_mass']>float(bound[0]))) if bound else None,
                conclusion='Final prior support not established by these metadata; no cut applied.'))
    np.savez_compressed(BASE/'representations.npz',**bank)
    (BASE/'representation-match.json').write_text(json.dumps(info,indent=2),encoding='utf-8')
    (BASE/'rift-prior-audit.json').write_text(json.dumps(rift,indent=2),encoding='utf-8')
    print(json.dumps(info,indent=2))

if __name__=='__main__':main()
