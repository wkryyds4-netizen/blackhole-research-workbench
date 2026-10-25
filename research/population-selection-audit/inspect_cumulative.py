"""Extract the fixed O1/O2/O3 cumulative release and audit its O3 row identity."""
from pathlib import Path
import hashlib
import json
import h5py
import numpy as np
from fetch_cumulative import NAME

BASE=Path(__file__).resolve().parent
SHA='89f912ba2e207f56455df9628d4b3fab76cd8d7924ea7f1263a16119bce1ab31'


def main():
    path=BASE/'sources'/NAME
    assert hashlib.sha256(path.read_bytes()).hexdigest()==SHA
    with h5py.File(path,'r') as f:
        a=f['events'][:]
        attrs={k:(v.tolist() if isinstance(v,np.ndarray) else v.item() if isinstance(v,np.generic) else v) for k,v in f.attrs.items()}
    # The release documents nonzero observed SNR exclusively for semianalytic O1/O2.
    early=a['semianalytic_observed_phase_maximized_snr_net']>0
    o3=~early
    assert a['time_geocenter'][early].max()<a['time_geocenter'][o3].min()
    far_names=[k for k in a.dtype.names if k.endswith('_far')]
    assert all(np.isinf(a[k][early]).all() for k in far_names)
    with np.load(BASE/'o3-bank.npz',allow_pickle=False) as b:
        old={k[5:]:b[k] for k in b.files if k.startswith('full_')}
    # Exact unique two-coordinate keys: no tolerance-based nearest-neighbor match.
    def keys(m,z):
        out=np.empty(len(m),dtype=[('m','f8'),('z','f8')]);out['m']=m;out['z']=z
        return out
    ko=keys(old['mass1_source'],old['redshift']);kn=keys(a['mass1_source'][o3],a['redshift'][o3])
    assert len(np.unique(ko))==len(ko) and len(np.unique(kn))==len(kn)
    common,io,jn=np.intersect1d(ko,kn,assume_unique=True,return_indices=True)
    new_ids=np.flatnonzero(o3)[jn]
    physical={k:dict(exact=bool(np.array_equal(old[k][io].astype(float),a[k][new_ids])),max_abs=float(np.max(abs(old[k][io].astype(float)-a[k][new_ids])))) for k in ['mass1_source','mass2_source','redshift']+[f'spin{i}{ax}' for i in [1,2] for ax in 'xyz']} if len(io) else {}
    far={}
    for old_name in sorted(k for k in old if k.startswith('ifar_')):
        new_name='o3_'+old_name[5:]+'_far'
        if new_name not in a.dtype.names:continue
        x=old[old_name][io].astype(float);y=a[new_name][new_ids]
        finite=x>0
        far[old_name]=dict(new_field=new_name,threshold_disagreements=int(np.count_nonzero((x>1)!=(y<1))),
            max_positive_ifar_times_far_error=float(np.max(abs(x[finite]*y[finite]-1))) if finite.any() else None,
            old_zero_new_infinite_agree=bool(np.isinf(y[~finite]).all()))
    columns={k:a[k] for k in a.dtype.names}
    columns.update(early=early,matched_old_o3_rows=io,matched_new_o3_rows=new_ids)
    np.savez_compressed(BASE/'cumulative-bank.npz',**columns)
    summary=dict(schema='blackhole-cumulative-source-audit/1',source_sha256=SHA,source_file=NAME,attributes=attrs,
        counts=dict(all=len(a),semianalytic_o1_o2=int(early.sum()),real_o3=int(o3.sum()),historical_full_bbh_o3=len(ko),exact_join=len(io)),
        time_ranges={label:[float(a['time_geocenter'][sel].min()),float(a['time_geocenter'][sel].max())] for label,sel in [('early',early),('o3',o3)]},
        physical_join=physical,search_join=far,
        ranges={k:[float(a[k].min()),float(a[k].max())] for k in ['redshift','mass1_source','mass2_source','weights']},
        cache_sha256=hashlib.sha256((BASE/'cumulative-bank.npz').read_bytes()).hexdigest(),
        limits='Exact identity of intersecting rows only. Mixture draw densities, time weights and source population coverage may differ. No assigned spin reference frequency.')
    (BASE/'cumulative-provenance.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
