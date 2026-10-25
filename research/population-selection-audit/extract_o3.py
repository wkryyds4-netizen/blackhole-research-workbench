"""Exactly join the combined O3 subset to the fixed full O3 BBH release."""
from pathlib import Path
import hashlib
import json
import h5py
import numpy as np

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]


def main():
    manifest=json.loads((BASE/'sources/manifest.json').read_text())
    source=next(x for x in manifest if x['file']=='endo3_bbhpop-LIGO-T2100113-v12.hdf5')
    path=BASE/'sources'/source['file']
    assert hashlib.sha256(path.read_bytes()).hexdigest()==source['sha256']
    old=ROOT/'research/population-stage2/source-data/o1+o2+o3_bbhpop_real+semianalytic-LIGO-T2100377-v2.hdf5'
    sha=hashlib.sha256(old.read_bytes()).hexdigest()
    assert sha=='81883c571a9fe0b30747c7bc9daca6d358fcd68ef375c7ace0b93aa9ea2e74fe'
    with h5py.File(path,'r') as h:
        full={k:h['injections'][k][:] for k in h['injections']}
        attrs={k:v.item() for k,v in h.attrs.items()}
        group_attrs={k:v.item() for k,v in h['injections'].attrs.items()}
    assert attrs==group_attrs
    with h5py.File(old,'r') as h:
        combined={k:h['injections'][k][:] for k in h['injections']}
        combined_attrs={k:v.item() for k,v in h.attrs.items()}
    mask=combined['name']==b'o3'
    join=np.rec.fromarrays([full['mass1_source'],full['redshift']])
    target=np.rec.fromarrays([combined['mass1_source'][mask],combined['redshift'][mask]])
    assert len(np.unique(join))==len(join)
    order=np.argsort(join);positions=np.searchsorted(join[order],target)
    assert np.all(positions<len(join))
    index=order[positions]
    assert np.array_equal(join[index],target) and len(np.unique(index))==len(index)
    fields={k:bool(np.array_equal(full[k][index],combined[k][mask],equal_nan=True))
            for k in combined if k in full and k!='mixture_weight'}
    assert all(fields.values())
    assert attrs['total_generated']==attrs['n_accepted']+attrs['n_rejected']
    assert attrs['n_accepted']==len(join)
    ranges={}
    for run in np.unique(combined['name']):
        take=combined['name']==run
        found=np.any(np.stack([combined[k]>1 for k in combined if k.startswith('ifar_')]),axis=0)
        if run in [b'o1',b'o2']:found|=(combined['optimal_snr_net']>10)&take
        ranges[run.decode()]=dict(rows=int(take.sum()),detected_rows=int((take&found).sum()),
            redshift_range=[float(combined['redshift'][take].min()),float(combined['redshift'][take].max())],
            detected_redshift_range=[float(combined['redshift'][take&found].min()),float(combined['redshift'][take&found].max())],
            accepted_above_1p9=int((take&(combined['redshift']>1.9)).sum()),
            detected_above_1p9=int((take&found&(combined['redshift']>1.9)).sum()))
    full_found=np.any(np.stack([full[k]>1 for k in full if k.startswith('ifar_')]),axis=0)
    cache=BASE/'o3-bank.npz'
    np.savez_compressed(cache,**{'full_'+k:v for k,v in full.items()},
                        **{'combined_'+k:v for k,v in combined.items()},combined_o3_to_full=index)
    report=dict(schema='blackhole-o3-injection-source/1',source=source,combined_sha256=sha,
        attributes=attrs,combined_attributes=combined_attrs,root_group_attributes_equal=True,
        full_accepted_rows=len(join),combined_o3_unique_matches=len(index),exact_shared_fields=fields,
        mixture_weight_full_range=[float(full['mixture_weight'].min()),float(full['mixture_weight'].max())],
        mixture_weight_combined_o3_range=[float(combined['mixture_weight'][mask].min()),float(combined['mixture_weight'][mask].max())],
        full_detected_rows=int(full_found.sum()),run_ranges=ranges,
        retained_fraction=float(len(index)/len(join)),
        nominal_time_trial_fraction=float(combined_attrs['total_generated']*attrs['analysis_time_s']/combined_attrs['analysis_time_s']/attrs['total_generated']),
        candidate_reference_fields=[k for k in list(full)+list(attrs) if any(t in k.lower() for t in ['freq','reference','approximant','waveform'])],
        cache_sha256=hashlib.sha256(cache.read_bytes()).hexdigest(),
        interpretation='Exact subset identity, declared proposal bounds and trial/time metadata; not an IID or common-spin-reference certification.')
    (BASE/'o3-provenance.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ['exact_shared_fields','source']},indent=2))


if __name__=='__main__':main()
