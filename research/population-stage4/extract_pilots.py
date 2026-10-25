"""Export safe numeric arrays and selected metadata from validated official HDF."""
from pathlib import Path
import hashlib
import json
import h5py
import numpy as np

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
LABELS={'phenom':'C01:IMRPhenomXPHM','mixed':'C01:Mixed','seob':'C01:SEOBNRv4PHM'}
MAPPING={'mass_1':'mass_1_source','mass_ratio':'mass_ratio','redshift':'redshift',
         'a_1':'a_1','a_2':'a_2','cos_tilt_1':'cos_tilt_1','cos_tilt_2':'cos_tilt_2',
         'chirp_mass_detector':'chirp_mass','luminosity_distance':'luminosity_distance'}

def read(group,path):
    if path not in group or not isinstance(group[path],h5py.Dataset):return None
    return [x.decode() if isinstance(x,bytes) else float(x) for x in np.atleast_1d(group[path][()])]

def main():
    match=json.loads((ROOT/'research/population-stage3/pe-pilot-match.json').read_text())
    bank={};records=[]
    for event in match['events']:
        index=event['original_index']
        path=ROOT/'research/population-stage3/source-data'/event['source_file']
        assert hashlib.sha256(path.read_bytes()).hexdigest()==event['source_sha256']
        with h5py.File(path,'r') as h:
            for short,label in LABELS.items():
                g=h[label];d=g['posterior_samples'][:]
                for k,v in MAPPING.items():
                    bank[f'{index}_{short}_{k}']=d[v]
                    assert np.isfinite(d[v]).all()
                records.append(dict(index=index,event=event['event'],label=label,short=short,n=len(d),
                    f_ref=read(g,'meta_data/meta_data/f_ref'),
                    cosmology=read(g,'meta_data/meta_data/cosmology'),
                    engine={k:read(g,'config_file/engine/'+k) for k in ['chirpmass-min','chirpmass-max','q-min','a_spin1-max','a_spin2-max','distance-max','fref']},
                    observed_chirp_mass_detector=[float(d['chirp_mass'].min()),float(d['chirp_mass'].max())],
                    source_sha256=event['source_sha256']))
    np.savez_compressed(BASE/'pilot-bank.npz',**bank)
    (BASE/'pilot-bank-metadata.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
    print([(r['index'],r['short'],r['n']) for r in records])

if __name__=='__main__':main()
