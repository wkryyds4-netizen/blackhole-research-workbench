"""Extract original azimuths and released infinity tilts without changing weights."""
from pathlib import Path
import hashlib
import json
import h5py
import numpy as np

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
FIELDS=['mass_1','mass_2','mass_1_source','mass_2_source','mass_ratio','redshift',
        'a_1','a_2','tilt_1','tilt_2','cos_tilt_1','cos_tilt_2','phi_12','chi_eff',
        'spin_1x','spin_1y','spin_1z','spin_2x','spin_2y','spin_2z',
        'tilt_1_infinity_only_prec_avg','tilt_2_infinity_only_prec_avg',
        'cos_tilt_1_infinity_only_prec_avg','cos_tilt_2_infinity_only_prec_avg',
        'chi_eff_infinity_only_prec_avg']


def main():
    records=json.loads((ROOT/'research/population-stage7/event-provenance.json').read_text())
    old=json.loads((ROOT/'research/population-stage5/download.json').read_text())
    oldmatch=json.loads((ROOT/'research/population-stage5/representation-match.json').read_text())
    results=[];payload={};historical=[]
    for index in [19,20,25,38,50]:
        if index==38:
            candidates=list((ROOT/'research/population-stage5/source-data').glob('*GW190805*nocosmo.h5'))
            assert len(candidates)==1
            path=candidates[0];sha=oldmatch['source_hashes']['nocosmo'];fref=oldmatch['metadata']['meta_data/meta_data/f_ref'][0]
            name='GW190805_211137';label=oldmatch['label']
        else:
            rec=next(x for x in records if x['index']==index);source=rec['sources']['nocosmo']
            path=ROOT/'research/population-stage7/source-data'/source['file'];sha=source['sha256']
            fref=rec['metadata']['meta_data/meta_data/f_ref'][0];name=rec['event'];label=rec['analysis_label']
        assert hashlib.sha256(path.read_bytes()).hexdigest()==sha
        with h5py.File(path,'r') as h:
            g=h[label];sample=g['posterior_samples'];missing=[k for k in FIELDS if k not in sample.dtype.names]
            data={k:sample[k].astype(float) for k in FIELDS if k in sample.dtype.names}
            packages={row[0].decode():row[1].decode() for row in h['version/packages'][()]}
            historical.append(dict(file=path.name,sha256=sha,
                backward_spin_evolution=[v.decode() for v in g['meta_data/meta_data/backward_spin_evolution'][()]],
                pesummary=[v.decode() for v in h['version/pesummary'][()]],lalsuite=packages.get('lalsuite')))
        n=len(data['mass_1']);payload.update({f'{index}_{k}':v for k,v in data.items()})
        payload[f'{index}_source_row']=np.arange(n,dtype=np.int64)
        finite={k:int(np.isfinite(v).sum()) for k,v in data.items()}
        current=(data['a_1']*data['cos_tilt_1']+data['mass_ratio']*data['a_2']*data['cos_tilt_2'])/(1+data['mass_ratio'])
        infkeys=['cos_tilt_1_infinity_only_prec_avg','cos_tilt_2_infinity_only_prec_avg']
        checks=dict(native_chieff_max_abs_error=float(np.max(abs(current-data['chi_eff']))))
        assert checks['native_chieff_max_abs_error']<1e-12
        phi=np.arctan2(data['spin_2y'],data['spin_2x'])-np.arctan2(data['spin_1y'],data['spin_1x'])
        checks['azimuth_cosine_max_abs_error']=float(np.max(abs(np.cos(phi)-np.cos(data['phi_12']))))
        assert checks['azimuth_cosine_max_abs_error']<1e-10
        if all(k in data for k in infkeys):
            infinity=(data['a_1']*data[infkeys[0]]+data['mass_ratio']*data['a_2']*data[infkeys[1]])/(1+data['mass_ratio'])
            valid=np.isfinite(infinity)
            checks['infinity_finite_rows']=int(valid.sum())
            checks['infinity_chieff_max_abs_drift']=float(np.max(abs(infinity[valid]-current[valid])))
            checks['stored_infinity_chieff_max_abs_error']=float(np.max(abs(infinity[valid]-data['chi_eff_infinity_only_prec_avg'][valid])))
            for i,k in enumerate(infkeys,1):
                checks[f'cos_tilt_{i}_max_abs_change']=float(np.max(abs(data[k][valid]-data[f'cos_tilt_{i}'][valid])))
        results.append(dict(index=index,event=name,analysis_label=label,source_file=path.name,source_sha256=sha,
            native_spin_reference_hz=fref,rows=n,missing_fields=missing,finite_counts=finite,checks=checks))
    cache=BASE/'spin-bank.npz';np.savez_compressed(cache,**payload)
    with np.load(ROOT/'data/posteriors.npz',allow_pickle=False) as tutorial:
        available=sorted({k.split('_',1)[1] for k in tutorial.files if k.split('_',1)[0].isdigit()})
    report=dict(schema='blackhole-spin-reference-source/1',events=results,tutorial_fields=available,
        tutorial_phi12_present='phi_12' in available,tutorial_reference_frequency_present=any('ref' in k or 'frequency' in k for k in available),
        cache_sha256=hashlib.sha256(cache.read_bytes()).hexdigest(),
        interpretation='Official derived infinity fields recovered for available source events. No new physical-reference population inference; injection reference remains unresolved.')
    (BASE/'spin-provenance.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    (BASE/'historical-metadata.json').write_text(json.dumps(historical,indent=2),encoding='utf-8')
    print(json.dumps(dict(events=[{k:e[k] for k in ['event','rows','missing_fields','checks']} for e in results],tutorial_fields=available),indent=2))


if __name__=='__main__':main()
