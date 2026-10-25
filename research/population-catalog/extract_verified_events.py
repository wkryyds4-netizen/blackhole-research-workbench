"""Match tutorial rows to verified original-prior catalog files, label by label."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np
import h5py

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
sys.path.insert(0,str(ROOT/'research/population-stage7'))
from extract_events import join,prior_audit,read,PARAMS,SOURCE
from diagnose_alias import NATIVE

EXTRA=['mass_1','mass_2','mass_2_source','luminosity_distance','chirp_mass','phi_12','tilt_1','tilt_2','chi_eff',
       'spin_1x','spin_1y','spin_1z','spin_2x','spin_2y','spin_2z',
       'cos_tilt_1_infinity_only_prec_avg','cos_tilt_2_infinity_only_prec_avg']


def main():
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--indices',nargs='+',type=int);group.add_argument('--available',action='store_true');args=parser.parse_args()
    downloads=json.loads((BASE/'download-results.json').read_text())
    output=BASE/'event-provenance.json';results=json.loads(output.read_text()) if output.exists() else []
    assert hashlib.sha256((ROOT/'data/posteriors.npz').read_bytes()).hexdigest()=='3aec709780ac373e485682323c20cbe6ed0c14028c66b030bb73d848739489c4'
    with np.load(ROOT/'data/posteriors.npz',allow_pickle=False) as archive:
        for index in ([r['index'] for r in downloads] if args.available else args.indices):
            rec=next(r for r in downloads if r['index']==index)
            existing=next((r for r in results if r['index']==index and r['source_sha256']==rec['sha256'] and r.get('native_source_identity_resolved')),None)
            if args.available and existing and (BASE/f'event-{index}.npz').exists():continue
            path=ROOT/rec.get('local_file','research/population-catalog/source-data/'+rec['source']['file'])
            assert hashlib.sha256(path.read_bytes()).hexdigest()==rec['sha256']
            local=[archive[f'{index}_{k}'] for k in PARAMS];labels=[];maps={}
            with h5py.File(path,'r') as h:
                for label,g in h.items():
                    if not isinstance(g,h5py.Group) or 'posterior_samples' not in g:continue
                    ds=g['posterior_samples']
                    if not isinstance(ds,h5py.Dataset) or not set(SOURCE)<=set(ds.dtype.names or []):continue
                    data=ds.fields(SOURCE)[:]
                    mapping,info=join(local,[data[k] for k in SOURCE]);labels.append(dict(label=label,**info));maps[label]=mapping
                full=[r for r in labels if r['matches']==len(local[0]) and r['ambiguous']==0]
                row=dict(index=index,event=rec['event'],source_file=rec['source']['file'],source_sha256=rec['sha256'],label_matches=labels,
                    uniquely_verified_tutorial_source=len(full)==1)
                resolved=len(full)==1;canonical=full[0]['label'] if resolved else None
                if len(full)>1 and index==14:
                    diagnosis=BASE/'alias-source-diagnosis.json';alias=json.loads(diagnosis.read_text())
                    assert alias['source_sha256']==rec['sha256'] and alias['native_posterior_equivalent'] and all(alias['native_columns_exact'].values())
                    assert {r['label'] for r in full}=={alias['canonical'],alias['alias']}
                    canonical=alias['canonical'];resolved=True
                    row['label_resolution']=dict(scope='Native posterior equivalence only; multiple labels remain. Canonical analysis has explicit waveform/frequency/prior metadata.',
                        diagnosis_sha256=hashlib.sha256(diagnosis.read_bytes()).hexdigest(),derived_columns_differ=True)
                elif len(full)>1 and index==48:
                    path=BASE/'cosmo-lineage-resolution.json'
                    if path.exists():
                        evidence=json.loads(path.read_text())
                        assert evidence['nocosmo_source']['sha256']==rec['sha256']
                        assert evidence['three_way_tutorial_rows_exact'] and all(evidence['native_columns_exact'].values())
                        assert evidence['config_equal'] and evidence['original_analytic_prior_declarations_equal']
                        canonical=evidence['canonical'];assert canonical in [r['label'] for r in full]
                        resolved=True
                        row['label_resolution']=dict(scope=evidence['resolution'],evidence_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                            containing_label_remains_ambiguous=True,mixed_full_posterior_equivalent=False)
                elif len(full)>1 and index==72:
                    path=BASE/'mixture-component-resolution.json'
                    if path.exists():
                        evidence=json.loads(path.read_text())
                        assert evidence['source_sha256']==rec['sha256'] and evidence['complete_mixed_partition']
                        assert evidence['tutorial_origin_waveform_unique'] and evidence['disjoint_analysis_input_points']
                        assert all(evidence['canonical_native_columns_exact'].values()) and all(evidence['other_native_columns_exact'].values())
                        canonical=evidence['canonical'];assert canonical in [r['label'] for r in full]
                        resolved=True
                        row['label_resolution']=dict(scope=evidence['resolution'],evidence_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                            containing_label_remains_ambiguous=True,mixed_full_posterior_equivalent=False)
                elif len(full)>1 and any(v['label']=='C01:IMRPhenomXPHM' for v in full):
                    canonical='C01:IMRPhenomXPHM';cds=h[canonical+'/posterior_samples'];equivalences=[]
                    if set(NATIVE)<=set(cds.dtype.names):
                        native=cds.fields(NATIVE)[:]
                        for candidate in full:
                            if candidate['label']==canonical:continue
                            other=h[candidate['label']+'/posterior_samples']
                            if len(other)!=len(cds) or not set(NATIVE)<=set(other.dtype.names):
                                equivalences.append(dict(label=candidate['label'],native_equivalent=False,reason='Shape or input fields differ'));continue
                            alt=other.fields(NATIVE)[:];m,joined=join([native[k] for k in SOURCE],[alt[k] for k in SOURCE])
                            complete=joined['matches']==len(native) and joined['ambiguous']==0 and len(np.unique(m))==len(m)
                            equality={k:bool(np.array_equal(native[k],alt[k][m],equal_nan=True)) for k in NATIVE} if complete else {}
                            equivalences.append(dict(label=candidate['label'],join=joined,native_columns_exact=equality,native_equivalent=complete and all(equality.values())))
                        resolved=len(equivalences)==len(full)-1 and all(v['native_equivalent'] for v in equivalences)
                    row['label_resolution']=dict(scope='Native inputs only; derived-field equivalence is not assumed.',equivalences=equivalences)
                row['native_source_identity_resolved']=resolved
                if resolved:
                    label=canonical;g=h[label];ds=g['posterior_samples'];fields=list(dict.fromkeys(SOURCE+[k for k in EXTRA if k in ds.dtype.names]))
                    samples=ds.fields(fields)[:];mapping=maps[label]
                    assert len(np.unique(mapping))==len(mapping)
                    exact={k:bool(np.array_equal(samples[k][mapping],v)) for k,v in zip(SOURCE,local)};assert all(exact.values())
                    metadata={k:read(g,k) for k in ['meta_data/meta_data/f_ref','config_file/config/reference-frequency','meta_data/meta_data/cosmology']}
                    audit=prior_audit(g,samples)
                    payload={'tutorial_source_rows':mapping}
                    for key in fields:
                        target='mass_1' if key=='mass_1_source' else 'mass_1_detector' if key=='mass_1' else key
                        payload[target]=samples[key]
                    target=BASE/f'event-{index}.npz';np.savez_compressed(target,**payload)
                    row.update(analysis_label=label,rows=len(samples),exact_coordinates=exact,metadata=metadata,prior_audit=audit,
                        phi12_recovered='phi_12' in fields,infinity_tilts_present=all(k in fields for k in ['cos_tilt_1_infinity_only_prec_avg','cos_tilt_2_infinity_only_prec_avg']),
                        cache_sha256=hashlib.sha256(target.read_bytes()).hexdigest())
            results=[r for r in results if r['index']!=index]+[row];output.write_text(json.dumps(results,indent=2),encoding='utf-8')
            print(json.dumps({k:v for k,v in row.items() if k not in ['prior_audit','exact_coordinates']},indent=2),flush=True)


if __name__=='__main__':main()
