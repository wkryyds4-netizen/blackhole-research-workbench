"""Auditable, no hidden time shifts/rotations, SXS waveform adapter.

Old NRAR products require numpy+h5py. RPDMB additionally requires sxs==2024.0.25
in simulator/.waveform-venv; Julia is not needed or imported by this code path.
Exports original archive times and all available complex modes, never a strain
projection. All subsequent peak shifts/alignment must be logged by the caller.
"""
from pathlib import Path
import argparse, hashlib, json, re
import numpy as np
import h5py

ROOT = Path(__file__).resolve().parent

def load_nr_waveform(path, group='Extrapolated_N2.dir', metadata_path=None):
    path=Path(path)
    with h5py.File(path) as f:
        compressed='sxs_format' in f.attrs
    if compressed:
        import sxs
        md=sxs.Metadata.from_file(metadata_path) if metadata_path else None
        w=sxs.rpdmb.load(path, ignore_validation=False, check_md5=True,
                        transform_to_inertial=True, drop_times_before='begin', metadata=md)
        t=np.array(w.t); modes={tuple(map(int,lm)):np.array(w.data[:,i]) for i,lm in enumerate(w.LM)}
        info={'reader':f'sxs {sxs.__version__} official RPDMB', 'format':'rotating_paired_diff_multishuffle_bzip2',
              'frame':'inertial (official corotating-to-inertial transform)',
              'archive_json':json.loads(path.with_suffix('.json').read_text())}
    else:
        with h5py.File(path) as f:
            g=f[group]; modes={}; t=None
            for key in g:
                match=re.fullmatch(r'Y_l(\d+)_m(-?\d+)\.dat',key)
                if not match: continue
                data=g[key][()]
                if t is None: t=data[:,0]
                else: assert np.array_equal(t,data[:,0]),'Mode times differ'
                modes[tuple(map(int,match.groups()))]=data[:,1]+1j*data[:,2]
            history=g['History.txt'][()].decode()
            attrs={k:np.asarray(v).tolist() for k,v in g.attrs.items()}
        assert attrs['FrameType']==1 and attrs['MIsScaledOut']==1 and attrs['RIsScaledOut']==1
        mass_match=re.search(r"D\['ChMass'\] = ([0-9.eE+-]+)",history)
        info={'reader':'direct documented NRAR columns [t,real,imag]', 'format':'NRAR',
              'group':group,'frame':'inertial','hdf5_attributes':attrs,
              'normalization_ChMass_code_units':float(mass_match.group(1)) if mass_match else None,
              'history':history}
    assert len(t)>100 and np.all(np.isfinite(t)) and np.all(np.diff(t)>0)
    assert all(np.all(np.isfinite(v)) and v.shape==t.shape for v in modes.values())
    expected={(ell,m) for ell in range(min(l for l,m in modes),max(l for l,m in modes)+1) for m in range(-ell,ell+1)}
    assert set(modes)==expected,'Incomplete mode multiplets'
    info.update(schema='blackhole-nr-waveform/1',source_file=str(path.relative_to(ROOT)),
        source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),time_unit='GM/c^3',mode_unit='r*h/M',
        spin_weight=-2,complex_convention='h_plus - i h_cross',time_shift_applied=0,
        extra_spatial_rotation_applied=False,samples=len(t),mode_count=len(modes),
        time_range_M=[float(t[0]),float(t[-1])],
        peak_all_modes_time_M=float(t[np.argmax(sum(abs(v)**2 for v in modes.values()))]),
        peak_lte4_time_M=float(t[np.argmax(sum(abs(v)**2 for (l,m),v in modes.items() if l<=4))]),
        peak_h22_time_M=float(t[np.argmax(abs(modes[(2,2)]))]),
        clock_warning='Radiation retarded-time coordinate; no exact mapping to horizon coordinate time is asserted.',
        training_membership='Not established; external archive comparison, not a held-out validation claim.')
    return t,modes,info

def export(path,output,group='Extrapolated_N2.dir',metadata_path=None):
    t,modes,info=load_nr_waveform(path,group,metadata_path)
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output,t=t,**{f'h_{l}_{m}':h for (l,m),h in modes.items()})
    info['npz_sha256']=hashlib.sha256(output.read_bytes()).hexdigest()
    output.with_suffix('.json').write_text(json.dumps(info,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps({k:info[k] for k in ['source_file','samples','mode_count','time_range_M','peak_lte4_time_M']},indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--dataset',choices=['0389','2327'],default='0389')
    args=parser.parse_args()
    if args.dataset=='0389':
        path=ROOT/'waveform-data/0389/Lev3/rhOverM_Asymptotic_GeometricUnits_CoM.h5'
        for n in [2,3,4]: export(path,ROOT/f'wave-source-results/sxs-0389-lev3-N{n}.npz',f'Extrapolated_N{n}.dir')
    else:
        for lev in [2,3]:
            export(ROOT/f'waveform-data/2327/Lev{lev}/Strain_N2.h5',ROOT/f'wave-source-results/sxs-2327-lev{lev}-N2.npz',metadata_path=ROOT/f'sxs-data/2327/Lev{lev}/metadata.json')
