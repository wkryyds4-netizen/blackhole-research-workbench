"""Source identity and rotation-invariant checks for actual public NR data."""
from pathlib import Path
import hashlib,json,sys,importlib.metadata
import numpy as np
import h5py
import sxs
ROOT=Path(__file__).resolve().parent
checks=[]
for order in [2,3,4]:
    exported=np.load(ROOT/f'wave-source-results/sxs-0389-lev3-N{order}.npz')
    with h5py.File(ROOT/'waveform-data/0389/Lev3/rhOverM_Asymptotic_GeometricUnits_CoM.h5') as f:
        group=f[f'Extrapolated_N{order}.dir']
        for key in exported.files:
            if key=='t':continue
            _,ell,m=key.split('_'); source=group[f'Y_l{ell}_m{m}.dat'][()]
            assert np.array_equal(source[:,0],exported['t'])
            assert np.array_equal(source[:,1]+1j*source[:,2],exported[key])
    checks.append(dict(dataset=f'0389 Lev3 N{order}',exact_source_modes=77))

for lev in [2,3]:
    path=ROOT/f'waveform-data/2327/Lev{lev}/Strain_N2.h5'
    metadata=sxs.Metadata.from_file(ROOT/f'sxs-data/2327/Lev{lev}/metadata.json')
    rotating=sxs.rpdmb.load(path,ignore_validation=False,transform_to_inertial=False,
                          drop_times_before='begin',metadata=metadata)
    exported=np.load(ROOT/f'wave-source-results/sxs-2327-lev{lev}-N2.npz')
    assert np.array_equal(rotating.t,exported['t'])
    errors=[]
    for ell in range(2,9):
        inertial_power=sum(abs(exported[f'h_{ell}_{m}'])**2 for m in range(-ell,ell+1))
        rotating_power=sum(abs(rotating.data[:,rotating.index(ell,m)])**2 for m in range(-ell,ell+1))
        # Absolute errors normalized to each multiplet's peak, robust to tiny late tails.
        errors.append(float(np.max(abs(inertial_power-rotating_power))/np.max(rotating_power)))
    assert max(errors)<1e-12,errors
    assert np.max(abs(np.linalg.norm(rotating.frame.ndarray,axis=1)-1))<1e-12
    checks.append(dict(dataset=f'2327 Lev{lev} N2',time_exact=True,
                       multiplet_power_rotation_max_relative_error=max(errors)))
assert 'juliacall' not in sys.modules and 'sxs.julia' not in sys.modules
report=dict(status='passed',checks=checks,julia_imported=False,
    environment={p:importlib.metadata.version(p) for p in ['sxs','numpy','scipy','h5py','quaternionic','spherical','numba']},
    sxs_wheel_sha256=hashlib.sha256((ROOT/'adapter-sources/sxs-2024.0.25-py3-none-any.whl').read_bytes()).hexdigest())
(ROOT/'wave-source-results/validation.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
