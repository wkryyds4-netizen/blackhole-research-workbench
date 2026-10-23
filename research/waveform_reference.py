"""Explicit parameter and rotation adapters for NR/NRSur waveform comparisons.

Coordinate orbital frames are a labelled approximation to waveform-defined
coorbital frames. This module never equates horizon and radiation clocks by
default and never rotates h22 alone as a spatial scalar.
"""
from pathlib import Path
import sys
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.spatial.transform import Rotation


def reference_masses(metadata):
    prefix = 'reference' if 'reference_mass1' in metadata else 'relaxed'
    m = np.array([metadata[prefix+'_mass1'], metadata[prefix+'_mass2']], float)
    if not np.all(np.isfinite(m)) or np.any(m <= 0):
        raise ValueError('Positive finite reference masses required')
    return m


def orbital_triad(separation, normal):
    """Columns are coorbital (n, L cross n, L) expressed in source axes.

    The normal is projected orthogonal to n before normalization. An actual
    radiation-frame normal can be supplied, but coordinate n is not thereby
    promoted to a radiation-defined orbital phase.
    """
    n = np.asarray(separation, float)
    z = np.asarray(normal, float)
    if n.shape != (3,) or z.shape != (3,) or not np.isfinite([n,z]).all():
        raise ValueError('Finite three-vectors required')
    if np.linalg.norm(n) < 1e-14:
        raise ValueError('Zero separation')
    n = n / np.linalg.norm(n)
    z = z - n*np.dot(n,z)
    if np.linalg.norm(z) < 1e-14:
        raise ValueError('Normal parallel to separation')
    z /= np.linalg.norm(z)
    return np.column_stack([n,np.cross(z,n),z])


def parameters_from_vectors(metadata, spin1, spin2, position1, position2, normal):
    """Sort body labels and project spins without dropping small components."""
    masses = reference_masses(metadata)
    heavy = int(np.argmax(masses))
    spins = np.asarray([spin1,spin2],float)
    positions = np.asarray([position1,position2],float)
    if spins.shape != (2,3) or not np.isfinite(spins).all() or np.any(np.linalg.norm(spins,axis=1)>1):
        raise ValueError('Invalid dimensionless spin vectors')
    basis = orbital_triad(positions[heavy]-positions[1-heavy],normal)
    chi = spins @ basis
    quat = Rotation.from_matrix(basis).as_quat()[[3,0,1,2]]
    return {
        'q':float(masses[heavy]/masses[1-heavy]),
        'chiA0':chi[heavy].tolist(),'chiB0':chi[1-heavy].tolist(),
        'heavy_source_body':heavy+1,'light_source_body':2-heavy,
        'reference_total_mass_source_units':float(masses.sum()),
        'coorbital_to_source_matrix':basis.tolist(),
        'coorbital_to_source_quaternion_wxyz':quat.tolist(),
        'reference_frame_kind':'coordinate orbital triad; not radiation-defined coorbital frame',
    }


def interpolate_horizon_reference(metadata, horizons, source_coordinate_time):
    """horizons: h5py.File. No extrapolation; native time derivatives via spline."""
    t=float(source_coordinate_time)
    positions=[]; velocities=[]; spins=[]
    for name in ['A','B']:
        values={}
        for field in ['CoordCenterInertial.dat','chiInertial.dat']:
            a=np.asarray(horizons[f'Ah{name}.dir/{field}'],float)
            if a.ndim!=2 or a.shape[1]!=4 or not np.isfinite(a).all() or np.any(np.diff(a[:,0])<=0):
                raise ValueError('Invalid horizon samples')
            if not a[0,0] <= t <= a[-1,0]:
                raise ValueError('Horizon extrapolation forbidden')
            values[field]=CubicSpline(a[:,0],a[:,1:],axis=0,extrapolate=False)
        positions.append(values['CoordCenterInertial.dat'](t))
        velocities.append(values['CoordCenterInertial.dat'](t,1))
        spins.append(values['chiInertial.dat'](t))
    r=positions[0]-positions[1];v=velocities[0]-velocities[1]
    omega=np.cross(r,v)/np.dot(r,r)
    result=parameters_from_vectors(metadata,*spins,*positions,omega)
    result.update(source_coordinate_time=t,
        coordinate_orbital_angular_frequency_source_units=omega.tolist(),
        coordinate_proxy_f22_cycles_per_M=float(np.linalg.norm(omega)*reference_masses(metadata).sum()/np.pi),
        interpolation='cubic spline on each native horizon time grid; no extrapolation')
    return result


def source_time_from_radiation_reference(peak_time_source_units, reference_time_M,
                                        total_mass_source_units, *, convention):
    """Reproduce the explicitly gauge-dependent identification in Varma2019.

    Input waveform times must be expressed in original simulation units, not
    already normalized to relaxed M. Output is a horizon coordinate time.
    """
    if convention != 'identify-retarded-and-coordinate-time':
        raise ValueError('Explicit clock-identification convention required')
    vals=np.asarray([peak_time_source_units,reference_time_M,total_mass_source_units])
    if not np.isfinite(vals).all() or total_mass_source_units<=0:
        raise ValueError('Invalid clock conversion')
    return float(peak_time_source_units+reference_time_M*total_mass_source_units)


def rotate_complete_modes(modes, quaternion_wxyz):
    """Rotate source mode components with an active source-to-target rotation.

    All m for every ell=2..ellMax are mandatory. q is [w,x,y,z], length 4
    or shape (4,N). For passive frame change from source axes to basis B,
    pass Rotation.from_matrix(B.T), not B. Delegates the Wigner-D convention
    to the same pinned gwsurrogate code as the evaluated model.
    """
    if not modes:
        raise ValueError('Empty mode set')
    ellmax=max(ell for ell,m in modes)
    keys=[(ell,m) for ell in range(2,ellmax+1) for m in range(-ell,ell+1)]
    if set(modes)!=set(keys):
        raise ValueError('Complete multiplets required; rotating h22 alone is invalid')
    h=np.asarray([modes[k] for k in keys],complex)
    if h.ndim!=2 or not np.isfinite(h).all():
        raise ValueError('Modes must share a finite one-dimensional time grid')
    q=np.asarray(quaternion_wxyz,float)
    if q.shape==(4,):q=np.repeat(q[:,None],h.shape[1],axis=1)
    if q.shape!=(4,h.shape[1]) or not np.isfinite(q).all() or not np.allclose(np.linalg.norm(q,axis=0),1,rtol=0,atol=1e-10):
        raise ValueError('Unit quaternions with shape (4,) or (4,N) required')
    runtime=Path(__file__).resolve().parent/'adapter-build/zig-runtime/gwsurrogate-1.2.0'
    if str(runtime) not in sys.path:sys.path.insert(0,str(runtime))
    from gwsurrogate.new.precessing_surrogate import rotateWaveform
    rotated=rotateWaveform(q,h)
    return dict(zip(keys,rotated))
