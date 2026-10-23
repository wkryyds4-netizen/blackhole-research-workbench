"""Independent invariants and convention checks for waveform_reference."""
import json
from pathlib import Path
import unittest
import numpy as np
import h5py
from scipy.spatial.transform import Rotation
from waveform_reference import (orbital_triad, parameters_from_vectors,
    interpolate_horizon_reference, source_time_from_radiation_reference,
    rotate_complete_modes)


class ReferenceTests(unittest.TestCase):
    def test_body_exchange_preserves_spins(self):
        meta={'reference_mass1':.4,'reference_mass2':.6}
        r=parameters_from_vectors(meta,[.1,.2,.3],[.2,.1,-.3],[1,0,0],[-1,0,0],[0,0,1])
        self.assertEqual(r['heavy_source_body'],2)
        self.assertAlmostEqual(r['q'],1.5)
        np.testing.assert_allclose(r['chiA0'],[-.2,-.1,-.3])
        np.testing.assert_allclose(np.asarray(r['coorbital_to_source_matrix'])@r['chiB0'],[.1,.2,.3])

    def test_triad_is_proper(self):
        b=orbital_triad([1,2,3],[-.2,.3,1])
        np.testing.assert_allclose(b.T@b,np.eye(3),atol=1e-14)
        self.assertAlmostEqual(np.linalg.det(b),1)

    def test_clock_requires_convention(self):
        with self.assertRaises(ValueError):source_time_from_radiation_reference(5000,-4300,1,convention='physical equality')
        self.assertEqual(source_time_from_radiation_reference(5000,-4300,.99,convention='identify-retarded-and-coordinate-time'),743.)

    def test_real_0389_metadata(self):
        folder=Path(__file__).resolve().parent/'sxs-data/0389/Lev3'
        meta=json.loads((folder/'metadata.json').read_text())
        with h5py.File(folder/'Horizons.h5') as f:
            r=interpolate_horizon_reference(meta,f,meta['reference_time'])
        self.assertEqual(r['heavy_source_body'],2)
        self.assertGreater(np.linalg.norm(r['chiA0']),5e-5)
        self.assertLess(np.linalg.norm(r['chiA0']),6e-5)

    def test_no_extrapolation(self):
        folder=Path(__file__).resolve().parent/'sxs-data/0389/Lev3'
        meta=json.loads((folder/'metadata.json').read_text())
        with h5py.File(folder/'Horizons.h5') as f:
            with self.assertRaises(ValueError):interpolate_horizon_reference(meta,f,-1e20)

    def test_rotation_rejects_incomplete_multiplet(self):
        with self.assertRaises(ValueError):rotate_complete_modes({(2,2):np.ones(3)},[1,0,0,0])

    def test_wigner_z_phase_sign(self):
        h={(2,m):np.ones(3,dtype=complex) for m in range(-2,3)}
        theta=.42;q=[np.cos(theta/2),0,0,np.sin(theta/2)]
        rotated=rotate_complete_modes(h,q)
        for (ell,m),v in rotated.items():
            np.testing.assert_allclose(v,np.exp(-1j*m*theta),atol=1e-14)

    def test_wigner_rotation_preserves_power_and_roundtrip(self):
        rng=np.random.default_rng(7)
        h={(l,m):rng.normal(size=11)+1j*rng.normal(size=11) for l in range(2,5) for m in range(-l,l+1)}
        q=Rotation.from_rotvec([.4,.2,-.3]).as_quat()[[3,0,1,2]]
        rotated=rotate_complete_modes(h,q)
        np.testing.assert_allclose(sum(abs(v)**2 for v in h.values()),sum(abs(v)**2 for v in rotated.values()),rtol=1e-14)
        q[1:]*=-1
        recovered=rotate_complete_modes(rotated,q)
        for k in h:np.testing.assert_allclose(recovered[k],h[k],atol=1e-13)


if __name__=='__main__':unittest.main()
