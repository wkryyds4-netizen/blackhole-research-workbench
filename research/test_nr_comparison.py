"""Independent checks against the stored raw NR, surrogate and study artifacts."""
from pathlib import Path
import hashlib,json,unittest
import numpy as np
import h5py
from scipy.interpolate import CubicSpline
from scipy.spatial.transform import Rotation
from waveform_reference import rotate_complete_modes

ROOT=Path(__file__).resolve().parent


class RealNRComparisonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.study=json.loads((ROOT/'stage3-results/waveform-comparison.json').read_text(encoding='utf-8'))
        cls.source=json.loads((ROOT/'wave-source-results/sxs-0389-lev3-N2.json').read_text(encoding='utf-8'))
        cls.meta=json.loads((ROOT/'sxs-data/0389/Lev3/metadata.json').read_text())
        cls.nr=np.load(ROOT/'wave-source-results/sxs-0389-lev3-N2.npz')
        cls.sur=np.load(ROOT/'stage3-results/matched-surrogate.npz')
        cls.full=json.loads((ROOT/'stage3-results/metrics-full.json').read_text())
        cls.base=next(x for x in cls.full if x['id']==cls.study['baseline_id'])
        cls.row=next(x for x in cls.study['rows'] if x['id']==cls.study['baseline_id'])

    @classmethod
    def tearDownClass(cls):
        cls.nr.close();cls.sur.close()

    def test_source_identity_and_normalization(self):
        path=ROOT/self.source['source_file']
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),self.study['provenance']['nr_raw_sha256'])
        mass=self.meta['reference_mass1']+self.meta['reference_mass2']
        self.assertAlmostEqual(mass,self.source['normalization_ChMass_code_units'],places=14)
        with h5py.File(path) as f:
            raw=f['Extrapolated_N2.dir/Y_l2_m2.dat'][()]
        np.testing.assert_array_equal(raw[:,0],self.nr['t'])
        np.testing.assert_array_equal(raw[:,1]+1j*raw[:,2],self.nr['h_2_2'])

    def test_reference_clock_and_labels(self):
        p=self.study['protocol'];pars=self.study['provenance']['parameters']
        mass=self.meta['reference_mass1']+self.meta['reference_mass2']
        self.assertAlmostEqual(p['reference_horizon_time_source'],(p['nr_peak_time_M']-4300)*mass,places=11)
        self.assertGreater(p['reference_horizon_time_source'],self.meta['relaxation_time'])
        self.assertEqual(pars['heavy_source_body'],2)
        self.assertAlmostEqual(pars['q'],self.meta['reference_mass2']/self.meta['reference_mass1'],places=14)

    def test_peak_is_interior_stationary_maximum(self):
        t=self.nr['t'];tp=self.study['protocol']['nr_peak_time_M']
        power=sum(abs(self.nr[f'h_{l}_{m}'])**2 for l in range(2,5) for m in range(-l,l+1))
        i=int(np.argmax(power))
        self.assertTrue(t[i-1]<=tp<=t[i+1])
        sp=CubicSpline(t,power)
        self.assertLess(sp(tp,2),0)
        # Newton correction estimates location error in M rather than an
        # amplitude-dependent raw derivative threshold. 1e-3 M is far below
        # the 0.5 M finest waveform grid used for the actual comparison.
        self.assertLess(abs(sp(tp,1)/sp(tp,2)),1e-3)
        self.assertGreaterEqual(sp(tp),power[i]-1e-11)

    def test_complete_rotation_power_invariant(self):
        modes={(l,m):self.nr[f'h_{l}_{m}'][::101] for l in range(2,5) for m in range(-l,l+1)}
        basis=np.asarray(self.study['provenance']['parameters']['coorbital_to_source_matrix'])
        q=Rotation.from_matrix(basis.T).as_quat()[[3,0,1,2]]
        result=rotate_complete_modes(modes,q)
        for l in range(2,5):
            np.testing.assert_allclose(sum(abs(modes[(l,m)])**2 for m in range(-l,l+1)),sum(abs(result[(l,m)])**2 for m in range(-l,l+1)),rtol=1e-13,atol=1e-20)

    def test_all_interpolation_supports(self):
        p=self.study['protocol'];tn=self.nr['t']-p['nr_peak_time_M'];ts=self.sur['t']
        lo,hi=p['interval_M']
        self.assertGreaterEqual(lo,max(tn[0],ts[0]))
        self.assertLessEqual(hi,min(tn[-1],ts[-1]))
        display=np.asarray(self.study['waveform']['time_M'])
        shift=self.row['time_shift_M']
        self.assertGreaterEqual(min(display-shift),ts[0]);self.assertLessEqual(max(display-shift),ts[-1])
        self.assertGreaterEqual(min(display),tn[0]);self.assertLessEqual(max(display),tn[-1])

    def test_independent_frequency_sum_at_reported_optimum(self):
        # Rebuild from raw arrays, without compare_waveforms/common_grid or its PSD/window functions.
        p=self.study['protocol'];r=self.base;cfg=r['config']
        basis=np.asarray(self.study['provenance']['parameters']['coorbital_to_source_matrix'])
        quat=Rotation.from_matrix(basis.T).as_quat()[[3,0,1,2]]
        nr=rotate_complete_modes({(2,m):self.nr[f'h_2_{m}'] for m in range(-2,3)},quat)[(2,2)]
        grid=np.linspace(*p['interval_M'],cfg['samples'])
        a=np.interp(grid,self.nr['t']-p['nr_peak_time_M'],nr)
        b=np.interp(grid,self.sur['t'],self.sur['h_2_2'])
        x=np.linspace(0,1,len(a));alpha=cfg['taper_fraction'];w=np.ones(len(a))
        left=x<alpha/2;right=x>1-alpha/2
        w[left]=.5*(1-np.cos(2*np.pi*x[left]/alpha))
        w[right]=.5*(1-np.cos(2*np.pi*(1-x[right])/alpha))
        self.assertAlmostEqual(float(np.interp(0,grid,w)),1,places=12)
        A=np.fft.fft(a*w,cfg['nfft']);B=np.fft.fft(b*w,cfg['nfft'])
        f=np.fft.fftfreq(cfg['nfft'],cfg['dt_s'])
        keep=(abs(f)>=cfg['fmin_hz'])&(abs(f)<=cfg['fmax_hz'])
        f=f[keep];A=A[keep];B=B[keep]
        y=abs(f)/215.;S=y**(-4.14)-5/y**2+111*(1-y*y+y**4/2)/(1+y*y/2)
        cross=A*np.conj(B)/S
        norm=np.sqrt(sum(abs(A)**2/S)*sum(abs(B)**2/S))
        def match(tau):return sum(cross*np.exp(2j*np.pi*f*tau))/norm
        z=match(r['time_shift_s'])
        self.assertAlmostEqual(abs(z),r['match'],places=13)
        self.assertAlmostEqual(np.angle(z),r['phase_shift_rad'],places=12)
        for delta in [-.1,.1]:
            self.assertGreaterEqual(abs(z)+1e-12,abs(match(r['time_shift_s']+delta*cfg['dt_s'])))

    def test_rows_have_unique_identities_and_no_boundary_optimum(self):
        rows=self.study['rows']
        ids=[r['id'] for r in rows]
        self.assertEqual(len(ids),len(set(ids)))
        self.assertEqual(set(ids),{r['id'] for r in self.full})
        self.assertTrue(all(0<=r['mismatch']<=1 for r in rows))
        self.assertFalse(any(r['at_time_boundary'] for r in rows))

    def test_interpretation_does_not_claim_holdout(self):
        self.assertFalse(self.study['provenance']['independent_validation'])
        self.assertEqual(self.study['provenance']['training_membership'],'unverified')
        self.assertFalse(next(r for r in self.study['rows'] if r['mass_msun']==40 and r['fmin_hz']==20)['start_band_supported'])


if __name__=='__main__':
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(RealNRComparisonTests)
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    (ROOT/'stage3-results/independent-nr-checks.json').write_text(json.dumps({'tests_run':result.testsRun,'errors':len(result.errors),'failures':len(result.failures),'passed':result.wasSuccessful()},indent=2),encoding='utf-8')
    raise SystemExit(0 if result.wasSuccessful() else 1)
