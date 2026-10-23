"""Numerical identities and failure-path tests, without independent NR claims."""
import json
from pathlib import Path
import unittest
import numpy as np
from waveform_metrics import compare_waveforms, common_grid, aligo_historical_psd


class WaveformMetricsTests(unittest.TestCase):
    def setUp(self):
        self.dt = 1/2048
        self.t = np.arange(4096)*self.dt
        self.base = self.packet(self.t)
        self.kw = dict(fmin_hz=20., fmax_hz=800., max_time_shift_s=.05)

    def packet(self, t, width=.15):
        return np.exp(-((t-1)/width)**2)*np.exp(-2j*np.pi*(80*t+15*t*t))

    def metric(self, a, b, **options):
        return compare_waveforms(a, b, self.dt, **{**self.kw, **options})

    def test_identity(self):
        for psd in ['white', 'aligo_historical']:
            r = self.metric(self.base, self.base, psd=psd)
            self.assertLess(r['mismatch'], 2e-15)
            self.assertAlmostEqual(r['time_shift_s'], 0, places=9)

    def test_fractional_delay_phase_amplitude_recovery(self):
        delay, phase, scale = .01337, .73, 2.7
        b = scale*np.exp(1j*phase)*self.packet(self.t-delay)
        for psd in ['white', 'aligo_historical']:
            r = self.metric(self.base, b, psd=psd)
            self.assertLess(r['mismatch'], 1e-11)
            self.assertAlmostEqual(r['time_shift_s'], -delay, places=7)
            self.assertAlmostEqual(r['phase_shift_rad'], -phase, places=5)
            self.assertAlmostEqual(r['amplitude_scale_b_to_a'], 1/scale, places=8)

    def test_real_strain_quadrature_phase(self):
        phase = .6
        r = self.metric(self.base.real, (self.base*np.exp(1j*phase)).real, signal_kind='real_strain')
        self.assertLess(r['mismatch'], 1e-12)
        self.assertIn('not general orbital', r['config']['phase_semantics'])

    def test_shape_difference(self):
        r = self.metric(self.base, self.packet(self.t, width=.04))
        self.assertGreater(r['mismatch'], .2)

    def test_negative_frequency_convention(self):
        # Complex h22 and its conjugate cannot be matched by a constant phase.
        r = self.metric(self.base, np.conj(self.base))
        self.assertGreater(r['mismatch'], .999)

    def test_lag_oversampling_convergence(self):
        b = self.packet(self.t-.01337)
        rs = [self.metric(self.base, b, oversample=k, refine=False) for k in [1, 4, 16]]
        self.assertLessEqual(rs[2]['mismatch'], rs[1]['mismatch']+1e-14)
        self.assertLessEqual(rs[1]['mismatch'], rs[0]['mismatch']+1e-14)
        self.assertLess(abs(rs[2]['time_shift_s']+.01337), self.dt/32)

    def test_padding_convergence(self):
        b = self.packet(self.t, width=.1)
        rs = [self.metric(self.base, b, padding_factor=k, psd='aligo_historical') for k in [2,4,8]]
        self.assertLess(max(r['mismatch'] for r in rs)-min(r['mismatch'] for r in rs), 1e-10)

    def test_amplitude_and_psd_scale_invariance(self):
        b = self.packet(self.t, width=.1)
        x = self.metric(self.base, b)
        y = self.metric(self.base*1e-21, b*3e-22, psd=lambda f: np.ones_like(f)*1e-47)
        self.assertAlmostEqual(x['mismatch'], y['mismatch'], places=13)

    def test_white_zero_lag_against_independent_parseval(self):
        # Narrowband packets have negligible excluded DC/out-of-band power.
        b = self.packet(self.t, width=.1)*np.exp(.4j)
        r = self.metric(self.base, b, taper_fraction=0, max_time_shift_s=0)
        direct = np.vdot(b,self.base)/np.sqrt(np.vdot(self.base,self.base).real*np.vdot(b,b).real)
        self.assertAlmostEqual(r['match'], abs(direct), places=12)
        self.assertAlmostEqual(r['phase_shift_rad'], np.angle(direct), places=12)

    def test_search_boundary_is_reported(self):
        r = self.metric(self.base, self.packet(self.t-.03), max_time_shift_s=.01)
        self.assertTrue(r['at_time_boundary'])
        self.assertAlmostEqual(r['time_shift_s'], -.01, places=8)
        self.assertGreater(r['mismatch'], 0)

    def test_no_circular_endpoint_wrap(self):
        a = np.exp(-((self.t-.1)/.02)**2)*np.exp(-2j*np.pi*100*self.t)
        b = np.exp(-((self.t-1.9)/.02)**2)*np.exp(-2j*np.pi*100*self.t)
        r = self.metric(a, b, max_time_shift_s=.3, taper_fraction=0)
        self.assertGreater(r['mismatch'], .999)

    def test_frequency_and_psd_failures(self):
        for options in [dict(fmin_hz=0), dict(fmax_hz=1024), dict(fmin_hz=900),
                        dict(psd=lambda f: f*0), dict(psd=lambda f: f*np.nan),
                        dict(psd='missing'), dict(psd='aligo_historical', fmin_hz=5)]:
            with self.subTest(options=options), self.assertRaises(ValueError):
                self.metric(self.base, self.base, **options)
        with self.assertRaises(ValueError):
            aligo_historical_psd([0,100])
        self.assertAlmostEqual(float(aligo_historical_psd(215)/1e-49), 33., places=12)

    def test_input_failures(self):
        for a,b,opts in [(self.base, self.base[:-1], {}), (self.base*0,self.base,{}),
                         (self.base,self.base,{'signal_kind':'real_strain'}),
                         (self.base,self.base,{'oversample':0}),
                         (self.base,self.base,{'max_time_shift_s':2})]:
            with self.subTest(opts=opts), self.assertRaises(ValueError):
                self.metric(a,b,**opts)

    def test_common_grid_no_extrapolation(self):
        t,a,b = common_grid([0,1,2,3,4],np.arange(5)+1j*np.arange(5),[1,2,3,4,5],np.arange(1,6),.5)
        np.testing.assert_allclose(t,[1,1.5,2,2.5,3,3.5,4])
        np.testing.assert_allclose(a,t+1j*t)
        np.testing.assert_allclose(b,t)
        with self.assertRaises(ValueError):
            common_grid([0,1,2,3],np.arange(4),[1,2,3,4],np.arange(4),.1,interval=(0,3))

    def test_result_json_serializable(self):
        result = self.metric(self.base, self.base)
        json.dumps(result,allow_nan=False)


if __name__ == '__main__':
    unittest.main(verbosity=2)
