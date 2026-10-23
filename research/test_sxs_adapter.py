"""Independent checks against original data; does not regenerate expectations."""
import unittest, json
import numpy as np
import h5py
from pathlib import Path
ROOT=Path(__file__).resolve().parent
PUBLIC=ROOT.parent/'simulator/public/datasets'

class SourceRoundtrip(unittest.TestCase):
    def test_all_native_knots(self):
        for item in json.loads((PUBLIC/'catalog.json').read_text(encoding='utf-8')):
            with self.subTest(dataset=item['id']):
                d=json.loads((PUBLIC/item['url'].split('/')[-1]).read_text(encoding='utf-8'))
                p=d['provenance']; mass=p['reference_total_mass_in_source_units']; zero=p['source_common_horizon_time']
                _,sim,level=item['id'].split('-'); level='Lev'+level[3:]
                rotation=np.asarray(p['rotation_matrix'])
                np.testing.assert_allclose(rotation@rotation.T,np.eye(3))
                self.assertAlmostEqual(np.linalg.det(rotation),1)
                self.assertGreater(len(d['frames']),100)
                self.assertEqual(next(f for f in d['frames'] if f['t']>=0)['bodies'][0]['id'],'C')
                with h5py.File(ROOT/'sxs-data'/sim/level/'Horizons.h5') as source:
                    for name in 'ABC':
                        frames=[f for f in d['frames'] if any(b['id']==name for b in f['bodies'])]
                        times=np.array([f['t'] for f in frames])*mass+zero
                        bodies=[next(b for b in f['bodies'] if b['id']==name) for f in frames]
                        for field,file,tol in [('position','CoordCenterInertial.dat',.003),('mass','ChristodoulouMass.dat',1e-5),('spin','chiInertial.dat',1e-4)]:
                            a=source[f'Ah{name}.dir/{file}'][:]
                            display=np.array([b[field] for b in bodies])
                            if field in ['position','spin']: display=display@rotation
                            else: display=display[:,None]
                            if field!='spin': display*=mass
                            # Source->export roundtrip at every retained frame.
                            expected=np.column_stack([np.interp(times,a[:,0],a[:,j]) for j in range(1,a.shape[1])])
                            np.testing.assert_allclose(display,expected,atol=2e-10,rtol=1e-10)
                            # Every source knot in retained interval, with the
                            # advertised bound evaluated in display coordinates.
                            mask=(a[:,0]>=times[0])&(a[:,0]<=times[-1]); native=a[mask]
                            reconstructed=np.column_stack([np.interp(native[:,0],times,display[:,j]) for j in range(display.shape[1])])
                            delta=reconstructed-native[:,1:]
                            if field!='spin':delta/=mass
                            self.assertLessEqual(float(np.abs(delta).max()),tol+2e-10)

if __name__=='__main__': unittest.main()
