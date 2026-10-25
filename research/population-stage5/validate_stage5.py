"""Independent full-likelihood product and distance-Jacobian checks."""
from pathlib import Path
import sys
import json
import hashlib
import numpy as np
from scipy.integrate import cumulative_trapezoid

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
sys.path.insert(0,str(ROOT/'research/population-stage4'))
from recover_prior import LAL

def main():
    frozen=json.loads((ROOT/'results/model.json').read_text())
    result=json.loads((BASE/'tail-results.json').read_text())
    expected=json.loads((ROOT/'research/population-stage3/validation.json').read_text())['frozen_hashes']
    for file,sha in expected.items():assert hashlib.sha256((ROOT/file).read_bytes()).hexdigest()==sha
    for x in json.loads((BASE/'sources/manifest.json').read_text()):
        assert hashlib.sha256((BASE/'sources'/x['file']).read_bytes()).hexdigest()==x['sha256']
    dl=json.loads((BASE/'download.json').read_text());path=BASE/'source-data'/dl['file']
    assert hashlib.sha256(path.read_bytes()).hexdigest()==dl['sha256']
    pos=next(i for i,e in enumerate(frozen['events']) if e['index']==38)
    errors=[]
    with np.load(BASE/'tail-checkpoints.npz',allow_pickle=False) as cp:
        r=cp['r'];fs=cp['fs'];means=cp['means']
        for ni,name in enumerate(cp['names'].tolist()):
            for ki,key in enumerate(result['keys']):
                rec=frozen['scenarios'][key]
                em=np.array(rec['event_moments'])[:,:2].copy();em[pos]=means[ni,ki]
                f=fs[ki]
                factors=(f[None,:]*em[:,0,None]+(1-f[None,:])*em[:,1,None])/em[:,1,None]
                xi=(f*rec['selection'][0]+(1-f)*rec['selection'][1])/rec['selection'][1]
                p=np.prod(factors.astype(np.longdouble),axis=0)/xi.astype(np.longdouble)**len(em)
                p=np.asarray(p/p.max(),dtype=float)
                c=cumulative_trapezoid(p,r,initial=0);c/=c[-1]
                q=np.interp([.05,.5,.95],c,r)
                errors.append(float(np.max(abs(q-np.array(result['scenarios'][name]['quantiles'][ki])))))
        assert max(errors)<1e-11
        for label,out in [('original','original'),('cosmo','cosmo_lal'),('nocosmo','nocosmo_lal')]:
            a=cp[f'sorted_fraction_{label}'];assert abs(a.sum()-1)<1e-12
            assert abs((1/np.sum(a*a))/result['weight_diagnostics'][out][0]['ess']-1)<1e-12
    z=np.linspace(.04,2.5,51);h=1e-5
    derivative=(LAL.luminosity_distance(z+h).value-LAL.luminosity_distance(z-h).value)/(2*h)
    analytic=LAL.luminosity_distance(z).value/(1+z)+(1+z)*LAL.hubble_distance.value/LAL.efunc(z)
    error=float(np.max(abs(derivative/analytic-1)));assert error<1e-8
    match=json.loads((BASE/'representation-match.json').read_text())
    assert all(match['exact_matches'].values()) and match['config_equal'] and match['analytic_prior_strings_equal']
    record=dict(schema='blackhole-population-stage5-validation/1',status='pass',
        frozen_hashes_verified=expected,source_and_download_hashes_verified=True,
        independent_full_product_max_quantile_error=max(errors),
        finite_difference_distance_jacobian_max_relative_error=error,
        sorted_weights_reproduce_ess=True,
        identity_error=result['checks']['change_of_measure_max_relative_error'],
        limit='Numerical/source validation; not independent posterior convergence or a calibrated waveform uncertainty.')
    (BASE/'validation.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    print(json.dumps(record,indent=2))

if __name__=='__main__':main()
