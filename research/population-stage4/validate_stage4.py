"""Independent full likelihood product, Jacobian and frozen-input checks."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.special import erf
from scipy.integrate import cumulative_trapezoid

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]

def main():
    expected=json.loads((ROOT/'research/population-stage3/validation.json').read_text())['frozen_hashes']
    for file,sha in expected.items():assert hashlib.sha256((ROOT/file).read_bytes()).hexdigest()==sha
    for rec in json.loads((BASE/'sources/manifest.json').read_text()):
        if 'sha256' in rec:assert hashlib.sha256((BASE/'sources'/rec['file']).read_bytes()).hexdigest()==rec['sha256']
    result=json.loads((BASE/'sensitivity-results.json').read_text())
    frozen=json.loads((ROOT/'results/model.json').read_text())
    with np.load(BASE/'sensitivity-checkpoints.npz',allow_pickle=False) as check:
        new=check['new'];names=check['names'].tolist();target=check['quantiles']
        indices=[next(i for i,e in enumerate(frozen['events']) if e['index']==ix) for ix in result['indices']]
        points=1001
        cdf=lambda x: (1+erf(x/np.sqrt(2)))/2
        negative=lambda s:(cdf(-1/s)-cdf(-2/s))/(.5-cdf(-2/s))
        r=np.linspace(negative(1.2),.5,points)
        errors=[];scaling_errors=[]
        for ni,name in enumerate(names):
            for ki,key in enumerate(result['keys']):
                s=[.3,.6,1.2][int(key[0])]
                f=(r-.5)/(negative(s)-.5)
                rec=frozen['scenarios'][key]
                means=np.array(rec['event_moments'])[:,:2].copy()
                for ei,pos in enumerate(indices):means[pos]=new[ni,ei,ki]
                xi=f*rec['selection'][0]+(1-f)*rec['selection'][1]
                def calc(em):
                    factors=(f[None,:]*em[:,0,None]+(1-f[None,:])*em[:,1,None])/em[:,1,None]
                    density=np.prod(factors.astype(np.longdouble),axis=0)/(xi/rec['selection'][1]).astype(np.longdouble)**len(means)
                    density=np.asarray(density/density.max(),dtype=float)
                    c=cumulative_trapezoid(density,r,initial=0);c/=c[-1]
                    return np.interp([.05,.5,.95],c,r)
                q=calc(means);errors.append(float(np.max(abs(q-target[ni,ki]))))
                scaled=means*np.exp(np.linspace(-4,4,len(means)))[:,None]
                scaling_errors.append(float(np.max(abs(calc(scaled)-q))))
        assert max(errors)<1e-11
        assert max(scaling_errors)<1e-12
    # Direct determinant of detector component masses vs source mass and q.
    m=np.array([3.,20.,80.]);q=np.array([.2,.5,.8]);z=np.array([.1,.8,2.])
    matrices=np.zeros((3,2,2));matrices[:,0,0]=1+z
    matrices[:,1,0]=q*(1+z);matrices[:,1,1]=m*(1+z)
    jac_error=float(np.max(abs(np.linalg.det(matrices)/(m*(1+z)**2)-1)))
    assert jac_error<1e-14
    recovered=json.loads((BASE/'prior-recovery.json').read_text())
    assert recovered['total_samples']==504649 and len(recovered['events'])==73
    assert recovered['maximum_relative_error']<1e-12
    assert all(v['maximum_relative_shape_residual']>1e-3 for e in recovered['events'] for v in e['negative_controls'].values())
    record=dict(schema='blackhole-population-stage4-validation/1',status='pass',
        frozen_input_hashes_verified=expected,source_hashes_verified=True,
        full_product_vs_replacement_log_max_quantile_error=max(errors),
        arbitrary_per_event_density_constant_max_quantile_error=max(scaling_errors),
        mass_coordinate_jacobian_determinant_max_relative_error=jac_error,
        all_73_event_negative_controls_fail=True,
        limit='Algorithmic/numerical validation, not validation of the assumed common Mixed prior or astrophysical calibration.')
    (BASE/'validation.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    print(json.dumps(record,indent=2))

if __name__=='__main__':main()
