"""Full-product check, original-row checksums and representation regression."""
from pathlib import Path
import argparse
import json
import hashlib
import sys
import numpy as np
from scipy.integrate import cumulative_trapezoid

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
sys.path.insert(0,str(ROOT/'research/population-stage4'))
from recover_prior import LAL


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--prefix',default='');args=parser.parse_args()
    result=json.loads((BASE/(args.prefix+'comparison-results.json')).read_text())
    frozen=json.loads((ROOT/'results/model.json').read_text());errors=[]
    for name,sha in result['frozen_hashes'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha
    provenance=json.loads((BASE/'event-provenance.json').read_text())
    for index in result['indices']:
        rec=next(x for x in provenance if x['index']==index)
        assert hashlib.sha256((BASE/f'event-{index}.npz').read_bytes()).hexdigest()==rec['cache_sha256']
        assert all(rec['exact_shared_fields'].values()) and rec['all_tutorial_rows_exact_in_both']
        assert all(rec['prior_audit']['shape_checks'].values())
    with np.load(BASE/(args.prefix+'comparison-checkpoints.npz'),allow_pickle=False) as cp:
        r=cp['r'];fs=cp['fs']
        for si,name in enumerate(cp['names'].tolist()):
            for ki,key in enumerate(result['keys']):
                means=cp['all_event_means'][si,:,ki,:]
                mix=(fs[ki,None,:]*means[:,0,None]+(1-fs[ki,None,:])*means[:,1,None])/means[:,1,None]
                selection=frozen['scenarios'][key]['selection']
                xi=(fs[ki]*selection[0]+(1-fs[ki])*selection[1])/selection[1]
                p=np.prod(mix.astype(np.longdouble),axis=0)/xi.astype(np.longdouble)**69
                p=np.asarray(p/p.max(),dtype=float)
                c=cumulative_trapezoid(p,r,initial=0);c/=c[-1]
                q=np.interp([.05,.5,.95],c,r)
                errors.append(float(np.max(abs(q-np.array(result['scenarios'][name]['quantiles'][ki])))))
        assert max(errors)<1e-11
    z=np.linspace(.03,2.5,53);h=1e-5
    numeric=(LAL.luminosity_distance(z+h).value-LAL.luminosity_distance(z-h).value)/(2*h)
    analytic=LAL.luminosity_distance(z).value/(1+z)+(1+z)*LAL.hubble_distance.value/LAL.efunc(z)
    derivative_error=float(np.max(abs(numeric/analytic-1)));assert derivative_error<1e-8
    report=dict(schema='blackhole-multievent-representation-validation/1',status='pass',
        indices=result['indices'],independent_full_product_max_quantile_error=max(errors),
        jacobian_finite_difference_max_relative_error=derivative_error,
        frozen_hashes_verified=result['frozen_hashes'],cache_hashes_and_provenance_flags_verified=True,
        weighted_prior_representation_quantile_error=result['weighted_representation_max_quantile_difference'],
        grid_max_quantile_error=result['max_grid_quantile_error'],
        interpretation='Numerical and provenance checks only; not a calibrated sampling-uncertainty result.')
    (BASE/(args.prefix+'validation.json')).write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
