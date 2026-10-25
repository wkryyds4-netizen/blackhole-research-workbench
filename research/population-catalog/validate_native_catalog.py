"""Independent vector-density and full-product verification for restored files."""
from pathlib import Path
import json
import hashlib
import numpy as np
from scipy.integrate import cumulative_simpson,cumulative_trapezoid
from scipy.special import beta,ndtr

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]


def main():
    result=json.loads((BASE/'native-catalog-results.json').read_text());integrals=json.loads((BASE/'native-event-integrals.json').read_text())
    frozen=json.loads((ROOT/'results/model.json').read_text())
    with np.load(BASE/'native-catalog-checkpoints.npz',allow_pickle=False) as cp:
        means=cp['baseline'].copy();selection=cp['selection'];r=cp['r'];fs=cp['fs'];zg=cp['zgrid'];zdv=cp['zdv'];znorm=float(cp['znorm'])
    grid=np.linspace(0,2.3,65537);E=np.sqrt(.3065*(1+grid)**3+.6935);I=cumulative_simpson(1/E,x=grid,initial=0);DH=299792.458/67.9
    checks=[]
    for index in result['newly_recomputed_source_indices']:
        rec=next(v for v in integrals['events'] if v['index']==index);path=BASE/f'event-{index}.npz'
        assert hashlib.sha256(path.read_bytes()).hexdigest()==rec['cache_sha256']
        with np.load(path,allow_pickle=False) as bank:d={k:bank[k] for k in bank.files}
        m=d['mass_1'];q=d['mass_ratio'];z=d['redshift'];a1=d['a_1'];a2=d['a_2']
        assert min(z)>0 and np.isfinite(z).all() and min(a1)>1e-12 and min(a2)>1e-12
        # The PE prior must be evaluated on every source row, including rows
        # outside the target population support. Those rows stay in the MC
        # denominator and receive zero target density, rather than being cut.
        if max(z)>=grid[-1]:
            endpoint=float(max(z))+.01
            points=int(np.ceil(endpoint/(2.3/65536)))+1
            grid=np.linspace(0,endpoint,points)
            E=np.sqrt(.3065*(1+grid)**3+.6935)
            I=cumulative_simpson(1/E,x=grid,initial=0)
        assert np.max(abs(d['mass_1_detector']/(1+z)/m-1))<1e-12
        assert np.max(abs(d['mass_2_source']/m/q-1))<1e-12
        phi=np.arctan2(d['spin_2y'],d['spin_2x'])-np.arctan2(d['spin_1y'],d['spin_1x'])
        phase_error=float(np.max(abs(np.cos(phi)-np.cos(d['phi_12']))));assert phase_error<1e-9
        effective=(a1*d['cos_tilt_1']+q*a2*d['cos_tilt_2'])/(1+q)
        effective_error=float(np.max(abs(effective-d['chi_eff'])));assert effective_error<1e-12
        distance=DH*np.interp(z,grid,I);dl=(1+z)*distance
        derivative=distance+(1+z)*DH/np.sqrt(.3065*(1+z)**3+.6935)
        prior=dl*dl*derivative*m*(1+z)**2/4
        mass=.9*m**-2/(1/3-1/100)+.1*np.exp(-.5*((m-35)/5)**2)/(5*np.sqrt(2*np.pi)*(ndtr(13)-ndtr(-6.4)))
        mass*=((m>=3)&(m<=100));ratio=np.where((q>=3/m)&(q<=1)&(m>3),2*q/np.maximum(1-(3/m)**2,1e-12),0)
        target_support=(z>=0)&(z<=2.3)
        redshift=np.interp(z,zg,zdv)*(1+z)/znorm*target_support;new=[]
        for aa,bb in [(1.5,8.5),(2,5)]:
            radial=(a1*a2)**(aa-1)*((1-a1)*(1-a2))**(bb-1)/beta(aa,bb)**2
            base=mass*ratio*redshift*radial/prior
            assert np.all(base[~target_support]==0)
            for sigma in [.3,.6,1.2]:
                norm=sigma*np.sqrt(2*np.pi)*(ndtr(0)-ndtr(-2/sigma))
                aligned=np.exp(-.5*((d['cos_tilt_1']-1)**2+(d['cos_tilt_2']-1)**2)/sigma**2)/norm**2
                new.append([float(np.mean(base*aligned)),float(np.mean(base/4))])
        new=np.array(new);error=float(np.max(abs(new/rec['mean']-1)));assert error<2e-7
        pos=next(i for i,e in enumerate(frozen['events']) if e['index']==index);means[pos]=new
        checks.append(dict(index=index,rows=len(m),mean_relative_error=error,azimuth_cosine_error=phase_error,chi_eff_error=effective_error,
            source_redshift_range=[float(min(z)),float(max(z))],prior_quadrature_maximum_redshift=float(grid[-1]),
            rows_outside_target_redshift_support=int(np.count_nonzero(~target_support)),outside_target_rows_retained_with_zero_target_density=True))
    quantiles=[]
    for k in range(6):
        f=fs[k];factors=f[None,:]*means[:,k,0,None]+(1-f[None,:])*means[:,k,1,None]
        p=np.prod(factors/factors[:,500,None],axis=0)
        efficiency=f*selection[k,0]+(1-f)*selection[k,1];p*=(efficiency[500]/efficiency)**len(means)
        c=cumulative_trapezoid(p,r,initial=0);c/=c[-1];quantiles.append(np.interp([.05,.5,.95],c,r))
    error=float(np.max(abs(np.array(quantiles)-result['quantiles'])));assert error<2e-8
    report=dict(schema='blackhole-native-catalog-independent-validation/2',new_source_events_checked=len(checks),checks=checks,full_product_quantile_error=error,
        result_sha256=hashlib.sha256((BASE/'native-catalog-results.json').read_bytes()).hexdigest(),
        scope='Independent source-coordinate density and flat-LambdaCDM prior quadrature, mass/azimuth/effective-spin identities, and full69-event products. Shares specified target cosmology grid, fixed selection and previous five-event baseline. Does not establish common-reference physics or likelihood beyond prior support.')
    (BASE/'native-catalog-validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(dict(events=len(checks),rows=sum(r['rows'] for r in checks),maximum_mean_relative_error=max(r['mean_relative_error'] for r in checks),full_product_quantile_error=error),indent=2))


if __name__=='__main__':main()
