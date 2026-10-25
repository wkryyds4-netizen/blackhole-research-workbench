"""Independent Cartesian-density summation directly from the official HDF.

Does not import the production physical(), weights(), or posterior() helpers.
Shares only the specified Planck15 interpolation table and fixed event means.
Run with the h5py adapter environment.
"""
from pathlib import Path
import hashlib
import json
import h5py
import numpy as np
from scipy.special import beta,ndtr
from scipy.integrate import cumulative_trapezoid

BASE=Path(__file__).resolve().parent


def main():
    prov=json.loads((BASE/'cumulative-provenance.json').read_text())
    result=json.loads((BASE/'cumulative-results.json').read_text())
    path=BASE/'sources'/prov['source_file']
    assert hashlib.sha256(path.read_bytes()).hexdigest()==prov['source_sha256']
    with h5py.File(path,'r') as source:
        x=source['events'][:];ng=source.attrs['total_generated']
    with np.load(BASE/'cumulative-checkpoints.npz',allow_pickle=False) as cp:
        zgrid=cp['zgrid'];zdv=cp['zdv'];znorm=float(cp['znorm']);r=cp['r'];fs=cp['fs'];means=cp['means']
    m=x['mass1_source'];m2=x['mass2_source'];q=m2/m;z=x['redshift']
    a=np.array([np.sqrt(x[f'spin{i}x']**2+x[f'spin{i}y']**2+x[f'spin{i}z']**2) for i in [1,2]])
    c=np.array([x[f'spin{i}z']/a[i-1] for i in [1,2]])
    power=m**-2/(1/3-1/100)
    peak=np.exp(-.5*((m-35)/5)**2)/(5*np.sqrt(2*np.pi)*(ndtr(13)-ndtr(-6.4)))
    # In dm1 dm2, pq(q|m1) dq/dm2 = 2 m2/(m1^2 - 3^2).
    mass=np.where((m>=3)&(m<=100)&(m2>=3)&(m2<=m),(.9*power+.1*peak)*2*m2/np.maximum(m*m-9,1e-12),0)
    redshift=np.interp(z,zgrid,zdv)*(1+z)/znorm*((z>=0)&(z<=2.3))
    prefactor=mass*redshift*x['weights']*np.exp(-x['lnpdraw_mass1_source_mass2_source_redshift_spin1x_spin1y_spin1z_spin2x_spin2y_spin2z'])/ng
    early=x['semianalytic_observed_phase_maximized_snr_net']>0
    real=np.logical_or.reduce([x[f'o3_{p}_far']<1 for p in ['gstlal','mbta','pycbc_bbh','pycbc_hyperbank']])
    assert not real[early].any()
    masks={'cumulative_optimal':real|(early&(x['estimated_optimal_snr_net']>10)),
           'cumulative_observed':real|(x['semianalytic_observed_phase_maximized_snr_net']>10)}
    sums={name:[] for name in masks}
    for alpha,b in [(1.5,8.5),(2,5)]:
        # p(a,c,phi)/(a^2) for each spin, uniform azimuth p(phi)=1/(2pi).
        radial=np.prod(a**(alpha-1)*(1-a)**(b-1)/beta(alpha,b)/(2*np.pi*a*a),axis=0)
        for sigma in [.3,.6,1.2]:
            aligned=np.prod(np.exp(-.5*((c-1)/sigma)**2)/(sigma*np.sqrt(2*np.pi)*(ndtr(0)-ndtr(-2/sigma))),axis=0)
            cart=radial*prefactor
            for name,mask in masks.items():sums[name].append([np.sum(cart[mask]*aligned[mask]),np.sum(cart[mask]/4)])
    errors={};qerrors={}
    for name,rows in sums.items():
        selection=np.array(rows);ref=np.array(result['selection_components'][name])
        errors[name]=float(np.max(abs(selection/ref-1)));assert errors[name]<1e-10
        qs=[]
        for k in range(6):
            f=fs[k]
            # Direct product of scaled event factors, separately from selection.
            event_factors=f[None,:]*means[:,k,0,None]+(1-f[None,:])*means[:,k,1,None]
            likelihood=np.prod(event_factors/event_factors[:,500,None],axis=0)
            efficiency=f*selection[k,0]+(1-f)*selection[k,1]
            likelihood*= (efficiency[500]/efficiency)**len(means)
            cdf=cumulative_trapezoid(likelihood,r,initial=0);cdf/=cdf[-1]
            qs.append(np.interp([.05,.5,.95],cdf,r))
        qerrors[name]=float(np.max(abs(np.array(qs)-result['scenarios'][name]['quantiles'])));assert qerrors[name]<1e-11
    report=dict(schema='blackhole-cumulative-independent-validation/1',source_sha256=prov['source_sha256'],
        direct_cartesian_selection_relative_errors=errors,direct_product_quantile_errors=qerrors,
        shared_inputs='Fixed event means, common r grid and Planck15 cosmology interpolation table. Density/Jacobian formula and likelihood product independently implemented.',
        not_tested='Astrophysical adequacy, source campaign sampling design, hidden support and common spin reference.')
    (BASE/'cumulative-validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
