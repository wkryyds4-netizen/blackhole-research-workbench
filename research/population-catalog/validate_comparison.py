"""Independent source checks, flat-LambdaCDM prior quadrature and full products."""
from pathlib import Path
import hashlib
import json
import h5py
import numpy as np
from scipy.integrate import cumulative_simpson,cumulative_trapezoid,quad
from scipy.special import ndtr,beta

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
SOURCE=['mass_1_source','mass_ratio','redshift','a_1','a_2','cos_tilt_1','cos_tilt_2']
PARAMS=['mass_1']+SOURCE[1:]


def main():
    result=json.loads((BASE/'catalog-comparison.json').read_text());provenance=json.loads((BASE/'event-provenance.json').read_text())
    with np.load(BASE/'catalog-comparison.npz',allow_pickle=False) as cp:
        r=cp['r'];fs=cp['fs'];selection=cp['selection'];means=cp['means'].copy();names=cp['names'].tolist();zg=cp['zgrid'];zdv=cp['zdv'];znorm=float(cp['znorm'])
    frozen=json.loads((ROOT/'results/model.json').read_text());checks=[]
    grid=np.linspace(0,2.3,65537);E=lambda z:np.sqrt(.3065*(1+z)**3+.6935)
    dimensionless_distance=cumulative_simpson(1/E(grid),x=grid,initial=0);DH=299792.458/67.9
    replacements={}
    with np.load(ROOT/'data/posteriors.npz',allow_pickle=False) as tutorial:
        for index in [0,14]:
            info=next(v for v in provenance if v['index']==index);path=BASE/'source-data'/info['source_file']
            assert hashlib.sha256(path.read_bytes()).hexdigest()==info['source_sha256']
            with h5py.File(path,'r') as h:
                ds=h[info['analysis_label']+'/posterior_samples'];x=ds.fields(SOURCE+['phi_12'])[:]
            with np.load(BASE/f'event-{index}.npz',allow_pickle=False) as bank:
                mapping=bank['tutorial_source_rows']
                assert len(np.unique(mapping))==len(mapping)
                for source,k in zip(SOURCE,PARAMS):
                    assert np.array_equal(x[source][mapping],tutorial[f'{index}_{k}'])
                    assert np.array_equal(x[source],bank[k])
                assert np.array_equal(x['phi_12'],bank['phi_12'])
            m=x['mass_1_source'];q=x['mass_ratio'];z=x['redshift'];a1=x['a_1'];a2=x['a_2']
            assert max(z)<grid[-1] and min(a1)>1e-12 and min(a2)>1e-12
            integral=np.interp(z,grid,dimensionless_distance);dl=(1+z)*DH*integral;ddz=DH*(integral+(1+z)/E(z))
            prior=dl**2*ddz*m*(1+z)**2/4
            # Check independent interpolation against direct quadrature on 32 rows.
            ids=np.linspace(0,len(z)-1,32,dtype=int)
            direct=np.array([quad(lambda zz:1/E(zz),0,z[i],epsabs=1e-13,epsrel=1e-13)[0] for i in ids])
            interpolation_error=float(np.max(abs(integral[ids]/direct-1)));assert interpolation_error<1e-8
            pm=.9*m**-2/(1/3-1/100)+.1*np.exp(-.5*((m-35)/5)**2)/(5*np.sqrt(2*np.pi)*(ndtr(13)-ndtr(-6.4)))
            pq=np.where((q>=3/m)&(q<=1)&(m>3),2*q/np.maximum(1-(3/m)**2,1e-12),0)
            pm*=((m>=3)&(m<=100));pz=np.interp(z,zg,zdv)*(1+z)/znorm
            new=[]
            for alpha,b in [(1.5,8.5),(2,5)]:
                magnitude=(a1*a2)**(alpha-1)*((1-a1)*(1-a2))**(b-1)/beta(alpha,b)**2
                base=pm*pq*pz*magnitude/prior
                for sigma in [.3,.6,1.2]:
                    normal=sigma*np.sqrt(2*np.pi)*(ndtr(0)-ndtr(-2/sigma))
                    aligned=np.exp(-.5*(((x['cos_tilt_1']-1)/sigma)**2+((x['cos_tilt_2']-1)/sigma)**2))/normal**2
                    new.append([np.mean(base*aligned),np.mean(base/4)])
            new=np.array(new);ref=np.array(next(e for e in result['event_diagnostics'] if e['index']==index)['new_mean'])
            error=float(np.max(abs(new/ref-1)));assert error<1e-7
            pos=next(i for i,e in enumerate(frozen['events']) if e['index']==index);replacements[index]=(pos,new)
            checks.append(dict(index=index,tutorial_rows_exact=len(mapping),phi12_all_rows_exact=len(x),source_prior_distance_quadrature_error=interpolation_error,new_mean_relative_error=error))
    baseline=means[names.index('five_event_baseline')].copy();case_means={name:baseline.copy() for name in names}
    for index,(pos,new) in replacements.items():
        case_means[f'plus_{index}'][pos]=new;case_means['all_new_joint'][pos]=new
    errors={}
    for name,em in case_means.items():
        qs=[]
        for k in range(6):
            f=fs[k];factors=f[None,:]*em[:,k,0,None]+(1-f[None,:])*em[:,k,1,None]
            likelihood=np.prod(factors/factors[:,500,None],axis=0)
            efficiency=f*selection[k,0]+(1-f)*selection[k,1];likelihood*=(efficiency[500]/efficiency)**len(em)
            cdf=cumulative_trapezoid(likelihood,r,initial=0);cdf/=cdf[-1]
            qs.append(np.interp([.05,.5,.95],cdf,r))
        errors[name]=float(np.max(abs(np.array(qs)-result['scenarios'][name]['quantiles'])));assert errors[name]<1e-8
    report=dict(schema='blackhole-catalog-comparison-independent-validation/1',source_checks=checks,direct_product_quantile_errors=errors,
        shared_inputs='Fixed five-event baseline means, common r grid, selection integrals and specified population cosmology interpolation. Source-prior distances/Jacobian, density formula, source-row equality and full likelihood products independently checked.',
        scope='Two newly restored native-reference posteriors. Does not validate common-reference physics, prior likelihood tails, full catalog restoration, or selection campaign coverage.')
    (BASE/'comparison-validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
