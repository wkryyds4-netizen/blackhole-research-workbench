"""Independent Cartesian densities, PE-prior quadrature and full products."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.integrate import cumulative_simpson,cumulative_trapezoid
from scipy.special import beta,ndtr

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
CAT=ROOT/'research/population-catalog'
SEL=ROOT/'research/population-selection-audit'


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    report=json.loads((BASE/'control-results.json').read_text());saved=json.loads((BASE/'event-integrals.json').read_text())
    frozen=json.loads((ROOT/'results/model.json').read_text())
    with np.load(BASE/'control-checkpoints.npz',allow_pickle=False) as cp:
        g=cp['g'];expected=cp['means'];zg=cp['zgrid'];dv=cp['zdv'];znorm=float(cp['znorm'])
        early_reference=cp['early_weights'];full_reference=cp['o3_union_weights'];mask_reference=cp['o3_matched_mask'];factor=float(cp['o3_factor'])
    DH=299792.458/67.9
    def target(d,prior):
        m=d['m'];q=d['q'];z=d['z'];a1=d['a1'];a2=d['a2'];x=(d['s1z']+q*d['s2z'])/(1+q)
        mass=.9*m**-2/(1/3-1/100)+.1*np.exp(-.5*((m-35)/5)**2)/(5*np.sqrt(2*np.pi)*(ndtr(13)-ndtr(-6.4)))
        mass*=((m>=3)&(m<=100))
        ratio=np.where((m>3)&(q>=3/m)&(q<=1),2*q/np.maximum(1-(3/m)**2,1e-12),0)
        rz=np.interp(z,zg,dv)*(1+z)/znorm*((z>=0)&(z<=2.3))
        base=mass*ratio*rz/prior;mom=(a1*a1+q*q*a2*a2)/(3*(1+q)**2);out=[]
        negative=np.where(x<0,2.,np.where(x>0,0.,1.));positive=2-negative
        for aa,bb in [(1.5,8.5),(2.,5.)]:
            radial=(a1*a2)**(aa-1)*((1-a1)*(1-a2))**(bb-1)/beta(aa,bb)**2
            for lam in [0.,4.,16.]:
                even=(1+lam*x*x)/(1+lam*mom)
                w=base*radial*even/4
                out.append(np.column_stack([w*negative,w*positive]))
        return np.stack(out,axis=1)
    def relative(actual,wanted):
        assert np.array_equal(actual==0,wanted==0),'Zero-support components differ'
        mask=wanted!=0
        return float(np.max(abs(actual[mask]/wanted[mask]-1))) if mask.any() else 0.
    checks=[];means=[]
    for event in frozen['events']:
        ix=event['index'];rec=next(r for r in saved['events'] if r['index']==ix);path=CAT/f'event-{ix}.npz'
        assert sha(path)==rec['source_cache_sha256']
        with np.load(path,allow_pickle=False) as bank:a={k:bank[k] for k in bank.files}
        z=a['redshift'];end=max(2.3,float(z.max())+.01)
        grid=np.linspace(0,end,int(np.ceil(end/(2.3/65536)))+1)
        integral=cumulative_simpson(1/np.sqrt(.3065*(1+grid)**3+.6935),x=grid,initial=0)
        comoving=DH*np.interp(z,grid,integral);dl=(1+z)*comoving
        derivative=comoving+(1+z)*DH/np.sqrt(.3065*(1+z)**3+.6935)
        prior=dl**2*derivative*a['mass_1']*(1+z)**2/4
        a1=np.sqrt(a['spin_1x']**2+a['spin_1y']**2+a['spin_1z']**2)
        a2=np.sqrt(a['spin_2x']**2+a['spin_2y']**2+a['spin_2z']**2)
        d=dict(m=a['mass_1'],q=a['mass_ratio'],z=z,a1=a1,a2=a2,s1z=a['spin_1z'],s2z=a['spin_2z'])
        w=target(d,prior);mean=w.mean(axis=0);error=relative(mean,np.array(rec['mean']));assert error<2e-7
        means.append(mean);checks.append(dict(index=ix,rows=len(z),mean_relative_error=error,source_rows_outside_target_z=int(np.count_nonzero(z>2.3))))
    means=np.array(means);assert relative(means,expected)<2e-7
    provenance=json.loads((SEL/'o3-provenance.json').read_text());assert sha(SEL/'o3-bank.npz')==provenance['cache_sha256']
    with np.load(SEL/'o3-bank.npz',allow_pickle=False) as bank:
        raw={k[5:]:bank[k] for k in bank.files if k.startswith('full_')}
        old={k[9:]:bank[k] for k in bank.files if k.startswith('combined_') and k!='combined_o3_to_full'}
    names=sorted(k for k in raw if k.startswith('ifar_'));assert len(names)==5
    union=np.any([raw[k]>1 for k in names],axis=0);mf=np.any([raw[k]>1 for k in names if k!='ifar_cwb'],axis=0)
    early=(old['name']!=b'o3')&(np.any([old[k]>1 for k in names],axis=0)|(old['optimal_snr_net']>10))
    def cartesian_weights(raw,mask):
        m=raw['mass1_source'][mask].astype(float);n=raw['mass2_source'][mask].astype(float)
        a1=np.sqrt(sum(raw[f'spin1{ax}'][mask].astype(float)**2 for ax in 'xyz'))
        a2=np.sqrt(sum(raw[f'spin2{ax}'][mask].astype(float)**2 for ax in 'xyz'))
        d=dict(m=m,q=n/m,z=raw['redshift'][mask].astype(float),a1=a1,a2=a2,
            s1z=raw['spin1z'][mask].astype(float),s2z=raw['spin2z'][mask].astype(float))
        # Cartesian-to-(m1,q,a1,c1,a2,c2) Jacobian with azimuths integrated.
        prior=raw['sampling_pdf'][mask].astype(float)*m*(2*np.pi*a1*a1)*(2*np.pi*a2*a2)
        return target(d,prior)
    early_w=cartesian_weights(old,early)*old['mixture_weight'][early,None,None]/provenance['combined_attributes']['total_generated']
    full_w=cartesian_weights(raw,union)
    early_error=relative(early_w,early_reference);full_error=relative(full_w,full_reference)
    assert max(early_error,full_error)<1e-9 and np.array_equal(mf[union],mask_reference)
    selections=dict(matched_filter=early_w.sum(axis=0)+full_w[mf[union]].sum(axis=0)*factor,
        with_cwb=early_w.sum(axis=0)+full_w.sum(axis=0)*factor,without_selection_diagnostic=np.ones((6,2)))
    product_checks={};middle=len(g)//2
    for name,selection in selections.items():
        qs=[]
        for k in range(6):
            terms=g[None,:]*means[:,k,0,None]+(1-g[None,:])*means[:,k,1,None]
            terms=terms/terms[:,middle,None]
            p=np.prod(terms.astype(np.longdouble),axis=0)
            xi=g*selection[k,0]+(1-g)*selection[k,1]
            p*=(xi[middle]/xi).astype(np.longdouble)**69
            p=np.asarray(p/p.max(),dtype=float)
            c=cumulative_trapezoid(p,g,initial=0);c/=c[-1];qs.append(np.interp([.05,.5,.95],c,g))
        error=float(np.max(abs(np.array(qs)-report['scenarios'][name]['quantiles'])));assert error<2e-8
        product_checks[name]=error
    output=dict(schema='blackhole-invariant-control-independent-validation/1',events=len(checks),rows=sum(r['rows'] for r in checks),checks=checks,
        maximum_event_mean_relative_error=max(r['mean_relative_error'] for r in checks),
        early_cartesian_weight_relative_error=early_error,full_o3_cartesian_weight_relative_error=full_error,
        full_product_quantile_errors=product_checks,result_sha256=sha(BASE/'control-results.json'),
        scope='Independent Cartesian spin inputs, explicit PE-prior distance quadrature, target formulas and full69 products. Shares specified cosmology grid, target model and injection masks. Does not prove preserved PN measure or physical campaign coverage.')
    (BASE/'control-validation.json').write_text(json.dumps(output,indent=2),encoding='utf-8')
    print(json.dumps({k:output[k] for k in ['events','rows','maximum_event_mean_relative_error','early_cartesian_weight_relative_error','full_o3_cartesian_weight_relative_error','full_product_quantile_errors']},indent=2))


if __name__=='__main__':main()
