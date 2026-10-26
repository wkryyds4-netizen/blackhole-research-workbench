"""Independent raw-source densities and full-product verification of54 models.

No imports of grid_model, run_grid, model or their weighting/posterior routines.
Uses the explicitly fixed target cosmology grid, while reconstructing PE priors
with an independent flat-LambdaCDM quadrature.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import hashlib
import json
import h5py
import numpy as np
from scipy.integrate import cumulative_simpson,cumulative_trapezoid
from scipy.special import beta,ndtr,betainc
from numpy.polynomial.legendre import leggauss

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
CAT=ROOT/'research/population-catalog'
SEL=ROOT/'research/population-selection-audit'
INV=ROOT/'research/population-invariant'
MIX=ROOT/'research/population-mixture'


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def relative(actual,expected):
    assert np.array_equal(actual==0,expected==0),'Zero-support components differ'
    take=expected!=0
    return float(np.max(abs(actual[take]/expected[take]-1))) if take.any() else 0.


def main():
    result=json.loads((BASE/'results.json').read_text())
    design=json.loads((BASE/'design.json').read_text());specs=result['models']
    gates=design['numerical_gates']
    for path,digest in result['input_sha256'].items():assert sha(ROOT/path)==digest,path
    assert sha(BASE/'checkpoints.npz')==result['checkpoint_sha256']
    with np.load(INV/'control-checkpoints.npz',allow_pickle=False) as cp:zgrid=cp['zgrid'];volume=cp['zdv']
    znorm={k:np.trapezoid(volume*(1+zgrid)**(k-1),zgrid) for k in [0,2,4]}
    assert len(specs)==54 and {r['index'] for r in specs}==set(range(54))
    assert {r['kappa'] for r in specs}=={0,2,4}
    assert {(r['alpha'],r['peak_fraction']) for r in specs}=={(2,.1),(3,.1),(2,0)}
    assert {tuple(r['magnitude']) for r in specs}=={(1.5,8.5),(2,5)}
    assert {r['shape'] for r in specs}=={0,4,16}
    assert len({(r['alpha'],r['peak_fraction'],r['kappa'],*r['magnitude'],r['shape']) for r in specs})==54

    def target(d,prior):
        m=d['m'];q=d['q'];z=d['z'];a1=d['a1'];a2=d['a2']
        x=(d['s1z']+q*d['s2z'])/(1+q)
        massratio=np.where((m>3)&(q>=3/m)&(q<=1),2*q/np.maximum(1-9/(m*m),1e-12),0.)
        gaussian=np.exp(-((m-35)/5)**2/2)/(5*np.sqrt(2*np.pi)*(ndtr(13)-ndtr(-6.4)))
        masses={}
        for alpha,peak in [(2,.1),(3,.1),(2,0)]:
            power=(1-alpha)*m**(-alpha)/(100**(1-alpha)-3**(1-alpha))
            masses[(alpha,peak)]=((1-peak)*power+peak*gaussian)*((m>=3)&(m<=100))
        radial={}
        for aa,bb in [(1.5,8.5),(2.,5.)]:
            radial[(aa,bb)]=(a1*a2)**(aa-1)*((1-a1)*(1-a2))**(bb-1)/beta(aa,bb)**2
        redshifts={k:np.interp(z,zgrid,volume)*(1+z)**(k-1)/znorm[k]*((z>=0)&(z<=2.3)) for k in [0,2,4]}
        even={lam:(1+lam*x*x)/(1+lam*(a1*a1+q*q*a2*a2)/(3*(1+q)**2)) for lam in [0,4,16]}
        negative=np.where(x<0,2.,np.where(x>0,0.,1.));positive=2-negative
        out=[]
        for row in specs:
            base=masses[(row['alpha'],row['peak_fraction'])]*massratio*redshifts[row['kappa']]*radial[tuple(row['magnitude'])]*even[row['shape']]/prior/4
            out.append(np.column_stack([base*negative,base*positive]))
        return np.stack(out,axis=1)

    def event_data(raw):
        z=raw['redshift'];m=raw['mass_1'];end=max(2.3,float(z.max())+.01)
        quadrature=np.linspace(0,end,int(np.ceil(end/(2.3/65536)))+1)
        expansion=np.sqrt(.3065*(1+quadrature)**3+.6935)
        radial=cumulative_simpson(1/expansion,x=quadrature,initial=0)
        dh=299792.458/67.9;comoving=dh*np.interp(z,quadrature,radial)
        luminosity=(1+z)*comoving
        derivative=comoving+(1+z)*dh/np.sqrt(.3065*(1+z)**3+.6935)
        prior=luminosity**2*derivative*m*(1+z)**2/4
        a1=np.sqrt(sum(raw['spin_1'+axis]**2 for axis in 'xyz'))
        a2=np.sqrt(sum(raw['spin_2'+axis]**2 for axis in 'xyz'))
        return dict(m=m,q=raw['mass_ratio'],z=z,a1=a1,a2=a2,s1z=raw['spin_1z'],s2z=raw['spin_2z']),prior

    saved_events=json.loads((BASE/'event-integrals.json').read_text())
    frozen=json.loads((ROOT/'results/model.json').read_text())
    source_records=json.loads((CAT/'event-provenance.json').read_text())
    event_means=[];event_checks=[]
    for position,event in enumerate(frozen['events']):
        index=event['index'];record=next(r for r in saved_events['events'] if r['index']==index)
        source=next(r for r in source_records if r['index']==index)
        path=CAT/f'event-{index}.npz'
        assert sha(path)==record['source_cache_sha256']==source['cache_sha256']
        fields=['mass_1','mass_ratio','redshift']+[f'spin_{i}{axis}' for i in [1,2] for axis in 'xyz']
        with np.load(path,allow_pickle=False) as bank:raw={k:bank[k] for k in fields}
        d,prior=event_data(raw);total=np.zeros((54,2))
        for start in range(0,len(prior),32768):
            sl=slice(start,start+32768)
            total+=target({k:v[sl] for k,v in d.items()},prior[sl]).sum(axis=0)
        mean=total/len(prior)
        error=relative(mean,np.array(record['mean']))
        assert error<gates['independent_event_integral_relative'],f'Event{index} relative difference {error}'
        event_means.append(mean);event_checks.append(dict(index=index,rows=len(prior),maximum_relative_error=error))
        if (position+1)%10==0 or position==68:print(f'Independent event densities: {position+1}/69',flush=True)
    means=np.array(event_means)

    def load_found(handle,mask):
        group=handle['injections'];keys=['mass1_source','mass2_source','redshift','sampling_pdf','mixture_weight']+[f'spin{i}{axis}' for i in [1,2] for axis in 'xyz']
        return {k:group[k][:][mask] for k in keys}

    def cartesian_sum(raw,mask,factor):
        positions=np.flatnonzero(mask);total=np.zeros((54,2))
        for start in range(0,len(positions),32768):
            take=positions[start:start+32768]
            m=raw['mass1_source'][take].astype(float);q=raw['mass2_source'][take].astype(float)/m
            a1=np.sqrt(sum(raw['spin1'+axis][take].astype(float)**2 for axis in 'xyz'))
            a2=np.sqrt(sum(raw['spin2'+axis][take].astype(float)**2 for axis in 'xyz'))
            # Explicit Cartesian -> reduced coordinate proposal Jacobian.
            prior=raw['sampling_pdf'][take].astype(float)*m*(2*np.pi*a1*a1)*(2*np.pi*a2*a2)
            d=dict(m=m,q=q,z=raw['redshift'][take].astype(float),a1=a1,a2=a2,s1z=raw['spin1z'][take].astype(float),s2z=raw['spin2z'][take].astype(float))
            total+=(target(d,prior)*raw['mixture_weight'][take,None,None]*factor).sum(axis=0)
        return total

    oldpath=ROOT/'research/population-stage2/source-data/o1+o2+o3_bbhpop_real+semianalytic-LIGO-T2100377-v2.hdf5'
    assert sha(oldpath)=='81883c571a9fe0b30747c7bc9daca6d358fcd68ef375c7ace0b93aa9ea2e74fe'
    with h5py.File(oldpath,'r') as handle:
        group=handle['injections'];name=group['name'][:]
        union=np.any([group[k][:]>1 for k in group if k.startswith('ifar_')],axis=0)
        early=(name!=b'o3')&(union|(group['optimal_snr_net'][:]>10))
        raw=load_found(handle,early);old_time=float(handle.attrs['analysis_time_s'])
        earlyterm=cartesian_sum(raw,np.ones(early.sum(),bool),1/float(handle.attrs['total_generated']))
    oldsource=json.loads((SEL/'o3-provenance.json').read_text())
    mixture=json.loads((MIX/'download.json').read_text())
    bankpaths={'pure_bbh':SEL/'sources'/oldsource['source']['file'],
               'official_mixture':MIX/'sources'/mixture['file']}
    expected_hashes={'pure_bbh':oldsource['source']['sha256'],'official_mixture':mixture['sha256']}
    selections={}
    for bank,path in bankpaths.items():
        assert sha(path)==expected_hashes[bank]
        with h5py.File(path,'r') as handle:
            group=handle['injections'];searches=[k for k in group if k.startswith('ifar_')]
            assert len(searches)==5
            union=np.any([group[k][:]>1 for k in searches],axis=0)
            matched=np.any([group[k][:]>1 for k in searches if k!='ifar_cwb'],axis=0)
            raw=load_found(handle,union)
            factor=float(handle.attrs['analysis_time_s'])/old_time/float(handle.attrs['total_generated'])
        for policy,mask in [('matched_filter',matched[union]),('with_cwb',np.ones(union.sum(),bool))]:
            selections[bank+':'+policy]=earlyterm+cartesian_sum(raw,mask,factor)

    checks={};peak_max=0.;quantile_max=0.;selection_max=0.
    with np.load(BASE/'checkpoints.npz',allow_pickle=False) as cp:
        assert relative(means,cp['means'])<gates['independent_event_integral_relative']
        g=cp['g'];mid=len(g)//2
        for name,selection in selections.items():
            se=relative(selection,cp[name+'_selection']);assert se<gates['independent_selection_relative'],(name,se)
            # Multiplicative69-factor posterior, normalized after relative
            # scaling at a fixed g=1/2. Not the producer's log-likelihood sum.
            likelihood=np.ones((54,len(g)))
            for event in means:
                mix=g*event[:,0,None]+(1-g)*event[:,1,None]
                likelihood*=mix/mix[:,mid,None]
            efficiency=g*selection[:,0,None]+(1-g)*selection[:,1,None]
            likelihood*=(efficiency[:,mid,None]/efficiency)**69
            likelihood/=np.trapezoid(likelihood,g,axis=1)[:,None]
            pe=float(np.max(abs(likelihood-cp[name+'_density']))/likelihood.max())
            assert pe<gates['independent_peak_scaled_density'],(name,pe)
            cdf=cumulative_trapezoid(likelihood,g,axis=1,initial=0);cdf/=cdf[:,-1,None]
            q=np.array([np.interp([.05,.5,.95],c,g) for c in cdf])
            qe=float(np.max(abs(q-result['cases'][name]['quantiles'])));assert qe<gates['independent_quantile_absolute'],(name,qe)
            checks[name]=dict(selection_relative_difference=se,peak_scaled_density_difference=pe,maximum_quantile_difference=qe)
            peak_max=max(peak_max,pe);quantile_max=max(quantile_max,qe);selection_max=max(selection_max,se)

    # Fixed-order Gaussian integration on each interpolation segment, separate
    # from the producer's adaptive quadrature with explicit breakpoints.
    nodes,gauss=leggauss(8)
    def integral(low,high,kappa):
        edges=np.r_[low,zgrid[(zgrid>low)&(zgrid<high)],high]
        x=(edges[1:,None]+edges[:-1,None])/2+(edges[1:,None]-edges[:-1,None])*nodes/2
        f=np.interp(x,zgrid,volume)*(1+x)**(kappa-1)
        return float(np.sum((edges[1:]-edges[:-1])/2*(f@gauss)))
    x,w=leggauss(128);m=95+5*x
    support_errors=[]
    for row,wanted in zip(specs,result['support']['models']):
        k=row['kappa'];tail=integral(1.9,2.3,k)/integral(0,2.3,k)
        alpha=row['alpha'];peak=row['peak_fraction']
        pdf=(1-peak)*(1-alpha)*m**(-alpha)/(100**(1-alpha)-3**(1-alpha))
        pdf+=peak*np.exp(-((m-35)/5)**2/2)/(5*np.sqrt(2*np.pi)*(ndtr(13)-ndtr(-6.4)))
        added=float(5*np.dot(w,pdf*(1-100/(m*m))/(1-9/(m*m))))
        aa,bb=row['magnitude'];spin=betainc(aa,bb,.998)**2
        for bank,missing in [('pure_bbh',1-spin*(1-tail)),('official_mixture',1-spin*(1-tail+tail*added))]:
            support_errors.append(abs(missing-wanted['banks'][bank]['uncovered_target_probability']))
    assert max(support_errors)<gates['independent_support_absolute'],max(support_errors)
    validation=dict(schema='blackhole-full69-mass-redshift-independent-validation/1',status='pass',
                    result_sha256=sha(BASE/'results.json'),validator_sha256=sha(Path(__file__)),
                    event_count=69,event_rows=sum(r['rows'] for r in event_checks),model_count=54,conditional_fits=216,
                    event_checks=event_checks,posterior_checks=checks,
                    maximum_event_integral_relative_difference=max(r['maximum_relative_error'] for r in event_checks),
                    maximum_selection_relative_difference=selection_max,maximum_peak_scaled_density_difference=peak_max,
                    maximum_quantile_difference=quantile_max,maximum_support_absolute_difference=max(support_errors),
                    independence='Independent raw-PE flat-LambdaCDM prior integral, Cartesian spin reconstruction,54 target density formulas, original HDF early/O3 sums and69 multiplicative event factors. Saved target cosmology grid intentionally reused.',
                    limits='Implementation consistency on the fixed public inputs. Does not validate omitted likelihood/sensitivity, chain correlations, actual campaign coverage, general parameter correlations or physical channels.')
    (BASE/'validation.json').write_text(json.dumps(validation,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in validation.items() if k not in ['event_checks','posterior_checks']},indent=2))


if __name__=='__main__':
    try:main()
    except Exception as error:
        path=BASE/'validation-failures.json'
        history=json.loads(path.read_text()) if path.exists() else []
        history.append(dict(result_sha256=sha(BASE/'results.json') if (BASE/'results.json').exists() else None,
                            validator_sha256=sha(Path(__file__)),error_type=type(error).__name__,message=str(error)))
        path.write_text(json.dumps(history,indent=2)+'\n',encoding='utf-8')
        raise
