"""Independent HDF-to-Cartesian-density integral and full-product check.

No imports from model.py, study_mixture.py or their weighting/posterior helpers.
Reuses the explicitly unchanged event-integral inputs and saved cosmology grid.
"""
from pathlib import Path
import hashlib
import json
import h5py
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.special import beta, ndtr
from numpy.polynomial.legendre import leggauss

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
INV=ROOT/'research/population-invariant'


def sha(path):
    with path.open('rb') as stream: return hashlib.file_digest(stream,'sha256').hexdigest()


def relative(a,b):
    assert np.array_equal(a==0,b==0)
    take=b!=0
    return float(np.max(np.abs(a[take]/b[take]-1))) if take.any() else 0.


def main():
    result=json.loads((BASE/'results.json').read_text(encoding='utf-8'))
    source=json.loads((BASE/'download.json').read_text(encoding='utf-8'))
    for name,digest in result['input_sha256'].items(): assert sha(ROOT/name)==digest, name
    assert sha(BASE/'sources'/source['file'])==source['sha256']==result['source_sha256']
    assert sha(BASE/'study-checkpoints.npz')==result['checkpoint_sha256']
    with np.load(INV/'control-checkpoints.npz',allow_pickle=False) as cp:
        zgrid=cp['zgrid']; volume=cp['zdv']; znorm=float(cp['znorm'])

    def cartesian(raw,metric):
        m=np.asarray(raw['mass1_source'],float); n=np.asarray(raw['mass2_source'],float)
        q=n/m; z=np.asarray(raw['redshift'],float)
        a1=np.sqrt(sum(np.asarray(raw['spin1'+axis],float)**2 for axis in 'xyz'))
        a2=np.sqrt(sum(np.asarray(raw['spin2'+axis],float)**2 for axis in 'xyz'))
        s1z=np.asarray(raw['spin1z'],float); s2z=np.asarray(raw['spin2z'],float)
        c1=s1z/a1; c2=s2z/a2
        gaussian=np.exp(-((m-35)/5)**2/2)/(5*np.sqrt(2*np.pi)*(ndtr(13)-ndtr(-6.4)))
        pm=(.9/(m*m*(1/3-1/100))+.1*gaussian)*((m>=3)&(m<=100))
        pq=2*q/np.maximum(1-9/(m*m),1e-12)*((m>3)&(n>=3)&(q<=1))
        pz=np.interp(z,zgrid,volume)*(1+z)/znorm*((z>=0)&(z<=2.3))
        # Target density in d m1 d m2 d^3 S1 d^3 S2 d z, with
        # azimuthal measures explicit. Joint proposal remains joint.
        jacobian=1/(m*(2*np.pi*a1*a1)*(2*np.pi*a2*a2))
        common=pm*pq*pz*jacobian/raw['sampling_pdf']
        x=(s1z+q*s2z)/(1+q)
        output=[]
        for aa,bb in [(1.5,8.5),(2.,5.)]:
            radial=(a1*a2)**(aa-1)*((1-a1)*(1-a2))**(bb-1)/(beta(aa,bb)**2)
            if metric=='g':
                for lam in [0.,4.,16.]:
                    tilt=(1+lam*x*x)/(1+lam*(a1*a1+q*q*a2*a2)/(3*(1+q)**2))/4
                    negative=np.where(x<0,2.,np.where(x>0,0.,1.))
                    v=common*radial*tilt
                    output.append(np.column_stack([v*negative,v*(2-negative)]))
            else:
                for sigma in [.3,.6,1.2]:
                    norm=sigma*np.sqrt(2*np.pi)*(ndtr(0)-ndtr(-2/sigma))
                    aligned=np.exp(-((c1-1)**2+(c2-1)**2)/(2*sigma*sigma))/(norm*norm)
                    v=common*radial
                    output.append(np.column_stack([v*aligned,v/4]))
        return np.stack(output,axis=1), dict(m=m,n=n,z=z,a1=a1,a2=a2)

    oldpath=ROOT/'research/population-stage2/source-data/o1+o2+o3_bbhpop_real+semianalytic-LIGO-T2100377-v2.hdf5'
    assert sha(oldpath)=='81883c571a9fe0b30747c7bc9daca6d358fcd68ef375c7ace0b93aa9ea2e74fe'
    def load(handle,mask):
        group=handle['injections']
        keys=['mass1_source','mass2_source','redshift','sampling_pdf','mixture_weight']+[f'spin{i}{ax}' for i in [1,2] for ax in 'xyz']
        return {k:group[k][:][mask] for k in keys}
    with h5py.File(oldpath,'r') as handle:
        group=handle['injections']; names=group['name'][:]
        detected=np.any([group[k][:]>1 for k in group if k.startswith('ifar_')],axis=0)
        early=(names!=b'o3')&(detected|(group['optimal_snr_net'][:]>10))
        earlyraw=load(handle,early)
        old_time=float(handle.attrs['analysis_time_s']); old_total=float(handle.attrs['total_generated'])
    with h5py.File(BASE/'sources'/source['file'],'r') as handle:
        group=handle['injections']; fields=[k for k in group if k.startswith('ifar_')]
        union=np.any([group[k][:]>1 for k in fields],axis=0)
        matched=np.any([group[k][:]>1 for k in fields if k!='ifar_cwb'],axis=0)
        mixed=load(handle,union); policy_masks={'matched_filter':matched[union],'with_cwb':np.ones(union.sum(),bool)}
        factor=float(handle.attrs['analysis_time_s'])/old_time/float(handle.attrs['total_generated'])
    checks={}; maximum_selection=0.;maximum_density=0.;maximum_quantile=0.
    with np.load(BASE/'study-checkpoints.npz',allow_pickle=False) as saved:
        for metric in ['g','r']:
            earlyweights,_=cartesian(earlyraw,metric)
            earlyterm=(earlyweights*earlyraw['mixture_weight'][:,None,None]/old_total).sum(axis=0)
            weights,coords=cartesian(mixed,metric)
            weights*=mixed['mixture_weight'][:,None,None]*factor
            old_domain=(coords['m']>=2)&(coords['m']<=100)&(coords['n']>=2)&(coords['z']<=1.9)&(coords['a1']<=.998)&(coords['a2']<=.998)
            terms=saved[metric+'_means']; grid=saved[metric+'_grid']
            if metric=='g': coefficients=np.broadcast_to(grid,(6,len(grid)))
            else:
                coefficients=np.array([(.5-grid)/(.5-(ndtr(-1/s)-ndtr(-2/s))/(ndtr(0)-ndtr(-2/s))) for s in [.3,.6,1.2]*2])
            early_error=relative(earlyterm,saved[metric+'_early']); assert early_error<1e-9
            checks[metric]={}
            for policy,mask in policy_masks.items():
                scenarios={'mixture_full_union':earlyterm+weights[mask].sum(axis=0),
                           'mixture_bbh_domain_only':earlyterm+weights[mask&old_domain].sum(axis=0)}
                for name,selection in scenarios.items():
                    expected=saved[f'{metric}_{policy}_{name}_selection']
                    se=relative(selection,expected);assert se<1e-9
                    # All 69 event factors, accumulated multiplicatively after
                    # per-event scaling; independent of log-product producer.
                    likelihood=np.ones((6,len(grid)))
                    for event in terms:
                        value=coefficients*event[:,0,None]+(1-coefficients)*event[:,1,None]
                        likelihood*=value/value.max(axis=1)[:,None]
                    efficiency=coefficients*selection[:,0,None]+(1-coefficients)*selection[:,1,None]
                    likelihood*=(efficiency.min(axis=1)[:,None]/efficiency)**69
                    likelihood/=np.trapezoid(likelihood,grid,axis=1)[:,None]
                    density_error=float(np.max(abs(likelihood-saved[f'{metric}_{policy}_{name}_density']))/likelihood.max())
                    assert density_error<1e-9
                    cdf=cumulative_trapezoid(likelihood,grid,axis=1,initial=0)
                    cdf/=cdf[:,-1,None]
                    quantiles=np.array([np.interp([.05,.5,.95],row,grid) for row in cdf])
                    qe=float(np.max(abs(quantiles-np.array(result['metrics'][metric]['scenarios'][policy][name]['quantiles']))));assert qe<2e-6
                    baseline=saved[f'{metric}_{policy}_pure_bbh_baseline_density']
                    tv=.5*np.trapezoid(abs(likelihood-baseline),grid,axis=1)
                    assert np.max(abs(tv-result['metrics'][metric]['scenarios'][policy][name]['posterior_total_variation_from_pure_bbh']))<1e-9
                    checks[metric][policy+'_'+name]=dict(early_relative_error=early_error,selection_relative_error=se,
                                                          maximum_relative_peak_density_error=density_error,maximum_quantile_error=qe)
                    maximum_selection=max(maximum_selection,se,early_error);maximum_density=max(maximum_density,density_error);maximum_quantile=max(maximum_quantile,qe)
    # Piecewise polynomial integration of the *coded* interpolated cosmology,
    # rather than the producer's adaptive quadrature.
    def zarea(low,high):
        z=np.r_[low,zgrid[(zgrid>low)&(zgrid<high)],high]
        dv=np.interp(z,zgrid,volume);h=np.diff(z);slope=np.diff(dv)/h
        return float(np.sum(dv[:-1]*(1+z[:-1])*h+(dv[:-1]+slope*(1+z[:-1]))*h*h/2+slope*h**3/3))
    tail=zarea(1.9,2.3)/znorm
    nodes,weights=leggauss(128);m=95+5*nodes
    mass=.9*m**-2/(1/3-1/100)+.1*np.exp(-.5*((m-35)/5)**2)/(5*np.sqrt(2*np.pi)*(ndtr(13)-ndtr(-6.4)))
    added=float(5*np.dot(weights,mass*(1-100/(m*m))/(1-9/(m*m))))
    support_error=max(abs(tail-result['support']['target_redshift_tail_mass']),abs(added-result['support']['conditional_mass_fraction_covered_by_imbh_at_high_z']))
    assert support_error<2e-8
    provenance=json.loads((BASE/'provenance.json').read_text())
    assert sha(BASE/'mixture-membership.npz')==provenance['mapping_sha256']
    old_provenance=json.loads((ROOT/'research/population-selection-audit/o3-provenance.json').read_text())
    pure_path=ROOT/'research/population-selection-audit/sources'/old_provenance['source']['file']
    assert sha(pure_path)==old_provenance['source']['sha256']
    with np.load(BASE/'mixture-membership.npz',allow_pickle=False) as mapping:join=mapping['pure_bbh_to_mixture']
    with h5py.File(pure_path,'r') as pure,h5py.File(BASE/'sources'/source['file'],'r') as mixture:
        take=(pure['injections/mass1_source'][:]<90)&(pure['injections/mass2_source'][:]>=3)
        ratio=mixture['injections/sampling_pdf'][:][join][take]/pure['injections/sampling_pdf'][:][take]
        ratio/=float(pure.attrs['total_generated'])/float(mixture.attrs['total_generated'])
        assert int(take.sum())==result['joint_proposal_bridge']['single_population_overlap_rows']
        bridge_error=float(np.max(abs(np.quantile(ratio,[0,.5,.99,1])-result['joint_proposal_bridge']['joint_pdf_ratio_0_50_99_100'])))
        assert bridge_error<1e-12
    report=dict(schema='blackhole-official-o3-mixture-independent-validation/1',status='pass',
                result_sha256=sha(BASE/'results.json'),validator_sha256=sha(Path(__file__)),checks=checks,
                maximum_selection_relative_error=maximum_selection,maximum_relative_peak_density_error=maximum_density,
                maximum_quantile_error=maximum_quantile,support_quadrature_absolute_error=support_error,
                joint_pdf_ratio_summary_error=bridge_error,
                independence='Independent Cartesian target/proposal ratios, direct HDF early/O3 sums, full 69-event multiplicative posterior and piecewise/Gauss-Legendre support integrals. Fixed event-integral source inputs reused.',
                scope='Implementation and saved-input verification. Not independent data, historical protocol recovery or calibrated astrophysical coverage.')
    (BASE/'validation.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='checks'},indent=2))


if __name__=='__main__':main()
