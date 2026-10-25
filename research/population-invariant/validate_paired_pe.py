"""Independent source-coordinate sums/covariances and joint full products."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.integrate import cumulative_simpson,cumulative_trapezoid
from scipy.special import beta,ndtr

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
CAT=ROOT/'research/population-catalog'


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    path=BASE/'paired-pe-results.json';report=json.loads(path.read_text())
    with np.load(BASE/'paired-pe-checkpoints.npz',allow_pickle=False) as cp:a={k:cp[k] for k in cp.files}
    with np.load(BASE/'paired-o3-checkpoints.npz',allow_pickle=False) as cp:
        b={k:cp[k] for k in cp.files if k in ['native_grid','native_fs','effective_grid','full','selections','full_quantiles','contrasts','analytic_contrast_covariance']}
    with np.load(BASE/'control-checkpoints.npz',allow_pickle=False) as cp:zg=cp['zgrid'];dv=cp['zdv'];znorm=float(cp['znorm'])
    interesting={0,127,255,383,511}
    for values in [a['contrasts_event'],a['contrasts_joint']]:
        for k in range(values.shape[1]):interesting.update([int(np.argmin(values[:,k])),int(np.argmax(values[:,k]))])
    details=[];reconstructed=[];bank_covariances=[]
    DH=299792.458/67.9
    for pos,info in enumerate(report['event_banks']):
        index=info['index'];source=CAT/f'event-{index}.npz';assert sha(source)==info['source_cache_sha256']
        with np.load(source,allow_pickle=False) as bank:d={k:bank[k] for k in bank.files}
        m=d['mass_1'];q=d['mass_ratio'];z=d['redshift']
        a1=np.sqrt(sum(d[f'spin_1{axis}']**2 for axis in 'xyz'));a2=np.sqrt(sum(d[f'spin_2{axis}']**2 for axis in 'xyz'))
        c1=d['spin_1z']/a1;c2=d['spin_2z']/a2;x=(d['spin_1z']+q*d['spin_2z'])/(1+q)
        end=max(2.3,float(z.max())+.01);grid=np.linspace(0,end,int(np.ceil(end/(2.3/65536)))+1)
        integral=cumulative_simpson(1/np.sqrt(.3065*(1+grid)**3+.6935),x=grid,initial=0)
        dc=DH*np.interp(z,grid,integral);dl=(1+z)*dc;ddz=dc+(1+z)*DH/np.sqrt(.3065*(1+z)**3+.6935)
        prior=dl**2*ddz*m*(1+z)**2/4
        pm=.9*m**-2/(1/3-1/100)+.1*np.exp(-.5*((m-35)/5)**2)/(5*np.sqrt(2*np.pi)*(ndtr(13)-ndtr(-6.4)))
        pm*=((m>=3)&(m<=100));pq=np.where((m>3)&(q>=3/m)&(q<=1),2*q/np.maximum(1-(3/m)**2,1e-12),0)
        pz=np.interp(z,zg,dv)*(1+z)/znorm*((z>=0)&(z<=2.3));base=pm*pq*pz/prior
        native=[];effective=[];moment=(a1*a1+q*q*a2*a2)/(3*(1+q)**2)
        negative=np.where(x<0,2.,np.where(x>0,0.,1.));positive=2-negative
        for aa,bb in [(1.5,8.5),(2.,5.)]:
            radial=(a1*a2)**(aa-1)*((1-a1)*(1-a2))**(bb-1)/beta(aa,bb)**2
            for sigma in [.3,.6,1.2]:
                norm=sigma*np.sqrt(2*np.pi)*(ndtr(0)-ndtr(-2/sigma))
                angle=np.exp(-.5*((c1-1)**2+(c2-1)**2)/sigma**2)/norm**2
                native.append(np.column_stack([base*radial*angle,base*radial/4]))
            for lam in [0.,4.,16.]:
                angle=(1+lam*x*x)/(1+lam*moment)/4
                effective.append(np.column_stack([base*radial*angle*negative,base*radial*angle*positive]))
        w=np.stack([np.stack(native,axis=1),np.stack(effective,axis=1)],axis=1).reshape(len(z),24)
        full=w.mean(axis=0);wanted=a['full_means'][pos].reshape(24)
        assert np.array_equal(full==0,wanted==0)
        nz=wanted!=0;mean_error=float(np.max(abs(full[nz]/wanted[nz]-1)));assert mean_error<2e-7
        n=report['retained_rows_per_event'];size=len(z)
        covariance=(1-n/size)/n*(w.T@w-size*np.outer(full,full))/(size-1)
        oldcov=a['event_mean_covariances'][pos]
        covariance_error=float(np.linalg.norm(covariance-oldcov)/np.linalg.norm(oldcov));assert covariance_error<2e-7
        rng=np.random.default_rng(np.random.SeedSequence([report['event_seed'],index]));fingerprint=hashlib.sha256();sum_errors=[]
        for rep in range(512):
            take=rng.choice(size,n,replace=False);fingerprint.update(take.astype('<i8').tobytes())
            if rep in interesting:
                sampled=w[take].mean(axis=0);expected=a['event_means'][rep,pos].reshape(24)
                assert np.array_equal(sampled==0,expected==0)
                nz=expected!=0;error=float(np.max(abs(sampled[nz]/expected[nz]-1)));assert error<2e-7;sum_errors.append(error)
        assert fingerprint.hexdigest()==info['all512_index_stream_sha256']
        reconstructed.append(full.reshape(2,6,2));bank_covariances.append(covariance)
        details.append(dict(index=index,rows=size,full_mean_relative_error=mean_error,covariance_relative_error=covariance_error,
            maximum_sampled_mean_relative_error=max(sum_errors),all512_index_fingerprints_checked=True,
            zero_weight_coordinates=np.flatnonzero(full==0).tolist()))
        if (pos+1)%10==0:print(f'Independent event/covariance checks: {pos+1}/69',flush=True)
    # Every event-only and joint posterior is checked with ordinary products,
    # without calling the producer's logarithmic posterior routine.
    grids=[b['native_grid'],b['effective_grid']];fractions=[b['native_fs'],np.broadcast_to(b['effective_grid'],(6,len(b['effective_grid'])))]
    max_event=0.;max_joint=0.;event_quantiles=[];joint_quantiles=[]
    for rep in range(512):
        eq=np.empty((2,2,6,3));jq=np.empty_like(eq)
        for t in range(2):
            grid=grids[t];fs=fractions[t];mu=a['event_means'][rep,:,t]
            factors=fs[None,:,:]*mu[:,:,0,None]+(1-fs[None,:,:])*mu[:,:,1,None]
            middle=len(grid)//2;numerator=np.prod((factors/factors[:,:,middle,None]).astype(np.longdouble),axis=0)
            for p in range(2):
                for selection,target in [(b['full'][p,t],eq),(b['selections'][rep,p,t],jq)]:
                    efficiency=fs*selection[:,0,None]+(1-fs)*selection[:,1,None]
                    density=numerator*(efficiency[:,middle,None]/efficiency).astype(np.longdouble)**69
                    density=np.asarray(density/density.max(axis=1)[:,None],dtype=float)
                    cdf=cumulative_trapezoid(density,grid,axis=1,initial=0);cdf/=cdf[:,-1,None]
                    target[p,t]=[np.interp([.05,.5,.95],c,grid) for c in cdf]
        max_event=max(max_event,float(np.max(abs(eq-a['quantiles_event'][rep]))));max_joint=max(max_joint,float(np.max(abs(jq-a['quantiles_joint'][rep]))))
        event_quantiles.append(eq);joint_quantiles.append(jq)
    assert max(max_event,max_joint)<2e-11
    event_values=np.array(event_quantiles)[...,1].reshape(512,24)@b['contrasts'].T
    joint_values=np.array(joint_quantiles)[...,1].reshape(512,24)@b['contrasts'].T
    error=float(max(np.max(abs(event_values-a['contrasts_event'])),np.max(abs(joint_values-a['contrasts_joint']))));assert error<3e-9
    propagated=np.zeros((24,24))
    for covariance,J in zip(bank_covariances,a['median_event_derivative']):propagated+=J@covariance@J.T
    predicted=b['contrasts']@propagated@b['contrasts'].T
    covariance_error=float(np.linalg.norm(predicted-a['analytic_event_contrast_covariance'])/np.linalg.norm(predicted));assert covariance_error<2e-7
    # Validate representative nondegenerate sampling-direction derivatives by
    # full-product finite differences on the unthinned banks.
    mu=a['full_means'];gradient_errors=[];original_errors=[];complex_errors=[];complex_step_changes=[]
    selected=[i for i,r in enumerate(report['event_banks']) if r['index'] in [0,19,50,72]]
    for t in range(2):
        grid=grids[t];fs=fractions[t]
        factors=fs[None,:,:]*mu[:,t,:,0,None]+(1-fs[None,:,:])*mu[:,t,:,1,None]
        middle=len(grid)//2;num=np.prod((factors/factors[:,:,middle,None]).astype(np.longdouble),axis=0)
        for p in range(2):
            eff=fs*b['full'][p,t,:,0,None]+(1-fs)*b['full'][p,t,:,1,None]
            base_density=num*(eff[:,middle,None]/eff).astype(np.longdouble)**69
            for i in selected:
                for k in range(6):
                    if np.any(mu[i,t,k]==0):continue
                    base=np.asarray(base_density[k],dtype=float)
                    cdf=cumulative_trapezoid(base,grid,initial=0);cdf/=cdf[-1];location=np.interp(.5,cdf,grid)
                    j=int(np.clip(np.searchsorted(grid,location)-1,0,len(grid)-2))
                    minor=int(np.argmin(mu[i,t,k]));minor_basis=fs[k] if minor==0 else 1-fs[k]
                    alpha=minor_basis*mu[i,t,k,minor]/factors[i,k];complex_values=[]
                    for eps in [1e-12,1e-20]:
                        density=base*(1+alpha*np.expm1(1j*eps))
                        cc=cumulative_trapezoid(density,grid,initial=0);cc/=cc[-1]
                        qc=grid[j]+(.5-cc[j])/(cc[j+1]-cc[j])*(grid[j+1]-grid[j])
                        complex_values.append(float(qc.imag/eps))
                    complex_step_changes.append(abs(complex_values[0]/complex_values[1]-1))
                    for c,basis in enumerate([fs[k],1-fs[k]]):
                        out=np.ravel_multi_index((p,t,k),(2,2,6));inp=np.ravel_multi_index((t,k,c),(2,6,2))
                        analytic=a['median_event_derivative'][i,out,inp]
                        expected_log=complex_values[-1] if c==minor else -complex_values[-1]
                        complex_errors.append(abs(expected_log/(analytic*mu[i,t,k,c])-1))
                        for eps,errors in [(1e-6,original_errors),(1e-5,gradient_errors)]:
                            step=mu[i,t,k,c]*eps;quantiles=[]
                            for sign in [1,-1]:
                                ratio=1+sign*step*basis/factors[i,k]
                                density=base*ratio;cc=cumulative_trapezoid(density,grid,initial=0);cc/=cc[-1]
                                quantiles.append(float(np.interp(.5,cc,grid)))
                            finite=(quantiles[0]-quantiles[1])/(2*step);errors.append(abs(finite/analytic-1))
    assert max(gradient_errors)<1e-4,max(gradient_errors)
    assert max(complex_errors)<1e-8 and max(complex_step_changes)<1e-8
    rows=[]
    for k,row in enumerate(report['contrasts']):
        lin=a['linear_event_contrast_perturbations'][:,k]
        rows.append(dict(name=row['name'],event_linear_empirical_over_analytic_variance=float(np.var(lin,ddof=1)/predicted[k,k]),
            event_nonlinear_over_linear_empirical_variance=float(np.var(event_values[:,k],ddof=1)/np.var(lin,ddof=1))))
    result=dict(schema='blackhole-joint-finite-bank-independent-validation/1',events=69,event_bank_checks=details,
        event_index_streams_checked=69,replicates_per_stream=512,sample_mean_replicates_checked=sorted(interesting),
        maximum_full_mean_relative_error=max(r['full_mean_relative_error'] for r in details),maximum_bank_covariance_relative_error=max(r['covariance_relative_error'] for r in details),
        maximum_sampled_mean_relative_error=max(r['maximum_sampled_mean_relative_error'] for r in details),
        all_event_only_full_product_replicates=512,all_joint_full_product_replicates=512,
        maximum_event_only_quantile_error=max_event,maximum_joint_quantile_error=max_joint,maximum_contrast_error_pp=error,
        linear_event_contrast_covariance_relative_error=covariance_error,maximum_selected_derivative_relative_error=float(max(gradient_errors)),
        original_relative_derivative_gate=1e-4,original_step=1e-6,original_maximum_derivative_relative_error=float(max(original_errors)),
        original_step_derivative_gate_passed=bool(max(original_errors)<1e-4),revised_finite_difference_step=1e-5,
        maximum_complex_log_odds_derivative_relative_error=float(max(complex_errors)),maximum_complex_step_pair_relative_change=float(max(complex_step_changes)),
        derivative_diagnosis='Original1e-6 finite difference loses precision on a weak log-mean derivative. Failure retained. Unchanged1e-4 relative gate passes at1e-5; independent fixed-bracket complex differentiation checks log odds at1e-12 and1e-20 without subtracting nearly equal real medians.',
        nonlinear_diagnostics=rows,result_sha256=sha(path),
        scope='Independent Cartesian source coordinates, explicit prior quadrature, all69 covariance matrices and512 draw fingerprints per source. Source-mean spotchecks include contrast extremes; all512 event/joint posterior products checked. Conditional finite-bank design only, not physical likelihood-tail or real-campaign calibration.')
    (BASE/'paired-pe-validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['events','maximum_full_mean_relative_error','maximum_bank_covariance_relative_error','maximum_sampled_mean_relative_error','maximum_event_only_quantile_error','maximum_joint_quantile_error','maximum_selected_derivative_relative_error']},indent=2))


if __name__=='__main__':main()
