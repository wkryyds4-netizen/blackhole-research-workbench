"""Shared finite-bank thinning for the full69 native and invariant controls."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
from scipy.integrate import cumulative_trapezoid

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
SEL=ROOT/'research/population-selection-audit'
CAT=ROOT/'research/population-catalog'
sys.path.insert(0,str(SEL))
from run_selection_study import physical,init_model,weights as native_weights,coordinates
from thinning_benchmark import srs_covariance
from run_control import weights as effective_weights


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def quantile_and_gradient(grid,fs,numerator,xi,n_events=69):
    """Quantiles and exact piecewise-grid median derivative wrt selection."""
    efficiency=fs*xi[:,0,None]+(1-fs)*xi[:,1,None]
    logs=numerator-n_events*np.log(efficiency)
    p=np.exp(logs-logs.max(axis=1)[:,None]);p/=np.trapezoid(p,grid,axis=1)[:,None]
    cdf=cumulative_trapezoid(p,grid,axis=1,initial=0);cdf/=cdf[:,-1,None]
    q=np.array([np.interp([.05,.5,.95],c,grid) for c in cdf]);grad=np.empty((6,2))
    for k in range(6):
        location=q[k,1];j=int(np.clip(np.searchsorted(grid,location)-1,0,len(grid)-2))
        slope=(cdf[k,j+1]-cdf[k,j])/(grid[j+1]-grid[j])
        for c,coefficient in enumerate([fs[k],1-fs[k]]):
            derivative=-n_events*coefficient/efficiency[k]
            integ=cumulative_trapezoid(p[k]*derivative,grid,initial=0)
            grad[k,c]=-(np.interp(location,grid,integ)-.5*integ[-1])/slope
    assert np.max(abs(np.sum(grad*xi,axis=1)))<1e-10,'Selection-scale invariance failed'
    return q,grad


def contrast_matrix():
    rows=[];names=[]
    def index(policy,metric,model):return np.ravel_multi_index((policy,metric,model),(2,2,6))
    def add(name,terms):
        row=np.zeros(24)
        for p,t,m,c in terms:row[index(p,t,m)]+=100*c
        rows.append(row);names.append(name)
    for p,policy in enumerate(['union','matched']):
        for a,label in enumerate(['low_spin','baseline_spin']):add(f'{policy}:native_width:{label}',[(p,0,3*a,1),(p,0,3*a+2,-1)])
        for l in range(3):add(f'{policy}:effective_amplitude:lambda{[0,4,16][l]}',[(p,1,3+l,1),(p,1,l,-1)])
        for a,label in enumerate(['low_spin','baseline_spin']):add(f'{policy}:effective_shape:{label}',[(p,1,3*a+2,1),(p,1,3*a,-1)])
    for m in range(6):add(f'effective_cwb_minus_matched:model{m}',[(0,1,m,1),(1,1,m,-1)])
    for a,label in enumerate(['low_spin','baseline_spin']):add(f'native_matched_minus_union_width:{label}',[(1,0,3*a,1),(1,0,3*a+2,-1),(0,0,3*a,-1),(0,0,3*a+2,1)])
    return names,np.array(rows)


def main():
    init_model();nr=json.loads((CAT/'native-catalog-results.json').read_text());nv=json.loads((CAT/'native-catalog-validation.json').read_text())
    er=json.loads((BASE/'control-results.json').read_text());ev=json.loads((BASE/'control-validation.json').read_text())
    assert nr['completed_full_source_events']==69 and nv['result_sha256']==sha(CAT/'native-catalog-results.json')
    assert ev['result_sha256']==sha(BASE/'control-results.json')
    source=json.loads((SEL/'o3-provenance.json').read_text());study=json.loads((SEL/'selection-results.json').read_text())
    assert sha(SEL/'o3-bank.npz')==source['cache_sha256']
    with np.load(SEL/'o3-bank.npz',allow_pickle=False) as bank:
        raw={k[5:]:bank[k] for k in bank.files if k.startswith('full_')};actual=bank['combined_o3_to_full']
    with np.load(CAT/'native-catalog-checkpoints.npz',allow_pickle=False) as cp:nmeans=cp['means'];r=cp['r'];nfs=cp['fs']
    with np.load(BASE/'control-checkpoints.npz',allow_pickle=False) as cp:emeans=cp['means'];g=cp['g'];early_effective=cp['early_weights'].sum(axis=0)
    with np.load(SEL/'paired-policy-checkpoints.npz',allow_pickle=False) as cp:old_native_xi=cp['selections']
    data=physical(raw);weights=np.stack([native_weights(data),effective_weights(data)],axis=1)
    matched=np.maximum.reduce([raw[k] for k in ['ifar_gstlal','ifar_mbta','ifar_pycbc_bbh','ifar_pycbc_hyperbank']])>1
    union=matched|(raw['ifar_cwb']>1)
    bank=np.stack([weights*union[:,None,None,None],weights*matched[:,None,None,None]],axis=1)
    m=len(bank);n=len(actual);assert (m,n)==(284672,128454)
    factor=source['attributes']['analysis_time_s']/source['combined_attributes']['analysis_time_s']/source['attributes']['total_generated']
    early_native=np.array(study['selection_components']['o1'])+np.array(study['selection_components']['o2'])
    fixed=np.stack([early_native,early_effective]);full=fixed[None,:,:,:]+factor*bank.sum(axis=0)
    historical=fixed[None,:,:,:]+factor*m/n*bank[actual].sum(axis=0)
    grids=[r,g];fractions=[nfs,np.broadcast_to(g,(6,len(g)))];means=[nmeans,emeans]
    numerators=[]
    for fs,mu in zip(fractions,means):
        with np.errstate(divide='ignore'):numerators.append(np.log(fs[None,:,:]*mu[:,:,0,None]+(1-fs[None,:,:])*mu[:,:,1,None]).sum(axis=0))
    def evaluate(xi,gradient=False):
        qs=[];derivatives=[]
        for p in range(2):
            qpolicy=[];gpolicy=[]
            for t in range(2):
                q,j=quantile_and_gradient(grids[t],fractions[t],numerators[t],xi[p,t]);qpolicy.append(q);gpolicy.append(j)
            qs.append(qpolicy);derivatives.append(gpolicy)
        return (np.array(qs),np.array(derivatives)) if gradient else np.array(qs)
    full_q,jac=evaluate(full,True);historical_q=evaluate(historical)
    assert np.max(abs(full_q[1,0]-nr['quantiles']))<1e-11
    for p,name in enumerate(['with_cwb','matched_filter']):assert np.max(abs(full_q[p,1]-er['scenarios'][name]['quantiles']))<1e-11
    # Independent local finite differences validate the posterior Jacobian.
    gradient_errors=[]
    for p in range(2):
        for t in range(2):
            for k in range(6):
                for c in range(2):
                    step=full[p,t,k,c]*1e-6;up=full[p,t].copy();down=up.copy();up[k,c]+=step;down[k,c]-=step
                    qu,_=quantile_and_gradient(grids[t],fractions[t],numerators[t],up)
                    qd,_=quantile_and_gradient(grids[t],fractions[t],numerators[t],down)
                    numerical=(qu[k,1]-qd[k,1])/(2*step)
                    gradient_errors.append(abs(numerical/jac[p,t,k,c]-1))
    assert max(gradient_errors)<1e-4,max(gradient_errors)
    covariance=(factor*m)**2*srs_covariance(bank.reshape(m,-1),n)
    derivative=np.zeros((24,48))
    for p in range(2):
        for t in range(2):
            for k in range(6):
                out=np.ravel_multi_index((p,t,k),(2,2,6))
                for c in range(2):derivative[out,np.ravel_multi_index((p,t,k,c),(2,2,6,2))]=jac[p,t,k,c]
    rng=np.random.default_rng(20261025);xis=[];qs=[];index_hashes=[]
    for i in range(512):
        take=rng.choice(m,n,replace=False)
        xi=fixed[None,:,:,:]+factor*m/n*bank[take].sum(axis=0)
        assert np.max(abs(xi[:,0]/old_native_xi[i]-1))<1e-11
        xis.append(xi);qs.append(evaluate(xi));index_hashes.append(hashlib.sha256(take.astype('<i8').tobytes()).hexdigest())
        if (i+1)%64==0:print(f'Full69 paired O3 thinning: {i+1}/512',flush=True)
    xis=np.array(xis);qs=np.array(qs)
    empirical=np.cov(xis.reshape(512,-1),rowvar=False,ddof=1)
    names,C=contrast_matrix();full_values=C@full_q[...,1].reshape(-1)
    values=qs[...,1].reshape(512,-1)@C.T
    predicted=C@derivative@covariance@derivative.T@C.T
    observed=np.cov(values,rowvar=False,ddof=1)
    contrasts=[]
    for k,name in enumerate(names):
        contrasts.append(dict(name=name,full_bank_value_pp=float(full_values[k]),conditional_thinning_5_50_95_pp=np.quantile(values[:,k],[.05,.5,.95]).tolist(),
            empirical_std_pp=float(np.std(values[:,k],ddof=1)),linearized_std_pp=float(np.sqrt(predicted[k,k])),
            nonpositive_fraction=float(np.mean(values[:,k]<=0))))
    report=dict(schema='blackhole-full69-paired-o3-control/1',events=69,replicates=512,seed=20261025,accepted_o3_rows=m,retained_rows=n,
        targets=['native_primary_negative_tilt','negative_effective_spin'],policies=['union_including_cwb','matched_filter_only'],
        full_quantiles=full_q.tolist(),historical_row_control_quantiles=historical_q.tolist(),contrasts=contrasts,
        empirical_selection_variance=np.diag(empirical).tolist(),
        selection_variance_ratios=(np.diag(empirical)/np.diag(covariance)).tolist(),
        contrast_variance_ratios=(np.diag(observed)/np.diag(predicted)).tolist(),
        median_shift_5_50_95_pp=np.moveaxis(np.quantile(100*(qs[...,1]-full_q[None,...,1]),[.05,.5,.95],axis=0),0,-1).tolist(),
        maximum_jacobian_finite_difference_relative_error=float(max(gradient_errors)),row_index_sha256=index_hashes,
        input_hashes=dict(native=nv['result_sha256'],effective=ev['result_sha256'],o3=source['cache_sha256'],design=sha(BASE/'paired-o3-design.json')),
        scope='Conditional SRSWOR artificial design on the fixed accepted O3 bank. Both metrics, policies and all models share each draw. Events and O1/O2 fixed. Previously used seed is reused, not additional independent trials.',
        limits='No actual historical selection law, PE-chain correlation, missing likelihood, astrophysical interval, p-value or calibrated real-campaign uncertainty is established. r and g are different functionals; no cross-target subtraction is performed.')
    (BASE/'paired-o3-results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    np.savez_compressed(BASE/'paired-o3-checkpoints.npz',selections=xis,quantiles=qs,full=full,historical=historical,full_quantiles=full_q,
        analytic_selection_covariance=covariance,empirical_selection_covariance=empirical,median_derivative=derivative,contrasts=C,contrast_values=values,
        analytic_contrast_covariance=predicted,empirical_contrast_covariance=observed,native_grid=r,native_fs=nfs,effective_grid=g,
        native_numerator=numerators[0],effective_numerator=numerators[1],native_means=nmeans,effective_means=emeans)
    print(json.dumps(dict(contrasts=contrasts,maximum_jacobian_error=max(gradient_errors),selection_variance_ratio_range=[min(report['selection_variance_ratios']),max(report['selection_variance_ratios'])]),indent=2))


if __name__=='__main__':main()
