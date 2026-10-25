"""Independent compact-bank covariance and all512 full-product propagation."""
from pathlib import Path
import hashlib
import json
import platform
import sys
import numpy as np
from scipy.integrate import cumulative_trapezoid

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
SEL=ROOT/'research/population-selection-audit'
sys.path.insert(0,str(SEL))
from run_selection_study import init_model,physical,weights


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    init_model();path=BASE/'paired-o3-results.json';report=json.loads(path.read_text())
    source=json.loads((SEL/'o3-provenance.json').read_text());assert sha(SEL/'o3-bank.npz')==source['cache_sha256']
    with np.load(BASE/'paired-o3-checkpoints.npz',allow_pickle=False) as cp:a={k:cp[k] for k in cp.files}
    with np.load(BASE/'control-checkpoints.npz',allow_pickle=False) as cp:effective_union=cp['o3_union_weights']
    with np.load(SEL/'o3-bank.npz',allow_pickle=False) as cp:raw={k[5:]:cp[k] for k in cp.files if k.startswith('full_')}
    union=np.any([raw[k]>1 for k in raw if k.startswith('ifar_')],axis=0)
    matched=np.any([raw[k]>1 for k in ['ifar_gstlal','ifar_mbta','ifar_pycbc_bbh','ifar_pycbc_hyperbank']],axis=0)
    count=len(union);retained=report['retained_rows'];assert count==report['accepted_o3_rows']
    native_union=weights(physical({k:v[union] for k,v in raw.items()}))
    assert len(effective_union)==union.sum()
    compact=np.stack([native_union,effective_union],axis=1)
    compact=np.stack([compact,compact*matched[union,None,None,None]],axis=1).reshape(int(union.sum()),48)
    sums=compact.sum(axis=0)
    # Independent centered cross-moment formula, explicitly including all
    # failed-detection zero rows in the finite-population denominator.
    population_cov=(compact.T@compact-np.outer(sums,sums)/count)/(count-1)
    factor=source['attributes']['analysis_time_s']/source['combined_attributes']['analysis_time_s']/source['attributes']['total_generated']
    analytic=(factor*count)**2*(1-retained/count)/retained*population_cov
    expected=a['analytic_selection_covariance']
    cov_error=float(np.linalg.norm(analytic-expected)/np.linalg.norm(expected));assert cov_error<1e-10
    flat_xi=a['selections'].reshape(512,48);centered=flat_xi-flat_xi.mean(axis=0)
    empirical=centered.T@centered/511
    empirical_error=float(np.linalg.norm(empirical-a['empirical_selection_covariance'])/np.linalg.norm(empirical));assert empirical_error<1e-12
    # Regenerate every draw fingerprint. Independently re-sum selected draws
    # spanning regular checkpoints and extrema of every reported contrast.
    interesting={0,63,127,255,383,511}
    for k in range(a['contrast_values'].shape[1]):
        interesting.add(int(np.argmin(a['contrast_values'][:,k])));interesting.add(int(np.argmax(a['contrast_values'][:,k])))
    full_values=np.zeros((count,48));full_values[union]=compact
    fixed=a['full'].reshape(48)-factor*sums
    rng=np.random.default_rng(report['seed']);sum_errors=[]
    for i in range(512):
        take=rng.choice(count,retained,replace=False)
        assert hashlib.sha256(take.astype('<i8').tobytes()).hexdigest()==report['row_index_sha256'][i]
        if i in interesting:
            xi=fixed+factor*count/retained*full_values[take].sum(axis=0)
            error=float(np.max(abs(xi/flat_xi[i]-1)));assert error<1e-10;sum_errors.append(error)
    # Full69-event products are precomputed independently of the selection.
    # Every simulated policy/metric/model posterior is then verified, not only
    # a subset of displayed quantiles.
    grids=[a['native_grid'],a['effective_grid']]
    fractions=[a['native_fs'],np.broadcast_to(grids[1],(6,len(grids[1])))]
    means=[a['native_means'],a['effective_means']];products=[]
    for grid,fs,mu in zip(grids,fractions,means):
        terms=fs[None,:,:]*mu[:,:,0,None]+(1-fs[None,:,:])*mu[:,:,1,None]
        middle=len(grid)//2
        product=np.prod((terms/terms[:,:,middle,None]).astype(np.longdouble),axis=0)
        products.append(product)
    def product_quantiles(xis):
        out=np.empty((len(xis),2,2,6,3))
        for i,xi in enumerate(xis):
            for p in range(2):
                for t in range(2):
                    grid=grids[t];fs=fractions[t];eff=fs*xi[p,t,:,0,None]+(1-fs)*xi[p,t,:,1,None]
                    density=products[t]*(eff[:,len(grid)//2,None]/eff).astype(np.longdouble)**69
                    density=np.asarray(density/density.max(axis=1)[:,None],dtype=float)
                    cdf=cumulative_trapezoid(density,grid,axis=1,initial=0);cdf/=cdf[:,-1,None]
                    out[i,p,t]=[np.interp([.05,.5,.95],v,grid) for v in cdf]
        return out
    verified=product_quantiles(a['selections']);qerr=float(np.max(abs(verified-a['quantiles'])));assert qerr<2e-11
    full=product_quantiles(a['full'][None])[0];ferr=float(np.max(abs(full-a['full_quantiles'])));assert ferr<2e-11
    contrasts=verified[...,1].reshape(512,24)@a['contrasts'].T
    cerr=float(np.max(abs(contrasts-a['contrast_values'])));assert cerr<3e-9
    # Separate finite512-draw covariance fluctuations from nonlinear posterior
    # propagation; these are distinct reasons to depart from a delta method.
    G=a['contrasts']@a['median_derivative']
    delta=(flat_xi-a['full'].reshape(48))@G.T
    linear_empirical=np.cov(delta,rowvar=False,ddof=1)
    nonlinear_cov=np.cov(contrasts,rowvar=False,ddof=1)
    analytic_contrast=G@analytic@G.T
    contrast_cov_error=float(np.linalg.norm(analytic_contrast-a['analytic_contrast_covariance'])/np.linalg.norm(analytic_contrast));assert contrast_cov_error<1e-10
    nonlinear=[]
    nominal=a['full_quantiles'][...,1].reshape(24)@a['contrasts'].T
    for k,row in enumerate(report['contrasts']):
        nonlinear.append(dict(name=row['name'],
            linear_empirical_over_analytic_variance=float(linear_empirical[k,k]/analytic_contrast[k,k]),
            nonlinear_over_linear_empirical_variance=float(nonlinear_cov[k,k]/linear_empirical[k,k]),
            nonlinear_minus_linear_5_50_95_pp=np.quantile(contrasts[:,k]-nominal[k]-delta[:,k],[.05,.5,.95]).tolist()))
    result=dict(schema='blackhole-full69-paired-o3-independent-validation/1',all_draw_hashes_checked=512,
        selection_sum_replicates_checked=sorted(interesting),maximum_spotcheck_selection_relative_error=max(sum_errors),
        compact_bank_analytic_covariance_relative_error=cov_error,empirical_covariance_relative_error=empirical_error,
        all_full_product_replicates_checked=512,full_product_quantile_error=qerr,full_bank_quantile_error=ferr,
        contrast_value_error_pp=cerr,contrast_covariance_relative_error=contrast_cov_error,nonlinear_diagnostics=nonlinear,
        runtime=dict(python=platform.python_version(),numpy=np.__version__),result_sha256=sha(path),
        scope='All512 draw fingerprints and posterior products; selection-sum spotchecks include every contrast extremum. Compact covariance includes the failed-search zero rows. Shares previously validated model weights and fixed event moments. No real-campaign or PE-chain coverage claim.')
    (BASE/'paired-o3-validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['all_draw_hashes_checked','selection_sum_replicates_checked','compact_bank_analytic_covariance_relative_error','full_product_quantile_error','contrast_value_error_pp']},indent=2))


if __name__=='__main__':main()
