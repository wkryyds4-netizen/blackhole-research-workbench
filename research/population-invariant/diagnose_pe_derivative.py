"""Diagnose the retained relative finite-difference derivative failure."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import cumulative_trapezoid

BASE=Path(__file__).resolve().parent


def main():
    with np.load(BASE/'paired-pe-checkpoints.npz',allow_pickle=False) as cp:mu=cp['full_means'];J=cp['median_event_derivative']
    with np.load(BASE/'paired-o3-checkpoints.npz',allow_pickle=False) as cp:
        grids=[cp['native_grid'],cp['effective_grid']];fractions=[cp['native_fs'],np.broadcast_to(grids[1],(6,len(grids[1])))];selection=cp['full']
    report=json.loads((BASE/'paired-pe-results.json').read_text());selected=[i for i,r in enumerate(report['event_banks']) if r['index'] in [0,19,50,72]]
    rows=[]
    for t in range(2):
        grid=grids[t];fs=fractions[t];factors=fs[None,:,:]*mu[:,t,:,0,None]+(1-fs[None,:,:])*mu[:,t,:,1,None]
        middle=len(grid)//2;num=np.prod((factors/factors[:,:,middle,None]).astype(np.longdouble),axis=0)
        for p in range(2):
            eff=fs*selection[p,t,:,0,None]+(1-fs)*selection[p,t,:,1,None]
            densities=num*(eff[:,middle,None]/eff).astype(np.longdouble)**69
            for i in selected:
                for k in range(6):
                    if np.any(mu[i,t,k]==0):continue
                    base=np.asarray(densities[k],dtype=float);cdf=cumulative_trapezoid(base,grid,initial=0);cdf/=cdf[-1]
                    q=np.interp(.5,cdf,grid);j=int(np.clip(np.searchsorted(grid,q)-1,0,len(grid)-2))
                    minor=int(np.argmin(mu[i,t,k]));basis=fs[k] if minor==0 else 1-fs[k]
                    alpha=basis*mu[i,t,k,minor]/factors[i,k]
                    eps=1e-20;density=base*(1+alpha*np.expm1(1j*eps));cc=cumulative_trapezoid(density,grid,initial=0);cc/=cc[-1]
                    qc=grid[j]+(.5-cc[j])/(cc[j+1]-cc[j])*(grid[j+1]-grid[j]);complex_log=float(qc.imag/eps)
                    for c,basis in enumerate([fs[k],1-fs[k]]):
                        estimates=[]
                        for eps in [1e-4,1e-5,1e-6,1e-7]:
                            values=[]
                            for sign in [1,-1]:
                                density=base*(1+sign*eps*basis*mu[i,t,k,c]/factors[i,k])
                                cc=cumulative_trapezoid(density,grid,initial=0);cc/=cc[-1];values.append(float(np.interp(.5,cc,grid)))
                            estimates.append(dict(step=eps,log_coordinate_derivative=(values[0]-values[1])/(2*eps),
                                upper_median_bin=int(np.searchsorted(grid,values[0])-1),lower_median_bin=int(np.searchsorted(grid,values[1])-1)))
                        out=np.ravel_multi_index((p,t,k),(2,2,6));inp=np.ravel_multi_index((t,k,c),(2,6,2))
                        analytic=float(J[i,out,inp]*mu[i,t,k,c]);finite=next(v['log_coordinate_derivative'] for v in estimates if v['step']==1e-6)
                        expected=complex_log if c==minor else -complex_log
                        rows.append(dict(source_index=report['event_banks'][i]['index'],policy=p,metric=t,model=k,component=c,
                            full_mean=float(mu[i,t,k,c]),baseline_median_bin=j,analytic_log_derivative=analytic,complex_minor_log_derivative=expected,
                            relative_error_at_original_step=abs(finite/analytic-1),finite_differences=estimates))
    worst=max(rows,key=lambda r:r['relative_error_at_original_step'])
    result=dict(schema='blackhole-pe-derivative-failure-diagnosis/1',original_relative_gate=1e-4,original_gate_passed=False,
        original_max_relative_error=worst['relative_error_at_original_step'],worst=worst,rows=rows,
        limits='Derivative diagnosis only; complex differentiation holds the baseline median interpolation bracket fixed. Does not change draw means, posterior quantiles, or scientific sampling design.')
    (BASE/'pe-derivative-diagnosis.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(worst,indent=2))


if __name__=='__main__':main()
