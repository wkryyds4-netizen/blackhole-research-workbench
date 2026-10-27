"""Independent prior/Cartesian/flow checks at fixed points from all69 sources.

This validates row weights, not the full integrals or astrophysical assumptions.
The scalar convolution, explicit standard-normal density and distance integral
are distinct from the production vector formulas.
"""
from pathlib import Path
import hashlib
import json
import sys
import h5py
import numpy as np
from scipy.integrate import quad
import torch
from conditional_channel_model import ConditionalChannels,FAMILIES,CHANNELS,A,UNIFORM_VOLUME
from projected_measure import effective_spin_prior,hidden_spin_density

BASE=Path(__file__).resolve().parent;ROOT=BASE.parents[1]
CAT=ROOT/'research/population-catalog';SEL=ROOT/'research/population-selection-audit'
sys.path.insert(0,str(SEL))
from run_selection_study import density_pair,physical,PARAMS

def sha(p):
    with p.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def independent_log_density4(f,theta):
    result=[]
    for name,flow in zip(CHANNELS,f.models):
        y=theta.copy();logjac=np.zeros(len(y))
        for col in [0,1,3]:
            bound=flow.upper[col];scale=flow.scale[col];v=theta[:,col]
            y[:,col]=np.log(v/(bound-v))/scale
            logjac-=np.log(abs(scale))+np.log(v)+np.log1p(-v/bound)
        y[:,2]=np.arctanh(theta[:,2]);logjac-=np.log1p(-theta[:,2]**2)
        raw=[.1,np.log(2)] if name=='CE' else [.1]
        with torch.no_grad():z,det=flow.network.forward(torch.as_tensor(y),conditional=flow.context(raw,len(y)))
        result.append(-.5*np.sum(z.numpy()**2,axis=1)-2*np.log(2*np.pi)+det.numpy()+logjac)
    return np.array(result).T


def source_prior(d):
    m=d['mass_1'];z=d['redshift'];H=299792.458/67.9
    ez=lambda v:np.sqrt(.3065*(1+v)**3+1-.3065)
    dc=np.array([H*quad(lambda t:1/ez(t),0,float(v),epsabs=1e-11,epsrel=1e-11)[0] for v in z])
    dl=dc*(1+z);derivative=dc+(1+z)*H/ez(z)
    return dl*dl*derivative*m*(1+z)**2/4


def audit(f,d,prior_override=None):
    m=d['mass_1'];q=d['mass_ratio'];z=d['redshift'];x=(d['a_1']*d['cos_tilt_1']+q*d['a_2']*d['cos_tilt_2'])/(1+q)
    valid=(m>=5)&(m<=100)&(m*q>=3)&(q>0)&(q<=1)&(z>0)&(z<=1.9)&(abs(x)<A)
    for i in [1,2]:valid&=(d[f'a_{i}']>0)&(d[f'a_{i}']<A)&(abs(d[f'cos_tilt_{i}'])<=1)
    prior=source_prior(d) if prior_override is None else prior_override
    _,production_prior,_=density_pair(d) if prior_override is None else (None,prior,None)
    prior_error=float(np.max(abs(production_prior/prior-1)))
    assert prior_error<1e-10,prior_error
    actual=f.weights(dict(d,prior=production_prior));expected=np.zeros_like(actual)
    # dMc/dm1 is obtained from Mc=m1*q^(3/5)/(1+q)^(1/5).
    derivative=q[valid]**.6/(1+q[valid])**.2
    theta=np.column_stack([m[valid]*derivative,q[valid],x[valid],z[valid]])
    if len(theta):
        densities=independent_log_density4(f,theta)
        for fi,family in enumerate(FAMILIES):
            projected=np.array([effective_spin_prior(v,r,A,family) for v,r in zip(x[valid],q[valid])])
            Q=np.array([hidden_spin_density(a,b,c,e,A,family) for a,b,c,e in zip(d['a_1'][valid],d['a_2'][valid],d['cos_tilt_1'][valid],d['cos_tilt_2'][valid])])
            expected[valid,fi,:5]=np.exp(densities+np.log(derivative[:,None])+np.log(Q[:,None])-np.log(projected[:,None])-np.log(prior[valid,None]))
            expected[valid,fi,5]=np.exp(np.log(Q)-np.log(projected)-np.log(prior[valid])-np.log(UNIFORM_VOLUME))
    assert np.array_equal(actual==0,expected==0),'Unexpected zero-weight disagreement'
    positive=expected>0;relative=float(np.max(abs(actual[positive]/expected[positive]-1))) if positive.any() else 0.
    assert relative<1e-5,relative
    return dict(rows=len(m),in_domain=int(valid.sum()),prior_relative_error=prior_error,weight_relative_error=relative)


def choose(d,seed):
    n=len(d['mass_1']);q=d['mass_ratio'];x=(d['a_1']*d['cos_tilt_1']+q*d['a_2']*d['cos_tilt_2'])/(1+q)
    extremes=[np.argmin(abs(x)),np.argmax(abs(x)),np.argmin(q),np.argmax(d['mass_1']),np.argmax(d['redshift'])]
    return np.unique(np.r_[np.random.default_rng(seed).choice(n,24,replace=False),extremes])


def main():
    torch.set_num_threads(1);f=ConditionalChannels();events=[]
    frozen=json.loads((ROOT/'results/model.json').read_text());provenance=json.loads((CAT/'event-provenance.json').read_text())
    for pos,e in enumerate(frozen['events']):
        ix=e['index'];p=CAT/f'event-{ix}.npz';info=next(r for r in provenance if r['index']==ix)
        assert sha(p)==info['cache_sha256']
        with np.load(p,allow_pickle=False) as bank:full={k:bank[k] for k in PARAMS}
        indices=choose(full,20261031+ix);d={k:v[indices] for k,v in full.items()}
        events.append(dict(index=ix,event=info['event'],selected_indices=indices.tolist(),source_cache_sha256=info['cache_sha256'],**audit(f,d)))
        if (pos+1)%10==0:print(f'Independent scalar/flow/prior point checks: {pos+1}/69 sources',flush=True)
    injections=[]
    with np.load(SEL/'o3-bank.npz',allow_pickle=False) as bank:
        raw={k[5:]:bank[k] for k in bank.files if k.startswith('full_')}
    sources=[('pure_bbh',raw)]
    rec=json.loads((ROOT/'research/population-mixture/download.json').read_text());p=ROOT/'research/population-mixture/sources'/rec['file']
    assert sha(p)==rec['sha256']
    keys=['mass1_source','mass2_source','redshift','sampling_pdf']+[f'spin{i}{ax}' for i in [1,2] for ax in 'xyz']
    with h5py.File(p,'r') as handle:sources.append(('official_mixture',{k:handle['injections'][k][:] for k in keys}))
    for name,raw in sources:
        physical_all=physical(raw);indices=choose(physical_all,20261101)
        d={k:v[indices] for k,v in physical_all.items()};selected={k:v[indices] for k,v in raw.items()}
        # Independent spherical-to-Cartesian angular integration and mass Jacobian.
        radius=[np.sqrt(sum(selected[f'spin{i}{ax}'].astype(float)**2 for ax in 'xyz')) for i in [1,2]]
        cart_to_seven=selected['mass1_source'].astype(float)*4*np.pi*np.pi*radius[0]**2*radius[1]**2
        prior=selected['sampling_pdf']*cart_to_seven
        error=float(np.max(abs(prior/d['prior']-1)));assert error<1e-12
        injections.append(dict(bank=name,selected_indices=indices.tolist(),cartesian_jacobian_relative_error=error,**audit(f,d,prior)))
    report=dict(schema='blackhole-channel-row-weight-independent-validation/1',status='pass',events=events,injections=injections,
                checked_events=len(events),checked_event_rows=sum(r['rows'] for r in events),
                scope='Fixed random/stress row checks on all69 sources and both O3 banks; not independent full-sum reintegration, posterior validation, or physical/campaign calibration.',
                input_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in [Path(__file__),BASE/'conditional_channel_model.py',BASE/'channel_flow.py',BASE/'projected_measure.py',BASE/'fast_measure.py',BASE/'likelihood-design.json',CAT/'event-provenance.json']})
    (BASE/'channel-weight-validation.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='pass',events=len(events),rows=report['checked_event_rows'],max_weight_relative_error=max(r['weight_relative_error'] for r in events+injections))),flush=True)


if __name__=='__main__':
    try:main()
    except Exception as error:
        p=BASE/'channel-weight-failures.json';history=json.loads(p.read_text()) if p.exists() else []
        history.append(dict(tool_sha256=sha(Path(__file__)),error=repr(error)))
        p.write_text(json.dumps(history,indent=2)+'\n',encoding='utf-8');raise
