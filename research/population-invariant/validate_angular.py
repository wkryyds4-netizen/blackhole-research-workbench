"""Independent split-domain quadrature of the proposed angular control."""
from pathlib import Path
import hashlib
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from invariant_model import angular_density

BASE=Path(__file__).resolve().parent


def integrate(q,a1,a2,lam,g,order=48):
    nodes,weights=leggauss(order)
    u=a1/(1+q);v=q*a2/(1+q)
    # Split at every intersection of x=0 with the c2 integration box. This
    # avoids counting quadrature leakage across the prescribed density jump.
    cuts=sorted(set([-1.,1.]+([float(t) for t in [-v/u,v/u] if -1<t<1] if u else [])))
    total=0.;negative=0.
    for left,right in zip(cuts[:-1],cuts[1:]):
        c1=(left+right)/2+(right-left)/2*nodes;w1=(right-left)/2*weights
        for c,weight in zip(c1,w1):
            boundary=np.clip(-u*c/v,-1,1) if v else (-1. if u*c>=0 else 1.)
            for lo,hi,is_negative in [(-1.,boundary,True),(boundary,1.,False)]:
                if hi<=lo:continue
                c2=(lo+hi)/2+(hi-lo)/2*nodes;w2=(hi-lo)/2*weights
                density=angular_density(q,a1,a2,c,np.array(c2),lam,g)
                integral=weight*float(w2@density);total+=integral
                if is_negative:negative+=integral
    return total,negative


def main():
    inputs=[(.05,.01,.99),(.2,.9,.02),(.65,.7,.6),(1.,.8,.8),(.999999,.99,.99),(.7,0.,.4),(.7,.4,0.),(.7,1e-9,.4)]
    rows=[]
    for q,a1,a2 in inputs:
        for lam in [0.,4.,16.]:
            for g in [0.,.1,.5,.9,1.]:
                norm,neg=integrate(q,a1,a2,lam,g)
                rows.append(dict(q=q,a1=a1,a2=a2,shape_lambda=lam,g=g,normalization=norm,negative_probability=neg))
    normerr=max(abs(r['normalization']-1) for r in rows);negerr=max(abs(r['negative_probability']-r['g']) for r in rows)
    assert normerr<1e-10 and negerr<1e-10,(normerr,negerr)
    # Both spins exactly zero have no physically meaningful negative fraction;
    # this zero-measure limit uses the documented midpoint convention.
    ties=[float(angular_density(.7,0.,0.,0.,0.,16.,g)) for g in [0.,.5,1.]]
    assert np.array_equal(ties,[.25]*3)
    report=dict(schema='blackhole-invariant-angular-validation/1',cases=len(rows),rows=rows,
        maximum_normalization_error=normerr,maximum_negative_probability_error=negerr,
        both_zero_spin_midpoint_density=ties,
        model_sha256=hashlib.sha256((BASE/'invariant_model.py').read_bytes()).hexdigest(),
        scope='Independent domain-split48-point Gaussian quadrature; validates the normalized angular model, not invariance of the physical PN measure, event likelihoods, or campaign coverage.')
    (BASE/'angular-validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['cases','maximum_normalization_error','maximum_negative_probability_error']},indent=2))


if __name__=='__main__':main()
