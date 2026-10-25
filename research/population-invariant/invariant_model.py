"""Normalized sign-asymmetry controls in conserved effective spin.

The control assumes2PN conservation of effective spin and preservation of the
isotropic spin measure. It is neither full-GR transport nor a formation model.
"""
import numpy as np

MAGNITUDES=[(1.5,8.5),(2.,5.)]
LAMBDAS=[0.,4.,16.]
KEYS=[f'beta-{a:g}-{b:g}_lambda-{lam:g}' for a,b in MAGNITUDES for lam in LAMBDAS]


def effective_spin(q,a1,a2,c1,c2):
    return (a1*c1+q*a2*c2)/(1+q)


def symmetric_ratio(q,a1,a2,x,lam):
    """Relative to an isotropic spin-angle density, conditioned on magnitudes."""
    u=a1/(1+q);v=q*a2/(1+q)
    return (1+lam*x*x)/(1+lam*(u*u+v*v)/3)


def component_factors(q,a1,a2,c1,c2,lam):
    x=np.asarray(effective_spin(q,a1,a2,c1,c2))
    ratio=symmetric_ratio(q,a1,a2,x,lam)
    negative=(2*(x<0)+(x==0)).astype(float)
    positive=(2*(x>0)+(x==0)).astype(float)
    return ratio[...,None]*np.stack([negative,positive],axis=-1)


def angular_density(q,a1,a2,c1,c2,lam,g):
    factors=component_factors(q,a1,a2,c1,c2,lam)
    return .25*(g*factors[...,0]+(1-g)*factors[...,1])
