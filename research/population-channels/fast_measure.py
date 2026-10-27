"""Vector positive convolutions for the explicitly declared spin conditionals.

This is numerical quadrature of a derived density, not interpolation of a fit.
The scalar independent implementation remains in projected_measure.py.
"""
from functools import lru_cache
import numpy as np
from numpy.polynomial.legendre import leggauss


@lru_cache(maxsize=8)
def rule(order, smoothed):
    u, w = leggauss(order)
    u = (u + 1) / 2
    w = w / 2
    if smoothed:
        denominator = u**4 + (1-u)**4
        t = u**4 / denominator
        reverse = (1-u)**4 / denominator
        w = w * 4*u**3*(1-u)**3 / denominator**2
    else:
        t, reverse = u, 1-u
    return t, reverse, w


def effective_spin_density(x, q, amax=.99, family='uniform_magnitude', order=48, chunk=4096):
    """Return f_chi(x|q), zero outside |x|<amax; inputs broadcast.

    q is in (0,1]. The validated production domain is q>=.03. Computation
    is dimensionless until the final 1/amax scaling. No clipping of samples.
    """
    x, q = np.broadcast_arrays(np.asarray(x, float), np.asarray(q, float))
    assert np.isfinite(x).all() and np.isfinite(q).all()
    assert np.all((q>0)&(q<=1)) and amax>0
    shape=x.shape; x=x.ravel(); q=q.ravel()
    out=np.zeros(len(x)); take=np.flatnonzero(abs(x)<amax)
    assert family in ['uniform_magnitude','beta2_5_magnitude']
    smoothed=family=='uniform_magnitude'
    t, reverse, weights=rule(order if smoothed else 6, smoothed)
    for start in range(0,len(take),chunk):
        ix=take[start:start+chunk]; qr=q[ix]; b=(1+qr)*abs(x[ix])/amax
        lo=np.maximum(-1., b-qr); hi=np.minimum(1., b+qr)
        cuts=np.stack([lo, np.clip(np.zeros(len(b)),lo,hi), np.clip(b,lo,hi), hi],axis=1)
        result=np.zeros(len(ix))
        for j in range(3):
            left=cuts[:,j,None]; right=cuts[:,j+1,None]; length=right-left
            # Distances to singular endpoints are evaluated without subtraction
            # of nearly equal numbers. Each interval has no internal singularity.
            s=left*reverse+right*t
            a=np.where(left==0,length*t,np.where(right==0,length*reverse,abs(s)))
            c=np.where(left==b[:,None],length*t,np.where(right==b[:,None],length*reverse,abs(b[:,None]-s))) / qr[:,None]
            if smoothed:
                value=(-np.log(np.maximum(a,np.finfo(float).tiny)))*(-np.log(np.maximum(c,np.finfo(float).tiny)))/4
            else:
                value=9*(1-a)**5*(1-c)**5
            result += (length[:,0] * np.sum(value*weights,axis=1))
        out[ix]=(1+qr)/qr*result/amax
    assert np.isfinite(out).all() and np.all(out>=0)
    return out.reshape(shape)


def hidden_spin_ratio(a1,a2,amax=.99,family='uniform_magnitude'):
    """Q_family/Q_uniform for isotropic directions, including radial support."""
    a1,a2=np.broadcast_arrays(np.asarray(a1,float),np.asarray(a2,float))
    support=(a1>0)&(a1<amax)&(a2>0)&(a2<amax)
    if family=='uniform_magnitude':value=np.ones(a1.shape)
    elif family=='beta2_5_magnitude':
        value=(30*(a1/amax)*(1-a1/amax)**4)*(30*(a2/amax)*(1-a2/amax)**4)
    else:raise ValueError(family)
    return np.where(support,value,0.)
