"""Explicit4D/7D measure bridge with declared hidden-spin conditionals.

No physical-channel prediction is added by this mathematical interface. It
makes the extra conditional distribution required by a4D model explicit.
"""
import numpy as np
from scipy.integrate import quad


def chirp_to_primary_jacobian(q):
    q=np.asarray(q,float)
    assert np.all((q>0)&(q<=1))
    return (1+q)**.2/q**.6


def projected_spin_z(s,amax,family):
    value=abs(float(s))
    if value>=amax:return 0.
    if family=='uniform_magnitude':
        return np.inf if value==0 else np.log(amax/value)/(2*amax)
    if family=='beta2_5_magnitude':return 3/amax*(1-value/amax)**5
    raise ValueError('Unknown hidden-spin conditional family')


def effective_spin_prior(x,q,amax=.99,family='uniform_magnitude'):
    assert 0<q<=1 and amax>0
    x=float(x)
    if abs(x)>=amax:return 0.
    if x==0 and family=='uniform_magnitude':return (1+q)*(2-np.log(q))/(2*amax)
    b=(1+q)*x
    low=max(-amax,b-q*amax);high=min(amax,b+q*amax)
    if low>=high:return 0.
    points=sorted(set(t for t in [0.,b] if low<t<high))
    f=lambda s:projected_spin_z(s,amax,family)*projected_spin_z((b-s)/q,amax,family)
    value,error=quad(f,low,high,points=points,epsabs=2e-9,epsrel=2e-9,limit=160)
    assert np.isfinite(value) and value>=0
    return (1+q)/q*value


def uniform_reduced_prior(prior7,q,x,amax=.99):
    """Given actual/scaled pi7=pi_mass,z*UniformSpin4, return pi4.

    Per-event constants omitted from pi7 stay omitted from pi4. The Jacobian
    is d m1/d chirp_mass, and f_chi is conditional on q.
    """
    return prior7*(4*amax**2)*chirp_to_primary_jacobian(q)*effective_spin_prior(x,q,amax)


def hidden_spin_density(a1,a2,c1,c2,amax,family):
    if not (0<a1<amax and 0<a2<amax and -1<=c1<=1 and -1<=c2<=1):return 0.
    if family=='uniform_magnitude':return 1/(4*amax**2)
    if family=='beta2_5_magnitude':
        radial=lambda a:30*(a/amax)*(1-a/amax)**4/amax
        return radial(a1)*radial(a2)/4
    raise ValueError('Unknown hidden-spin conditional family')


def lift_channel_density(density4,q,a1,a2,c1,c2,amax=.99,family='uniform_magnitude'):
    x=(a1*c1+q*a2*c2)/(1+q)
    spin=hidden_spin_density(a1,a2,c1,c2,amax,family)
    if spin==0:return 0.
    marginal=effective_spin_prior(x,q,amax,family)
    assert marginal>0
    return density4/chirp_to_primary_jacobian(q)*spin/marginal
