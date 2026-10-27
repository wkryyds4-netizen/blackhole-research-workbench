"""Explicit restricted-domain original lift of published four-dimensional flows."""
from pathlib import Path
import json
import numpy as np
from channel_flow import ChannelFlow,CHANNELS
from fast_measure import effective_spin_density,hidden_spin_ratio

BASE=Path(__file__).resolve().parent
FAMILIES=['uniform_magnitude','beta2_5_magnitude']
A=.99
MASS_AREA=95-3*np.log(20)
UNIFORM_VOLUME=MASS_AREA*1.9*(2*A)


def coordinates(d):
    m=np.asarray(d['mass_1'],float);q=np.asarray(d['mass_ratio'],float)
    x=(d['a_1']*d['cos_tilt_1']+q*d['a_2']*d['cos_tilt_2'])/(1+q)
    jac=(1+q)**.2/q**.6
    theta=np.stack([m/jac,q,x,d['redshift']],axis=1)
    support=(m>=5)&(m<=100)&(m*q>=3)&(q>0)&(q<=1)&(d['redshift']>0)&(d['redshift']<=1.9)&(abs(x)<A)
    return theta,jac,support


def physical_domain(theta):
    mc,q,x,z=theta.T
    jac=(1+q)**.2/q**.6;m=mc*jac
    return np.isfinite(theta).all(axis=1)&(m>=5)&(m<=100)&(m*q>=3)&(q>0)&(q<=1)&(z>0)&(z<=1.9)&(abs(x)<A)


class ConditionalChannels:
    def __init__(self):
        self.design=json.loads((BASE/'likelihood-design.json').read_text())
        self.models=[ChannelFlow(c) for c in CHANNELS]

    def density4(self,theta):
        values=[];birth=self.design['context']['birth_spin']
        for channel,model in zip(CHANNELS,self.models):
            context=[birth,np.log(self.design['context']['ce_efficiency'])] if channel=='CE' else [birth]
            values.append(np.exp(model.log_density(theta,context)))
        return np.stack(values,axis=1)

    def weights(self,d):
        """N x two hidden-spin laws x (five raw channels + normalized U).

        Raw channels still need their separately recorded domain normalizer.
        d['prior'] is pi7 in (m1,q,z,a1,a2,c1,c2), exact or event-scaled.
        """
        theta,jac,support=coordinates(d)
        support &= (d['a_1']>0)&(d['a_1']<A)&(d['a_2']>0)&(d['a_2']<A)
        support &= (abs(d['cos_tilt_1'])<=1)&(abs(d['cos_tilt_2'])<=1)
        output=np.zeros((len(theta),2,6));indices=np.flatnonzero(support)
        if not len(indices):return output
        selected=theta[indices];raw=self.density4(selected)
        for j,family in enumerate(FAMILIES):
            f=effective_spin_density(selected[:,2],selected[:,1],A,family)
            ratio=hidden_spin_ratio(d['a_1'][indices],d['a_2'][indices],A,family)
            factor=ratio/(4*A*A*f*d['prior'][indices])
            output[indices,j,:5]=raw/jac[indices,None]*factor[:,None]
            output[indices,j,5]=factor/UNIFORM_VOLUME
        assert np.isfinite(output).all() and np.all(output>=0)
        return output
