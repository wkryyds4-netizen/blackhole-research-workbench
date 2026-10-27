"""Importance-corrected five-channel simplex posterior; no raw-data IO."""
import numpy as np
from scipy.special import gammaincinv,gammaln,logsumexp
from scipy.stats import qmc


def log_likelihood(beta,means,efficiency,chunk=16384):
    beta=np.asarray(beta,float);means=np.asarray(means,float);efficiency=np.asarray(efficiency,float)
    assert beta.ndim==2 and beta.shape[1]==5 and means.ndim==2 and means.shape[1]==5
    assert np.all(beta>0) and np.all(means>=0) and np.all(efficiency>0)
    output=[]
    for start in range(0,len(beta),chunk):
        b=beta[start:start+chunk]
        with np.errstate(divide='ignore'):
            output.append(np.log(b@means.T).sum(axis=1)-len(means)*np.log(b@efficiency))
    return np.concatenate(output)


def dirichlet_logpdf(beta,alpha):
    return gammaln(np.sum(alpha))-np.sum(gammaln(alpha))+np.log(beta)@(np.asarray(alpha)-1)


def dirichlet_draws(u,alpha):
    # Scrambled Sobol points are in(0,1); guard endpoints explicitly rather
    # than allow gamma(0) or gamma(1) to create invalid simplex points.
    u=np.clip(u,np.nextafter(0.,1.),np.nextafter(1.,0.))
    g=gammaincinv(np.asarray(alpha)[None,:],u)
    return g/g.sum(axis=1)[:,None]


def weighted_quantiles(values,weights):
    out=[]
    for v in np.asarray(values).T:
        order=np.argsort(v);w=weights[order];c=np.cumsum(w)-.5*w
        out.append(np.interp([.05,.5,.95],c/c[-1],v[order]))
    return np.array(out)


def sample(means,efficiency,pilot_power=16,final_power=18,seed=20261101,keep=False):
    means=np.asarray(means,float);efficiency=np.asarray(efficiency,float)
    assert np.all(means.max(axis=1)>0)
    # This known scale is restored in the evidence; it cancels in posterior weights.
    scale=means.max(axis=1);m=means/scale[:,None]
    pilot=dirichlet_draws(qmc.Sobol(5,scramble=True,seed=seed+100).random_base2(pilot_power),np.ones(5))
    ll=log_likelihood(pilot,m,efficiency);w=np.exp(ll-logsumexp(ll))
    center=w@pilot;var=w@((pilot-center)**2)
    concentration=float(np.clip(np.median(center*(1-center)/np.maximum(var,1e-20)-1),2,500))
    alpha=np.maximum(center*concentration,.2)
    u=qmc.Sobol(6,scramble=True,seed=seed).random_base2(final_power)
    prior=u[:,0]<.2;beta=np.empty((len(u),5))
    beta[prior]=dirichlet_draws(u[prior,1:],np.ones(5))
    beta[~prior]=dirichlet_draws(u[~prior,1:],alpha)
    proposal=np.logaddexp(np.log(.2)+gammaln(5),np.log(.8)+dirichlet_logpdf(beta,alpha))
    ll=log_likelihood(beta,m,efficiency)
    logs=ll+gammaln(5)-proposal
    weights=np.exp(logs-logsumexp(logs))
    detected=beta*efficiency[None,:];detected/=detected.sum(axis=1)[:,None]
    out=dict(seed=seed,draws=len(beta),proposal_alpha=alpha.tolist(),pilot_concentration_ess=float(1/np.sum(w*w)),
             concentration_ess=float(1/np.sum(weights*weights)),maximum_weight_fraction=float(weights.max()),
             intrinsic_quantiles=weighted_quantiles(beta,weights).tolist(),detected_quantiles=weighted_quantiles(detected,weights).tolist(),
             intrinsic_mean=(weights@beta).tolist(),log_evidence=float(logsumexp(logs)-np.log(len(logs))+np.log(scale).sum()))
    if keep:out.update(samples=beta,weights=weights,log_sampling_density=proposal)
    return out
