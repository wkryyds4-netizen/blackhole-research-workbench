"""A second posterior algorithm, independent of branching_posterior.py."""
import numpy as np
from scipy.optimize import minimize
from scipy.special import ndtri
from scipy.stats import rankdata


def to_simplex(z):
    full=np.column_stack([z,np.zeros(len(z))]);full-=full.max(axis=1)[:,None]
    b=np.exp(full);return b/b.sum(axis=1)[:,None]


def target(z,means,efficiency):
    b=to_simplex(z)
    with np.errstate(divide='ignore',invalid='ignore'):
        # Uniform intrinsic Dirichlet prior, plus determinant of the log-ratio
        # map. Computing a plain likelihood in z would give a different prior.
        value=np.log(b@means.T).sum(axis=1)-len(means)*np.log(b@efficiency)+np.log(b).sum(axis=1)
    return np.where(np.isfinite(value),value,-np.inf)


def ranks(a):
    r=rankdata(a.ravel()).reshape(a.shape)
    return ndtri((r-.375)/(a.size+.25))


def basic_rhat(a):
    n=a.shape[1];within=np.var(a,axis=1,ddof=1).mean()
    between=n*np.var(a.mean(axis=1),ddof=1)
    return float(np.sqrt(((n-1)*within/n+between/n)/within))


def rhat(a):
    n=a.shape[1]//2;split=np.concatenate([a[:,:n],a[:,-n:]],axis=0)
    return max(basic_rhat(ranks(split)),basic_rhat(ranks(abs(split-np.median(split)))))


def ess(a):
    m,n=a.shape
    centered=a-a.mean(axis=1)[:,None]
    transformed=np.fft.rfft(centered,n=2*n,axis=1)
    covariance=np.fft.irfft(transformed*np.conj(transformed),n=2*n,axis=1)[:,:n]/n
    within=np.var(a,axis=1,ddof=1).mean()
    varplus=(n-1)*within/n+np.var(a.mean(axis=1),ddof=1)
    rho=1-(within-covariance.mean(axis=0))/varplus;rho[0]=1
    pair_sum=0.;previous=np.inf
    for i in range(0,n-1,2):
        pair=float(rho[i]+rho[i+1])
        if pair<0:break
        pair=min(pair,previous);pair_sum+=pair;previous=pair
    tau=max(1.,2*pair_sum-1)
    return float(min(m*n,m*n/tau))


def sample(means,efficiency,design,seed=None):
    rng=np.random.default_rng(design['seed'] if seed is None else seed)
    means=np.asarray(means,float);efficiency=np.asarray(efficiency,float)
    means=means/means.max(axis=1)[:,None]
    def objective(z):return -target(np.asarray(z)[None,:],means,efficiency)[0]
    candidates=[]
    for b in [np.ones(5)/5,*rng.dirichlet(np.ones(5),size=2)]:
        result=minimize(objective,np.log(b[:4]/b[4]),method='BFGS',options={'gtol':1e-7,'maxiter':1000})
        assert np.isfinite(result.fun)
        candidates.append(result)
    mode=min(candidates,key=lambda r:r.fun).x
    h=1e-3;eye=np.eye(4)*h;H=np.zeros((4,4));center=objective(mode)
    for i in range(4):
        H[i,i]=(objective(mode+eye[i])-2*center+objective(mode-eye[i]))/(h*h)
        for j in range(i):
            H[i,j]=H[j,i]=(objective(mode+eye[i]+eye[j])-objective(mode+eye[i]-eye[j])-objective(mode-eye[i]+eye[j])+objective(mode-eye[i]-eye[j]))/(4*h*h)
    eigen,vectors=np.linalg.eigh(H)
    # This only defines a proposal. Positive regularization does not alter
    # the exact target or reject any posterior point.
    root=vectors@np.diag(1/np.sqrt(np.maximum(eigen,1e-5)))
    chains=design['chains'];z=mode+rng.normal(size=(chains,4))@root.T*2
    current=target(z,means,efficiency);scale=2.38/2
    warm=design['warmup_steps'];production=design['production_steps'];thin=design['save_every']
    draws=[];accepted=np.zeros(chains,int);block=0
    for step in range(warm+production):
        proposed=z+rng.normal(size=z.shape)@root.T*scale
        values=target(proposed,means,efficiency)
        take=np.log(rng.uniform(size=chains))<values-current
        z[take]=proposed[take];current[take]=values[take]
        if step<warm:
            block+=int(take.sum())
            if (step+1)%100==0:
                scale*=np.exp(block/(100*chains)-.25);block=0
        else:
            accepted+=take
            if (step-warm+1)%thin==0:draws.append(to_simplex(z))
    draws=np.transpose(np.array(draws),(1,0,2))
    flat=draws.reshape((-1,5));checks=[]
    for i in range(5):
        a=draws[:,:,i];lo,hi=np.quantile(a,[.05,.95])
        checks.append(dict(component=i,rank_folded_split_rhat=rhat(a),bulk_ess=ess(ranks(a)),tail_ess=min(ess((a<=lo).astype(float)),ess((a>=hi).astype(float)))))
    return dict(intrinsic_quantiles=np.quantile(flat,[.05,.5,.95],axis=0).T.tolist(),intrinsic_mean=flat.mean(axis=0).tolist(),
                chains=chains,saved_per_chain=draws.shape[1],acceptance_fraction=(accepted/production).tolist(),proposal_scale=float(scale),
                minimum_hessian_eigenvalue=float(eigen.min()),diagnostics=checks,
                maximum_rhat=max(r['rank_folded_split_rhat'] for r in checks),minimum_bulk_ess=min(r['bulk_ess'] for r in checks),minimum_tail_ess=min(r['tail_ess'] for r in checks))
