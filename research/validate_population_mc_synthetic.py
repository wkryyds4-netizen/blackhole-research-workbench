"""Small, independent calibration of paired empirical bootstrap mathematics."""
import json
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def exact_covariance(weights,total):
    mean=weights.sum(axis=0)/total
    return weights.T@weights/total**2-np.outer(mean,mean)/total


def draw_means(rng,weights,total,repeats):
    m=len(weights)
    probs=np.r_[np.full(m,1/total),1-m/total]
    output=[]
    for start in range(0,repeats,500):
        counts=rng.multinomial(total,probs,size=min(500,repeats-start))
        assert np.all(counts.sum(axis=1)==total)
        output.append(counts[:,:m]@weights/total)
    return np.concatenate(output)


def compare_covariance(draws,truth,relative_allowance=0):
    empirical=np.cov(draws,rowvar=False,ddof=1)
    blocks=np.array([np.cov(part,rowvar=False,ddof=1) for part in np.array_split(draws,20)])
    se=blocks.std(axis=0,ddof=1)/np.sqrt(len(blocks))
    allowed=6*se+relative_allowance*np.sqrt(np.outer(np.diag(truth),np.diag(truth)))+1e-14
    assert np.all(np.abs(empirical-truth)<allowed),(empirical,truth,allowed)
    return {'empirical':empirical.tolist(),'predicted':truth.tolist(),
            'max_absolute_error':float(np.max(abs(empirical-truth))),
            'max_error_over_allowed':float(np.max(abs(empirical-truth)/allowed))}


def run():
    rng=np.random.default_rng(20261024)
    repeats=20000
    t=np.linspace(0,2*np.pi,400,endpoint=False)
    a=1+.25*np.sin(t)
    weights=np.column_stack((a,1+.25*np.cos(t),1-.25*np.sin(t),3*a))
    event=draw_means(rng,weights,400,repeats)
    injection=draw_means(rng,weights,2000,repeats)
    ce=exact_covariance(weights,400);ci=exact_covariance(weights,2000)
    raw_event=compare_covariance(event,ce)
    raw_injection=compare_covariance(injection,ci)
    mu=weights.mean(axis=0);xi=weights.sum(axis=0)/2000
    # K=3 independent events, and one common injection campaign.
    e2=draw_means(rng,weights,400,repeats)
    e3=draw_means(rng,weights,400,repeats)
    logl=np.log(event)+np.log(e2)+np.log(e3)-3*np.log(injection)
    delta_cov=3*ce/np.outer(mu,mu)+9*ci/np.outer(xi,xi)
    delta_check=compare_covariance(logl,delta_cov,relative_allowance=.03)
    # Proportional columns have identical normalized population likelihood.
    cancellation=float(np.max(abs(logl[:,3]-logl[:,0])))
    assert cancellation<1e-12
    # Explicit full-zero vector vs collapsed zero bucket, both exact covariance.
    small=np.array([[1.,2.],[2.,1.],[3.,5.]])
    full=np.vstack((small,np.zeros((17,2))))
    full_draws=draw_means(rng,full,20,repeats)
    collapsed_draws=draw_means(rng,small,20,repeats)
    expected=exact_covariance(full,20)
    assert np.allclose(expected,exact_covariance(small,20),rtol=0,atol=1e-16)
    explicit_check=compare_covariance(full_draws,expected)
    collapsed_check=compare_covariance(collapsed_draws,expected)
    # Duplicating each observed row creates no independent draws. Treating the
    # copies as n*k IID observations spuriously divides mean covariance by k.
    copies=7
    original=np.array([[.1,3.],[.4,2.],[1.,1.],[3.,.4],[9.,.1]])
    duplicated=np.repeat(original,copies,axis=0)
    original_cov=exact_covariance(original,len(original))
    naive_cov=exact_covariance(duplicated,len(duplicated))
    assert np.allclose(naive_cov,original_cov/copies,rtol=1e-14,atol=1e-14)
    duplicate_rng=np.random.default_rng(2610247)
    original_counts=duplicate_rng.multinomial(len(original),np.full(len(original),1/len(original)),size=repeats)
    original_mean=original_counts@original/len(original)
    # Correct cluster resampling repeats every member of an original row's
    # perfect-correlation cluster together, using the same original count.
    cluster_mean=np.repeat(original_counts,copies,axis=1)@duplicated/len(duplicated)
    cluster_cancellation=float(np.max(abs(cluster_mean-original_mean)))
    assert cluster_cancellation<1e-13
    cluster_cov=np.cov(cluster_mean,rowvar=False,ddof=1)
    assert np.allclose(cluster_cov,np.cov(original_mean,rowvar=False,ddof=1),rtol=1e-14,atol=1e-14)
    naive_draws=draw_means(duplicate_rng,duplicated,len(duplicated),repeats)
    naive_check=compare_covariance(naive_draws,naive_cov)
    duplicate_check={'copies_per_original_row':copies,'original_rows':len(original),
                     'naive_exact_covariance_ratio':(np.diag(naive_cov)/np.diag(original_cov)).tolist(),
                     'cluster_variance_ratio':(np.diag(cluster_cov)/np.diag(np.cov(original_mean,rowvar=False,ddof=1))).tolist(),
                     'same_draw_cluster_mean_max_error':cluster_cancellation,
                     'naive_empirical_covariance_validation':naive_check,
                     'interpretation':'The 1/7 variance reduction is spurious when copied rows are mislabeled independent. Correct perfectly correlated cluster replication leaves every bootstrap mean unchanged.'}
    report={'status':'pass','date':'2026-10-24','seed':20261024,'replicates':repeats,
            'event_exact_covariance':raw_event,'injection_exact_covariance':raw_injection,
            'log_likelihood_delta_covariance':delta_check,'proportional_model_cancellation_max':cancellation,
            'explicit_zero_vector_covariance':explicit_check,'collapsed_zero_bucket_covariance':collapsed_check,
            'duplicate_rows_do_not_create_information':duplicate_check,
            'scope':'Independent synthetic mathematics validation only; does not validate unknown sampling design or missed tails in real inputs.',
            'tolerance':'Raw covariance: six block-based Monte Carlo standard errors. Log covariance: same plus 3% delta-approximation allowance.'}
    path=ROOT/'research/population-stage2/independent-synthetic-validation.json'
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({'status':'pass','replicates':repeats,'proportional_cancellation':cancellation,
                      'log_delta_max_error_over_allowed':delta_check['max_error_over_allowed']},indent=2))


if __name__=='__main__':run()
