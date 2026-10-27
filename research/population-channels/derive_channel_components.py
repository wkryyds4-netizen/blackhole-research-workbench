"""Separate physical-model and common-regularizer terms of a fitted mixture.

The uniform regularizer is a mathematical density term, not a measured extra
astrophysical population. Original detected_quantiles label regularized channels.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import hashlib
import json
import numpy as np
from branching_posterior import sample,weighted_quantiles

BASE=Path(__file__).resolve().parent
def sha(p):
    with p.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def save(p,d):
    tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8');tmp.replace(p)


def main():
    fits=json.loads((BASE/'branching-results.json').read_text());assert len(fits['cases'])==24
    design=json.loads((BASE/'likelihood-design.json').read_text())
    sampler=json.loads((BASE/'branching-sampler-design.json').read_text())
    events=json.loads((BASE/'channel-log-event-integrals.json').read_text())['events']
    selection=json.loads((BASE/'channel-log-selection-integrals.json').read_text())
    normalizers=json.loads((BASE/'channel-log-normalization.json').read_text())
    norm=np.array([next(r['fraction'] for r in normalizers['channels'] if r['channel']==c) for c in design['channels']])
    report=dict(schema='blackhole-channel-detected-component-decomposition/1',cases={},
                quantities=dict(regularized_label='Original detected_quantiles: beta_c*((1-eps)*xi_c+eps*xi_U)/xi_total. Common U is allocated to channel labels according to beta; not an independent attribution by the data.',
                                model_term='(1-eps)*beta_c*xi_c/xi_total; these five fractions plus uniform_term sum to1.',
                                uniform_term='eps*xi_U/xi_total; mathematical regularizer responsibility, not a physical outlier fraction.',
                                model_term_conditional='beta_c*xi_c/sum_d(beta_d*xi_d); normalized within the non-regularizer model term.'),
                input_sha256={p.name:sha(p) for p in [Path(__file__),BASE/'branching_posterior.py',BASE/'branching-results.json',BASE/'channel-log-integrals-complete.json',BASE/'branching-sampler-design.json']})
    path=BASE/'detected-component-results.json'
    for key,case in fits['cases'].items():
        fi=design['hidden_spin_conditionals'].index(case['family']);eps=case['epsilon']
        raw=np.array([r['mean'][fi] for r in events])
        A0=raw[:,:5]/norm[None,:];U=raw[:,5]
        rs=np.array(selection['early']['total'])[fi]+np.array(selection['banks'][case['bank']]['policies'][case['policy']]['total'])[fi]
        xi0=rs[:5]/norm;xiU=rs[5]
        means=(1-eps)*A0+eps*U[:,None];efficiency=(1-eps)*xi0+eps*xiU
        reps=[]
        for seed,original in zip(sampler['replicate_seeds'],case['replicates']):
            r=sample(means,efficiency,sampler['pilot_power'],sampler['final_power'],seed,keep=True)
            assert np.max(abs(np.array(r['intrinsic_quantiles'])-original['intrinsic_quantiles']))<1e-9
            b=r['samples'];w=r['weights'];total=b@efficiency
            model=(1-eps)*b*xi0[None,:]/total[:,None]
            uniform=eps*xiU/total
            conditional=b*xi0[None,:];conditional/=conditional.sum(axis=1)[:,None]
            label=b*efficiency[None,:]/total[:,None]
            assert np.max(abs(model.sum(axis=1)+uniform-1))<1e-12
            assert np.max(abs(label-model-b*uniform[:,None]))<1e-12
            reps.append(dict(seed=seed,model_term_quantiles=weighted_quantiles(model,w).tolist(),
                             uniform_term_quantiles=weighted_quantiles(uniform[:,None],w)[0].tolist(),
                             model_term_conditional_quantiles=weighted_quantiles(conditional,w).tolist()))
        beta=np.mean([r['intrinsic_mean'] for r in case['replicates']],axis=0)
        event_uniform=eps*U/(means@beta)
        largest=np.argsort(event_uniform)[-5:][::-1]
        report['cases'][key]=dict(model_term_quantiles=np.mean([r['model_term_quantiles'] for r in reps],axis=0).tolist(),
            uniform_term_quantiles=np.mean([r['uniform_term_quantiles'] for r in reps],axis=0).tolist(),
            model_term_conditional_quantiles=np.mean([r['model_term_conditional_quantiles'] for r in reps],axis=0).tolist(),
            replicate_maximum_model_term_quantile_span=float(np.ptp(np.array([r['model_term_quantiles'] for r in reps]),axis=0).max()),replicates=reps,
            posterior_mean_point_diagnostic=dict(uniform_selection_fraction=float(eps*xiU/(efficiency@beta)),
                summed_event_uniform_responsibility=float(event_uniform.sum()),
                largest_event_uniform_responsibilities=[dict(event=events[i]['event'],index=events[i]['index'],responsibility=float(event_uniform[i])) for i in largest],
                note='Responsibilities evaluated at posterior mean beta; not posterior-averaged event classifications.'),
            uniform_efficiency=float(xiU),unregularized_channel_efficiency=xi0.tolist())
        save(path,report);print(f'Decomposed detected labels: {key}',flush=True)
    report.update(status='computed',scope='Algebraic decomposition of24 conditional fitted models. No additional physical channel, population membership or missing-support estimate is measured.')
    save(path,report)


if __name__=='__main__':main()
