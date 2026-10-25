"""Translate an assumed selection budget into a future direct-draw design.

Counts are planning numbers for an ideal restored O3 detection experiment,
not results from an executed campaign or importance-weighted samples.
"""
from pathlib import Path
import hashlib
import json
import math
import numpy as np

BASE=Path(__file__).resolve().parent


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    r=json.loads((BASE/'support-results.json').read_text())
    v=json.loads((BASE/'support-validation.json').read_text())
    assert v['result_sha256']==sha(BASE/'support-results.json')
    # 0.5% is a declared planning tolerance, not a measured truncation error.
    relative_budget=.005;alpha_family=.05;comparisons=12;alpha=alpha_family/comparisons
    with np.load(BASE/'control-checkpoints.npz',allow_pickle=False) as cp:
        matched=cp['selection_matched'];union=cp['selection_union']
    rows=[]
    for model,key in enumerate(r['keys']):
        domain=r['known_o3_omitted_domain'][model//3]
        mass=domain['known_omitted_union_mass_in_coded_selection_measure']
        minimum=float(min(matched[model].min(),union[model].min()))
        cap=relative_budget*minimum/(r['o3_time_fraction']*mass)
        n=math.ceil(math.log(alpha)/math.log1p(-cap))
        upper=-math.expm1(math.log(alpha)/n)
        assert upper<=cap
        rows.append(dict(model=model,key=key,sign_components=['negative','positive'],
            maximum_component_averaged_tail_detection_probability=cap,
            minimum_independent_direct_draws_per_component_if_zero_detected=n,
            zero_detection_binomial_upper_bound=upper))
    plan=dict(schema='blackhole-tail-validation-planning/1',timing='Post-result resource planning, not a preregistered or executed injection campaign.',
        relative_selection_completion_budget=relative_budget,
        resulting_selection_log_shape_budget=69*math.log1p(relative_budget),
        minimum_two_policy_amplitude_ordering_B_threshold=min(x['critical_log_correction_oscillation'] for x in r['ordering_certificate_thresholds']),
        family_confidence=1-alpha_family,component_confidence=1-alpha,component_count=comparisons,
        design='For each of6 invariant models and each sign component, draw iid sources from its normalized distribution conditional on the known omitted O3 domain. Evaluate the restored O3 union search; its probability upper bound also bounds the matched-filter subset. Require both component bounds at once. Counts below assume zero detections and use exact one-sided binomial bounds with Bonferroni over12 components.',
        model_requirements=rows,total_direct_draws_if_all_zero=2*sum(x['minimum_independent_direct_draws_per_component_if_zero_detected'] for x in rows),
        required_before_execution=['Recover or otherwise validate the actual O3 reference/measurement and detection protocol for the quantity to be bounded; a new arbitrary protocol cannot bound the historical one.',
            'Validate target direct sampling, including conditional sign/even-shape density, cosmology, redshift/spin union support, extrinsic averaging and time/noise sampling.',
            'Record a fixed stopping rule and pipeline thresholds. Nonzero detections require updated exact binomial bounds, not declaring success from zero-count planning.',
            'Keep represented-bank MC uncertainty and all event-integral omissions in a separate budget. The observed zero event component does not satisfy a finite uniform relative-completion bound if its missing component is positive.'],
        limitations='Planning counts apply to independent direct target draws with the correct historical detection law. Weighted injection reuse, correlated noise/time draws or unknown search configurations do not satisfy this binomial design automatically. Meeting selection-only tolerance would not establish total inference robustness.',
        support_results_sha256=sha(BASE/'support-results.json'),support_validation_sha256=sha(BASE/'support-validation.json'),producer_sha256=sha(Path(__file__)))
    (BASE/'tail-validation-plan.json').write_text(json.dumps(plan,indent=2),encoding='utf-8')
    print(json.dumps({k:plan[k] for k in ['relative_selection_completion_budget','resulting_selection_log_shape_budget','total_direct_draws_if_all_zero']},indent=2))


if __name__=='__main__':main()
