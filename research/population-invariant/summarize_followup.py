"""Collect verified follow-up evidence; never refit or overwrite the baseline."""
from pathlib import Path
import hashlib
import json

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
FROZEN={
    'results/model.json':'6dcf3c9f4fd51b2007977d131c61f69c7d861afb986551ea82b70a7667a51083',
    'data/posteriors.npz':'3aec709780ac373e485682323c20cbe6ed0c14028c66b030bb73d848739489c4',
    'data/injections.npz':'00ce46283aaf5c1e3df1a8b479a2f3c5bf904e4aef6aa0cf641f93df54bf7a13'}


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    for name,digest in FROZEN.items():assert sha(ROOT/name)==digest,name
    reports={}
    for stem in ['control','paired-o3','paired-pe','support']:
        result_path=BASE/f'{stem}-results.json'
        validation=json.loads((BASE/f'{stem}-validation.json').read_text())
        assert validation['result_sha256']==sha(result_path)
        reports[stem]=json.loads(result_path.read_text())
    protocol=hashlib.sha256(b'invariant-control/1'+(BASE/'invariant_model.py').read_bytes()+
        (BASE/'run_control.py').read_bytes()+(ROOT/'model.py').read_bytes()+
        (ROOT/'research/population-stage5/run_tail_study.py').read_bytes()).hexdigest()
    assert json.loads((BASE/'event-integrals.json').read_text())['protocol']==protocol
    pe=reports['paired-pe'];support=reports['support']
    assert pe['input_hashes']['paired_o3']==sha(BASE/'paired-o3-results.json')
    assert pe['input_hashes']['paired_o3_validation']==sha(BASE/'paired-o3-validation.json')
    assert pe['input_hashes']['design']==sha(BASE/'paired-pe-design.json')
    plan=json.loads((BASE/'tail-validation-plan.json').read_text())
    assert plan['support_results_sha256']==sha(BASE/'support-results.json')
    assert plan['support_validation_sha256']==sha(BASE/'support-validation.json')
    assert plan['producer_sha256']==sha(BASE/'plan_tail_validation.py')
    derivative=json.loads((BASE/'paired-pe-validation.json').read_text())
    assert derivative['original_step_derivative_gate_passed'] is False
    files=['paired-o3-design.json','paired_o3.py','paired-o3-results.json','validate_paired_o3.py','paired-o3-validation.json',
        'paired-pe-design.json','paired_pe.py','paired-pe-results.json','validate_paired_pe.py','paired-pe-validation.json',
        'diagnose_pe_derivative.py','pe-derivative-diagnosis.json','plot_finite_banks.py','finite-bank-comparison.svg',
        'support-design.json','support_sensitivity.py','support-results.json','validate_support.py','support-validation.json',
        'plot_support.py','support-sensitivity.svg','plan_tail_validation.py','tail-validation-plan.json','summarize_followup.py']
    milestone=dict(schema='blackhole-full69-finite-bank-support-milestone/1',date='2026-10-25',
        events=69,full_source_rows=reports['control']['full_source_rows'],replicates=512,
        matched_contrasts=[row for row in pe['contrasts'] if row['name'].startswith('matched:')],
        matched_support_ordering_thresholds=[row for row in support['ordering_certificate_thresholds'] if row['policy']==0],
        completion_event=support['completion_event'],
        maximum_matched_infinite_completion_shift_pp=max(support['event_completions'][-1]['median_shift_pp'][0]),
        known_omitted_o3_redshift_mass=support['known_o3_omitted_domain'][0]['population_redshift_tail_mass'],
        tail_planning_direct_draws=plan['total_direct_draws_if_all_zero'],
        original_derivative_step_passed=False,derivative_failure_and_followup_preserved=True,
        file_sha256={name:sha(BASE/name) for name in files},frozen_hashes=FROZEN,event_integral_protocol=protocol,
        completed=['all69 paired O3/event/joint finite-bank experiment and independent validation',
            'conditional log-shape envelopes, single-zero-event integral-completion family and independent numerical validation',
            'known omitted O3 domain upper bound and unexecuted direct-draw resource plan',
            'two scientific figures generated and visually checked; research documentation updated'],
        still_required=['actual missed-likelihood/selection evidence is not recovered by these conditional diagnostics',
            'original tilt common reference and historical injection execution metadata remain unresolved',
            'coherent bilingual LaTeX/PDF rewrite, workbench adaptation and publication/reproduction/license audit'],
        scope='Conditional sensitivity results, not real-campaign/PE-chain coverage or formation-channel identification.',
        new_git_publication=False,manuscript_updated_to_this_milestone=False)
    (BASE/'finite-bank-milestone.json').write_text(json.dumps(milestone,indent=2),encoding='utf-8')
    print(json.dumps({k:milestone[k] for k in ['events','full_source_rows','replicates','maximum_matched_infinite_completion_shift_pp','tail_planning_direct_draws','original_derivative_step_passed','new_git_publication']},indent=2))


if __name__=='__main__':main()
