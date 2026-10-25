"""Replay historical defaults on the previously fixed 70-row validation set.

This is a follow-up diagnostic, not a new held-out accuracy test. The original
release-agreement gate and its failure are retained independently.
"""
from pathlib import Path
import hashlib
import json
import platform
import warnings
import numpy as np
import scipy
from replay_historical_precision import load_historical

BASE=Path(__file__).resolve().parent


def main():
    function,constants=load_historical()
    inputs_path=BASE/'alias-heldout-input.json'
    digest=hashlib.sha256(inputs_path.read_bytes()).hexdigest()
    inputs=json.loads(inputs_path.read_text())
    original_path=BASE/'alias-heldout-validation.json'
    original=json.loads(original_path.read_text())
    assert original['input_sha256']==digest and len(inputs['rows'])==70
    checkpoint=BASE/'historical-defaults-checkpoint.json'
    protocol=dict(input_sha256=digest,commit='af4bc21a9f9246ef5403052ec3b66ed2b52d90ed',
        constants=constants,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        options=dict(LPNorder=2,LPNspins=False),scope='Unmodified historical Python routine; only lal physical constants supplied; current NumPy/SciPy runtime.')
    previous=json.loads(checkpoint.read_text()) if checkpoint.exists() else dict(protocol=protocol,rows=[])
    assert previous['protocol']==protocol,'Protocol changed: preserve and version previous checkpoint.'
    rows=previous['rows']
    for record in inputs['rows']:
        index=record['source_row'];v=record['inputs'];q=v['q']
        m1=v['total_detector_mass']/(1+q);m2=q*m1
        for frequency in [16.,20.]:
            if any(r['source_row']==index and r['frequency']==frequency for r in rows):continue
            target=record['canonical_cosines'] if frequency==16 else record['alias_cosines']
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                result=function(m1*constants['MSUN_SI'],m2*constants['MSUN_SI'],v['a1'],v['a2'],v['theta1'],v['theta2'],v['phi12'],frequency,LPNorder=2,LPNspins=False)
            cosines=np.cos([result['tilt1_inf'],result['tilt2_inf']]);assert np.isfinite(cosines).all()
            old=next(r for r in original['rows'] if r['source_row']==index and r['frequency']==frequency)
            rows.append(dict(source_row=index,reason=record['reason'],frequency=frequency,
                cosines=cosines.tolist(),difference=(cosines-target).tolist(),original_independent_difference=old['difference'],
                warnings=[str(w.message) for w in caught if not issubclass(w.category,DeprecationWarning)]))
            temporary=checkpoint.with_suffix('.tmp');temporary.write_text(json.dumps(dict(protocol=protocol,rows=rows),indent=2),encoding='utf-8');temporary.replace(checkpoint)
        if len(rows)%20==0:print(f'Historical-default comparisons completed: {len(rows)}/140',flush=True)
    assert len(rows)==140
    summaries=[]
    for frequency in [16.,20.]:
        subset=[r for r in rows if r['frequency']==frequency]
        values=np.array([r['difference'] for r in subset]);before=np.array([r['original_independent_difference'] for r in subset])
        summaries.append(dict(frequency=frequency,rows=len(subset),
            maximum_released_cosine_residual=float(np.max(abs(values))),
            median_row_maximum_released_cosine_residual=float(np.median(abs(values).max(axis=1))),
            original_independent_maximum_residual=float(np.max(abs(before))),
            rows_with_warnings=sum(bool(r['warnings']) for r in subset)))
    report=dict(schema='blackhole-historical-defaults-followup/1',protocol=protocol,summaries=summaries,rows=rows,
        original_validation_sha256=hashlib.sha256(original_path.read_bytes()).hexdigest(),
        original_numerical_checks_passed=original['numerical_checks_passed'],
        limits='Follow-up on the same fixed 70 validation rows after seeing one failure. Not a new held-out test, full historical execution reconstruction, whole-posterior transport, or an injection-reference determination.')
    (BASE/'historical-defaults-validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(summaries,indent=2),flush=True)


if __name__=='__main__':main()
