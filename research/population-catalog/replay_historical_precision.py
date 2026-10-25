"""Evaluate unmodified published 7.0 Python tilt routines in isolation.

Only lal's three physical constants are provided, from the same fixed header.
This is a controlled Python-algorithm replay, not a reconstructed LAL install.
No copied upstream code is published by this adapter.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import re
import sys
import types
import warnings
import numpy as np

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
SOURCES=ROOT/'research/population-selection-audit/sources'
PREFIX='gitlab-1438-af4bc21a9f92-'


def load_historical():
    manifest={r['file']:r for r in json.loads((SOURCES/'manifest.json').read_text())}
    def verified(suffix):
        path=SOURCES/(PREFIX+suffix);assert hashlib.sha256(path.read_bytes()).hexdigest()==manifest[path.name]['sha256'];return path
    header=verified('lal--lib--std--LALConstants.h').read_text();constants={}
    for name in ['G_SI','C_SI','MSUN_SI','MTSUN_SI']:
        match=re.search(r'^#define LAL_'+name+r'\s+([\d.eE+\-]+)',header,re.M);assert match;constants[name]=float(match[1])
    assert 'lal' not in sys.modules,'Use the isolated runtime without an installed/imported lal module'
    lal=types.ModuleType('lal')
    for name in ['G_SI','C_SI','MSUN_SI']:setattr(lal,name,constants[name])
    sys.modules['lal']=lal
    package=types.ModuleType('historical_laltilts');package.__path__=[];sys.modules[package.__name__]=package
    loaded={}
    for short in ['tilts_at_infinity_utils','calc_tilts_prec_avg_regularized']:
        path=verified('lalsimulation--python--lalsimulation--tilts_at_infinity--'+short+'.py')
        name=package.__name__+'.'+short;spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module);loaded[short]=module
    return loaded['calc_tilts_prec_avg_regularized'].prec_avg_tilt_comp,constants


def main():
    function,constants=load_historical()
    sys.path.insert(0,str(ROOT/'research/population-reference'));from transport import to_infinity
    inputs=json.loads((BASE/'alias-heldout-input.json').read_text());validation=json.loads((BASE/'alias-heldout-validation.json').read_text())
    worst=max(validation['rows'],key=lambda r:max(map(abs,r['difference'])))['source_row']
    near_equal=next(r['source_row'] for r in inputs['rows'] if r['reason']=='maximum q')
    rows=[]
    for index in sorted({worst,near_equal}):
        record=next(r for r in inputs['rows'] if r['source_row']==index);v=record['inputs'];q=v['q'];m1=v['total_detector_mass']/(1+q);m2=q*m1
        for frequency in [16.,20.]:
            target=record['canonical_cosines'] if frequency==16 else record['alias_cosines'];answers={};messages=[]
            for label,options in [('historical_defaults',{}),('historical_tight',dict(ode_atol_base=1e-12,ode_atol_floor=1e-15,ode_rtol=1e-13,lin_tol=1e-12,root_tol=1e-10,method='DOP853',use_solve_ivp=True,use_fallback=False))]:
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always')
                    result=function(m1*constants['MSUN_SI'],m2*constants['MSUN_SI'],v['a1'],v['a2'],v['theta1'],v['theta2'],v['phi12'],frequency,LPNorder=2,LPNspins=False,**options)
                cosines=np.cos([result['tilt1_inf'],result['tilt2_inf']]);assert np.isfinite(cosines).all()
                answers[label]=dict(cosines=cosines.tolist(),released_difference=(cosines-target).tolist())
                messages += [dict(method=label,message=str(w.message)) for w in caught if not issubclass(w.category,DeprecationWarning)]
            eta=q/(1+q)**2
            old_mtsun=constants['MSUN_SI']*constants['G_SI']/constants['C_SI']**3
            speed=(np.pi*v['total_detector_mass']*old_mtsun*frequency)**(1/3)
            correction=1+(1.5+eta/6)*speed**2+(27/8-19*eta/8+eta*eta/24)*speed**4
            ours=to_infinity(**dict(v,reference_hz=frequency),initial_separation=(correction/speed)**2,solver='DOP853',rtol=3e-14,atol=3e-16)
            oc=np.array([ours[f'cos_tilt_{i}_infinity'] for i in [1,2]])
            answers['independent_same_constants']=dict(cosines=oc.tolist(),released_difference=(oc-target).tolist())
            rows.append(dict(source_row=index,frequency=frequency,answers=answers,warnings=messages,
                historical_tight_minus_independent=(np.array(answers['historical_tight']['cosines'])-oc).tolist()))
            print(f'Completed source row {index}, {frequency:g} Hz',flush=True)
    summaries={name:float(max(max(map(abs,r['answers'][name]['released_difference'])) for r in rows)) for name in ['historical_defaults','historical_tight','independent_same_constants']}
    report=dict(schema='blackhole-historical-python-precision-replay/1',commit='af4bc21a9f9246ef5403052ec3b66ed2b52d90ed',constants=constants,
        runtime=dict(python=sys.version,numpy=np.__version__,dependency_scope='Current isolated NumPy/SciPy; lal module provides ONLY three constants from the pinned 7.0 header. Unmodified historical pure-Python routines loaded by file path; no hybrid or compiled LAL simulation.'),
        maximum_released_cosine_residual=summaries,maximum_tight_implementation_difference=float(max(max(map(abs,r['historical_tight_minus_independent'])) for r in rows)),rows=rows,
        original_numerical_checks_passed=validation['numerical_checks_passed'],original_gate_and_failure_preserved=True,
        original_validation_sha256=hashlib.sha256((BASE/'alias-heldout-validation.json').read_bytes()).hexdigest(),
        limits='A controlled algorithm/tolerance/unit-convention diagnosis on two source rows, not a complete historical execution environment or full-domain validation.')
    (BASE/'historical-precision-replay.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps({k:report[k] for k in ['maximum_released_cosine_residual','maximum_tight_implementation_difference','original_numerical_checks_passed']},indent=2))


if __name__=='__main__':main()
