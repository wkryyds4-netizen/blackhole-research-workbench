"""Fixed independent scalar, moment and boundary gates for the fast bridge."""
from pathlib import Path
import hashlib
import json
import time
import warnings
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import IntegrationWarning
from projected_measure import effective_spin_prior,chirp_to_primary_jacobian,hidden_spin_density
from fast_measure import effective_spin_density,hidden_spin_ratio

BASE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    design=json.loads((BASE/'bridge-design.json').read_text());gates=design['gates']; A=.99
    rng=np.random.default_rng(20261029)
    qs=np.r_[np.geomspace(.03,.99,21),1.]
    stress=[]
    for q in qs:
        surfaces=[0.,q/(1+q),(1-q)/(1+q),1/(1+q),1.]
        stress.extend((q,float(np.clip(v+delta,0,1))) for v in surfaces for delta in [-1e-8,0,1e-8])
    q=np.r_[rng.uniform(.03,1,384),[s[0] for s in stress]]
    x=A*np.r_[rng.uniform(-1,1,384),[s[1] for s in stress]]
    reports=[]
    for family in design['families']:
        start=time.perf_counter();fast=effective_spin_density(x,q,A,family);elapsed=time.perf_counter()-start
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always',IntegrationWarning)
            reference=np.array([effective_spin_prior(v,z,A,family) for v,z in zip(x,q)])
        positive=reference>1e-12
        relative=float(np.max(abs(fast[positive]/reference[positive]-1)))
        absolute=float(np.max(abs(fast-reference)))
        doubled=effective_spin_density(x,q,A,family,order=96)
        order_error=float(np.max(abs(fast[positive]/doubled[positive]-1)))
        assert relative<gates['density_relative_where_reference_above_1e_minus_12'],(family,'relative',relative)
        assert absolute<gates['density_absolute_everywhere'],(family,'absolute',absolute)
        assert order_error<gates['uniform_order_doubling_relative'],(family,'order',order_error)
        assert np.array_equal(effective_spin_density(-x,q,A,family),fast)
        assert np.all(effective_spin_density(np.array([-2*A,-A,A,2*A]),.5,A,family)==0)
        nodes,w=leggauss(512);xp=A*(1+nodes)/2;checks=[]
        for qr in [.03,.1,.5,1.]:
            f=effective_spin_density(xp,qr,A,family)
            norm=A*np.dot(w,f);moment=A*np.dot(w,xp*xp*f)
            expected=A*A*(1+qr*qr)/(9 if family=='uniform_magnitude' else 28)/(1+qr)**2
            assert abs(norm-1)<gates['normalization_absolute'],(family,'normalization',qr,norm)
            assert abs(moment-expected)<gates['moment_absolute'],(family,'moment',qr,moment)
            checks.append(dict(q=qr,normalization=float(norm),second_moment=float(moment),analytic=float(expected)))
        # Both lifts checked against the independently implemented scalar prior.
        qi=rng.uniform(.03,1,128);a1=rng.uniform(.001,A,128);a2=rng.uniform(.001,A,128)
        c1=rng.uniform(-1,1,128);c2=rng.uniform(-1,1,128);xx=(a1*c1+qi*a2*c2)/(1+qi)
        marginal=effective_spin_density(xx,qi,A,family)
        p4=rng.uniform(.001,1,128);pi7=rng.uniform(.1,1,128)
        ratio=hidden_spin_ratio(a1,a2,A,family)
        fast_weight=p4/(pi7*4*A*A*chirp_to_primary_jacobian(qi))*ratio/marginal
        scalar_spin=np.array([hidden_spin_density(a,b,c,d,A,family) for a,b,c,d in zip(a1,a2,c1,c2)])
        equivalent=p4/chirp_to_primary_jacobian(qi)*scalar_spin/marginal/pi7
        identity=float(np.max(abs(fast_weight/equivalent-1)))
        assert identity<gates['lift_weight_identity_relative']
        reports.append(dict(family=family,points=len(x),max_scalar_relative_error=relative,max_scalar_absolute_error=absolute,
                            max_order_doubling_relative_error=order_error,scalar_integration_warning_count=len(caught),moments=checks,
                            lift_identity_error=identity,vector_seconds=elapsed))
        print(json.dumps(reports[-1]),flush=True)
    report=dict(schema='blackhole-channel-fast-bridge-validation/1',status='pass',checks=reports,
                input_sha256={p.name:sha(p) for p in [BASE/'bridge-design.json',BASE/'fast_measure.py',BASE/'projected_measure.py',Path(__file__)]},
                scope=design['limits'])
    (BASE/'fast-measure-validation.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    try:main()
    except Exception as error:
        p=BASE/'fast-measure-failures.json';history=json.loads(p.read_text()) if p.exists() else []
        history.append(dict(tool_sha256=sha(Path(__file__)),module_sha256=sha(BASE/'fast_measure.py'),error=repr(error)))
        p.write_text(json.dumps(history,indent=2)+'\n',encoding='utf-8')
        raise
