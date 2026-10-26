"""Analytic moments and quadrature checks of the declared4D/7D interface."""
from pathlib import Path
import hashlib
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from projected_measure import effective_spin_prior,uniform_reduced_prior,lift_channel_density

BASE=Path(__file__).resolve().parent


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    # Gates fixed before these new prior integrations:2e-5 normalization,
    #5e-6 second-moment error, and1e-12 algebraic change-of-measure identity.
    nodes,weights=leggauss(512);amax=.99;x=amax*nodes
    history_path=BASE/'measure-history.json'
    history=json.loads(history_path.read_text()) if history_path.exists() else []
    # Preserve the first terminal run's observed failure. Refine the integral
    # by splitting the known |x| cusp at zero; do not relax the2e-5 gate.
    coarse_pdf=np.array([effective_spin_prior(value,.1,amax) for value in x])
    coarse_norm=float(amax*np.dot(weights,coarse_pdf))
    initial=dict(family='uniform_magnitude',q=.1,nodes=512,rule='one Gauss-Legendre interval across[-A,A]',normalization=coarse_norm,normalization_gate=2e-5,status='pass' if abs(coarse_norm-1)<2e-5 else 'failed')
    if initial not in history:history.append(initial)
    rows=[]
    for family,variance_denominator in [('uniform_magnitude',9),('beta2_5_magnitude',28)]:
        for q in [.1,.5,1.]:
            positive=amax*(nodes+1)/2
            pp=np.array([effective_spin_prior(v,q,amax,family) for v in positive])
            pn=np.array([effective_spin_prior(-v,q,amax,family) for v in positive])
            norm=float(amax/2*np.dot(weights,pp+pn));mean=float(amax/2*np.dot(weights,positive*(pp-pn)))
            second=float(amax/2*np.dot(weights,positive**2*(pp+pn)))
            expected=amax**2*(1+q*q)/(variance_denominator*(1+q)**2)
            assert abs(norm-1)<2e-5,(family,q,norm)
            assert abs(mean)<1e-10
            assert abs(second-expected)<5e-6,(family,q,second,expected)
            rows.append(dict(family=family,q=q,normalization=norm,mean=mean,second_moment=second,analytic_second_moment=expected))
            print(f'{family}, q={q}: normalized scalar spin prior checked',flush=True)
    checks=[]
    for q in [.1,.5,1.]:
        prior7=.37;density4=.012;a1=.3;a2=.5;c1=.2;c2=-.4
        x=(a1*c1+q*a2*c2)/(1+q)
        reduced=uniform_reduced_prior(prior7,q,x,amax)
        lifted=lift_channel_density(density4,q,a1,a2,c1,c2,amax)
        error=abs(lifted/prior7/(density4/reduced)-1)
        assert error<1e-12
        checks.append(dict(q=q,full_to_reduced_weight_identity_relative_error=error,
                           naive_density4_over_prior7_ratio_to_correct_weight=float(reduced/prior7)))
    report=dict(schema='blackhole-physical-channel-measure-bridge/1',status='pass',
                tool_sha256=sha(Path(__file__)),module_sha256=sha(BASE/'projected_measure.py'),checks=rows,
                algebraic_checks=checks,
                quadrature='512-point Gaussian rules separately on[-A,0] and[0,A]; original unsplit512-node failure retained in measure-history.json. Original tolerances unchanged.',
                formula_reduced='pi4(Mc,q,x,z)=pi7(m1,q,a1,a2,c1,c2,z)*4*A^2*(dm1/dMc)*f_U(x|q), when pi7=pi_mass,z*Q_U on its support.',
                formula_lift='p7=p4(Mc,q,x,z)*(dMc/dm1)*Q_spin(a1,a2,c1,c2|q)/f_Q(x|q).',
                limits='Exact declared independent-magnitude/isotropic-angle conditionals, not proof that they are physical formation-channel spin conditionals or match every historical PE/injection prior. No GW data fit or selection calibration.')
    history.append(dict(rule='split at zero,512 nodes per half',status='pass',maximum_normalization_error=max(abs(r['normalization']-1) for r in rows),tool_sha256=sha(Path(__file__))))
    history_path.write_text(json.dumps(history,indent=2)+'\n',encoding='utf-8')
    (BASE/'measure-validation.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='checks'},indent=2))


if __name__=='__main__':main()
