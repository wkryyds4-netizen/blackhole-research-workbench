"""Deterministic precession-averaged transport with an explicit native reference."""
import numpy as np
import precession
from scipy.integrate import solve_ivp


def to_infinity(theta1,theta2,phi12,q,a1,a2,total_detector_mass,reference_hz,
                rtol=1e-9,atol=1e-11,roundtrip=False,solver='upstream',initial_separation=None):
    values=np.array([theta1,theta2,phi12,q,a1,a2,total_detector_mass,reference_hz])
    if not np.isfinite(values).all():raise ValueError('Every transport input, including native reference, must be finite.')
    if not (0<q<1 and 0<a1<=1 and 0<a2<=1 and total_detector_mass>0 and reference_hz>0):
        raise ValueError('This adapter requires unequal positive masses/spins and a positive native reference; no default is supplied.')
    if not (0<=theta1<=np.pi and 0<=theta2<=np.pi):raise ValueError('Tilt outside [0,pi].')
    angles=dict(theta1=theta1,theta2=theta2,deltaphi=phi12,q=q,chi1=a1,chi2=a2)
    r0=(float(precession.gwfrequency_to_pnseparation(**angles,fGW=reference_hz,M_msun=total_detector_mass).item())
        if initial_separation is None else float(initial_separation))
    if not np.isfinite(r0) or r0<=0:raise ValueError('Nonphysical PN separation.')
    # The default angular outputs resample precession phases at finite separation.
    # Request only deterministic conserved quantities and derive infinity tilts.
    if solver=='upstream':
        # 2.1.2 checks misspelled 'rol'/'aol' keys and overwrites rtol/atol.
        # Do not pretend a passed tolerance changes this branch's integration.
        answer=precession.inspiral_precav(**angles,r=[r0,np.inf],requested_outputs=['kappa','chieff'])
        kin,kout=map(float,answer['kappa'][0]);ce=float(answer['chieff'].item())
    elif solver=='DOP853':
        _,k0,chieff=precession.angles_to_conserved(**angles,r=r0)
        kin=float(k0.item());ce=float(chieff.item());u0=float(precession.eval_u(r0,q).item())
        def rhs(u,k):return [float(np.asarray(precession.rhs_precav(k,u,ce,q,a1,a2)).item())]
        integration=solve_ivp(rhs,(u0,0.),[kin],method='DOP853',rtol=rtol,atol=atol)
        if not integration.success:raise RuntimeError(integration.message)
        kout=float(integration.y[0,-1])
    else:raise ValueError('Unknown solver.')
    dc=precession.eval_deltachiinf(kout,ce,q,a1,a2)
    c1=float(precession.eval_costheta1(dc,ce,q,a1).item());c2=float(precession.eval_costheta2(dc,ce,q,a2).item())
    if not np.isfinite([c1,c2]).all() or max(abs(c1),abs(c2))>1+1e-7:raise ValueError('Unphysical transported cosine.')
    result=dict(pn_separation=r0,cos_tilt_1_infinity=c1,cos_tilt_2_infinity=c2,
                chieff=ce,chieff_transport_error=float(abs((a1*c1+q*a2*c2)/(1+q)-ce)),
                kappa_initial=kin,kappa_infinity=kout,pn_separation_below_6=r0<6,solver=solver,
                initial_separation_override=initial_separation is not None)
    if roundtrip:
        if solver=='upstream':
            back=precession.inspiral_precav(kappa=kout,chieff=ce,q=q,chi1=a1,chi2=a2,r=[np.inf,r0],requested_outputs=['kappa'])
            back_value=float(back['kappa'][0,-1])
        else:
            back=solve_ivp(rhs,(0.,u0),[kout],method='DOP853',rtol=rtol,atol=atol)
            if not back.success:raise RuntimeError(back.message)
            back_value=float(back.y[0,-1])
        result['roundtrip_kappa_abs_error']=float(abs(back_value-kin))
    return result
