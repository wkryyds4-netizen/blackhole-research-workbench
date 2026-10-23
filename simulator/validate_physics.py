"""Physical invariants and signal processing checks; writes numerical evidence."""
from pathlib import Path
import json
import numpy as np
from physics import waveform,validate,DEFAULT,fit_event,fd_waveform,FREQ

def main():
    p={**DEFAULT,'inclination':0.}
    a=waveform(p);b=waveform({**p,'distance':880.})
    x=np.asarray(a['hp']);y=np.asarray(b['hp'])
    distance_error=float(np.max(np.abs(x-2*y))/np.max(np.abs(x)))
    assert distance_error<1e-6
    edge=waveform({**p,'inclination':90.})
    inclination_error=float(np.max(np.abs(x-2*np.asarray(edge['hp'])))/np.max(np.abs(x)))
    assert inclination_error<1e-6 and np.max(np.abs(edge['hc']))<1e-12
    t=np.array(a['t']);freq=np.array([np.nan if x is None else x for x in a['frequency']])
    early=float(np.nanmedian(freq[(t<-.3)&(t>-.5)]));late=float(np.nanmedian(freq[(t<-.02)&(t>-.06)]))
    assert late>early>15, (early,late)
    assert abs(t[np.argmax(a['envelope'])])<1/a['sample_rate']
    spun=waveform({**p,'s1':.6,'s2':.6})
    assert not np.allclose(np.interp(t,spun['t'],spun['hp']),a['hp'],atol=1e-4)
    # Scale-free vacuum BBH mass relation tested in frequency domain away from tapers.
    h1,_=fd_waveform(validate(p));h2,_=fd_waveform(validate({**p,'m1':72.,'m2':58.}))
    idx=np.where((FREQ>=30)&(FREQ<=150))[0]
    # H(2M,f)=4 H(M,2f) up to a linear time/phase offset (f_ref convention).
    ratio=np.abs(h2[idx])/(4*np.maximum(np.abs(h1[2*idx]),1e-40))
    mass_error=float(np.max(np.abs(ratio-1)))
    assert mass_error<1e-5,mass_error
    alt=waveform({**p,'model':'IMRPhenomXAS'})
    assert np.all(np.isfinite(alt['hp']))
    for bad in [{'m1':10},{'m2':90},{'s1':2},{'distance':float('nan')},{'model':'unknown'}]:
        try:validate({**DEFAULT,**bad})
        except ValueError:pass
        else:raise AssertionError(f'Accepted {bad}')
    fits={}
    for det in ['H1','L1']:
        fit=fit_event({**DEFAULT,'m1':39.,'m2':32.},det)
        assert 0<=fit['correlation']<=1 and len(fit['residual'])==len(fit['observed'])
        assert np.allclose(np.array(fit['observed'])-np.array(fit['template']),fit['residual'],atol=2e-6)
        fits[det]={k:fit[k] for k in ['correlation','time_shift_ms']}
    result={'passed':True,'distance_relative_error':distance_error,'inclination_relative_error':inclination_error,'mass_FD_amplitude_scaling_error':mass_error,'early_frequency_hz':early,'late_frequency_hz':late,'version':a['meta'],'event_shape_fits':fits,'limitations':'No independent LALSuite or NR accuracy comparison in this environment. These are invariants, direction and pipeline checks.'}
    path=Path(__file__).parent/'validation.json';path.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
