"""Aligned-spin IMR waveforms. Rendering never feeds back into the calculation."""
from pathlib import Path
from functools import lru_cache
import importlib.metadata
import threading
import numpy as np
from scipy import signal
import jax
import jax.numpy as jnp
import ripplegw
from prepare_event import whiten_filter, FS as EVENT_FS

jax.config.update('jax_enable_x64', True)
ROOT=Path(__file__).resolve().parent
FS=4096
N=FS*32
MTSUN=4.9254909476412675e-6
FREQ=np.fft.rfftfreq(N,1/FS)
LOCK=threading.Lock()
DEFAULT={'m1':36.,'m2':29.,'s1':0.,'s2':0.,'distance':440.,'inclination':35.,'model':'IMRPhenomD'}
MODELS={}

def validate(p):
    p={**DEFAULT,**p}
    ranges={'m1':(20,100),'m2':(10,100),'s1':(-.8,.8),'s2':(-.8,.8),'distance':(100,2000),'inclination':(0,180)}
    for k,(lo,hi) in ranges.items():
        p[k]=float(p[k])
        if not np.isfinite(p[k]) or not lo<=p[k]<=hi:raise ValueError(f'{k} must be in [{lo}, {hi}]')
    if p['m1']<p['m2'] or p['m1']/p['m2']>4:raise ValueError('Require m1 >= m2 and m1/m2 <= 4')
    if p['model'] not in ['IMRPhenomD','IMRPhenomXAS']:raise ValueError('Unknown waveform model')
    return {k:p[k] for k in DEFAULT}

def fd_waveform(p):
    name=p['model']
    if name not in MODELS:
        # Pin the published 0.2.1 interface, not the development documentation API.
        MODELS[name]=jax.jit(getattr(ripplegw,name)(f_ref=20.))
    total=p['m1']+p['m2'];eta=p['m1']*p['m2']/total**2
    pars={'M_c':total*eta**.6,'eta':eta,'s1_z':p['s1'],'s2_z':p['s2'],'d_L':p['distance'],'phase_c':0.,'iota':0.}
    # Ripple uses f[1]-f[0] to set its high-frequency cutoff: never clamp the
    # first bins to a repeated frequency. Evaluate the uniform positive grid.
    start=int(np.searchsorted(FREQ,20.))
    out=MODELS[name](jnp.asarray(FREQ[start:]),{k:jnp.asarray(v,dtype=jnp.float64) for k,v in pars.items()})
    # Smooth analysis-band edges. They are documented processing, not source physics.
    low=np.sin(np.clip((FREQ-20)/5,0,1)*np.pi/2)**2
    high=np.cos(np.clip((FREQ-900)/100,0,1)*np.pi/2)**2
    hp=np.zeros(len(FREQ),complex);hc=np.zeros(len(FREQ),complex)
    hp[start:]=np.asarray(out['p']);hc[start:]=np.asarray(out['c'])
    return hp*low*high,hc*low*high

def compute(p):
    p=validate(p)
    with LOCK:
        hf,cf=fd_waveform(p)
    hp=np.fft.fftshift(np.fft.irfft(hf,n=N))*FS
    hc=np.fft.fftshift(np.fft.irfft(cf,n=N))*FS
    env=np.hypot(hp,hc);peak=int(np.argmax(env))
    t=(np.arange(N)-peak)/FS
    phase=np.unwrap(np.angle(hp+1j*hc))
    freq=np.abs(np.gradient(phase)*FS/(2*np.pi))
    mass=p['m1']+p['m2'];eta=p['m1']*p['m2']/mass**2;mc=mass*eta**.6
    tau=5/256*(mc*MTSUN)**(-5/3)*(np.pi*20)**(-8/3)
    start=-min(12,max(.55,1.5*tau))
    take=(t>=start)&(t<=.12)
    ci=np.cos(np.deg2rad(p['inclination']))
    hplus=hp*(1+ci**2)/2;hcross=hc*ci
    # Approximate separation inferred only in inspiral, in units GM/c^2.
    valid=(t<-.012)&(env>env.max()*.025)&(freq>18)&(freq<900)
    freq=np.where(valid,freq,np.nan)
    sep=np.where(valid,(1/(np.pi*mass*MTSUN*np.nan_to_num(freq,nan=1)))**(2/3),np.nan)
    def arr(a,scale=1):return [None if not np.isfinite(v) else float(f'{v*scale:.8g}') for v in a[take]]
    result={'parameters':p,'t':arr(t),'hp':arr(hplus,1e21),'hc':arr(hcross,1e21),'envelope':arr(env,1e21),'phase':arr(phase-phase[peak]),'frequency':arr(freq),'separation':arr(sep),'sample_rate':FS,'meta':{'backend':'rippleGW','version':importlib.metadata.version('rippleGW'),'mass_frame':'detector','strain_unit':'1e-21','f_ref_hz':20,'taper_hz':[20,25,900,1000],'peak_alignment':'face-on quadrature-envelope maximum','chirp_mass':mc,'eta':eta,'q':p['m2']/p['m1'],'chi_eff':(p['m1']*p['s1']+p['m2']*p['s2'])/mass,'peak_strain':float(np.max(np.hypot(hplus,hcross))),'duration':float(t[take][-1]-t[take][0])}}
    if not np.all(np.isfinite(hplus)) or not np.all(np.isfinite(hcross)):raise RuntimeError('Nonfinite waveform')
    return result

@lru_cache(maxsize=16)
def cached_wave(key):return compute(dict(key))

def waveform(p):return cached_wave(tuple(validate(p).items()))

def fit_event(p,detector):
    """Independent time/phase/amplitude display fit, NOT coherent inference."""
    if detector not in ['H1','L1']:raise ValueError('Unknown detector')
    w=waveform(p);data=np.load(ROOT/'data'/f'{detector}.npz',allow_pickle=False)
    t=data['time'];y=data['white'];psd=data['psd']
    # Face-on independent quadratures avoid rank loss at an edge-on inclination.
    base=waveform({**p,'inclination':0.})
    a=np.interp(t,base['t'],np.asarray(base['hp'])*1e-21,left=0,right=0)
    b=np.interp(t,base['t'],np.asarray(base['hc'])*1e-21,left=0,right=0)
    a=whiten_filter(a,psd);b=whiten_filter(b,psd)
    fit=(t>=-.2)&(t<=.08);idx=np.where(fit)[0];target=y[fit]
    best=None
    for shift in range(-164,165,2):
        design=np.column_stack([a[idx-shift],b[idx-shift]])
        coeff=np.linalg.lstsq(design,target,rcond=None)[0]
        pred=design@coeff;loss=float(np.sum((target-pred)**2))
        if best is None or loss<best[0]:best=(loss,shift,coeff,pred)
    loss,shift,coeff,pred=best
    full=np.roll(a,shift)*coeff[0]+np.roll(b,shift)*coeff[1]
    corr=float(np.dot(target,pred)/np.sqrt(np.dot(target,target)*np.dot(pred,pred)))
    take=(t>=-.5)&(t<=.25)
    return {'detector':detector,'t':np.round(t[take],7).tolist(),'observed':np.round(y[take],6).tolist(),'template':np.round(full[take],6).tolist(),'residual':np.round((y-full)[take],6).tolist(),'correlation':corr,'time_shift_ms':shift/EVENT_FS*1000,'fit_window':[-.2,.08],'fit_type':'Independent amplitude/phase/time least-squares display fit; not posterior inference or coherent detector projection','parameters':w['parameters']}
