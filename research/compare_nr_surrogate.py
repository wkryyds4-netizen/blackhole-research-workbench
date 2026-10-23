"""Fixed-protocol SXS0389 / NRSur7dq4 mode experiment; no held-out claim."""
from pathlib import Path
import sys, json, hashlib
import numpy as np
import h5py
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar
from scipy.spatial.transform import Rotation
from scipy.signal.windows import tukey
from waveform_reference import interpolate_horizon_reference, rotate_complete_modes
from waveform_metrics import compare_waveforms, common_grid, MTSUN_S

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'stage3-results'
OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT/'adapter-build/zig-runtime/gwsurrogate-1.2.0'))
import gwsurrogate as gw

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read_nr(order):
    p=ROOT/f'wave-source-results/sxs-0389-lev3-N{order}.npz'
    with np.load(p) as z:
        return z['t'], {(l,m):z[f'h_{l}_{m}'] for l in range(2,5) for m in range(-l,l+1)}

def peak(t, modes):
    power=np.sum([np.abs(h)**2 for h in modes.values()],axis=0)
    i=int(np.argmax(power)); sp=CubicSpline(t,power)
    return float(minimize_scalar(lambda x:-sp(x),bounds=(t[i-1],t[i+1]),method='bounded',options={'xatol':1e-10}).x)

def main():
    meta=json.loads((ROOT/'sxs-data/0389/Lev3/metadata.json').read_text())
    source=json.loads((ROOT/'wave-source-results/sxs-0389-lev3-N2.json').read_text())
    mass=meta['reference_mass1']+meta['reference_mass2']
    assert abs(source['normalization_ChMass_code_units']-mass)<1e-12
    tn,hn=read_nr(2); tp=peak(tn,hn)
    reference_time=(tp-4300)*mass
    assert reference_time>meta['relaxation_time']
    with h5py.File(ROOT/'sxs-data/0389/Lev3/Horizons.h5') as f:
        params=interpolate_horizon_reference(meta,f,reference_time)
    quat=Rotation.from_matrix(np.array(params['coorbital_to_source_matrix']).T).as_quat()[[3,0,1,2]]
    nr=rotate_complete_modes(hn,quat)[(2,2)]
    times=tn-tp
    model=ROOT/'surrogate-data/NRSur7dq4.h5'
    assert sha(model)=='3fc2f58f6ed7a4fff58856b6c01a7dc9155954d2b9fa5d58126a7807826d45d5'
    sur=gw.LoadSurrogate(str(model))
    ts,hs,_=sur(params['q'],params['chiA0'],params['chiB0'],f_low=0.,f_ref=0.,dt=.5,units='dimensionless',ellMax=4)
    np.savez_compressed(OUT/'matched-surrogate.npz',t=ts,**{f'h_{l}_{m}':v for (l,m),v in hs.items()})
    # The common interval is fixed before evaluation; all variants use N2's peak.
    interval=[-4290.,90.]
    assert ts[0]<=interval[0] and ts[-1]>=interval[1]
    t,a,b=common_grid(times,nr,ts,hs[(2,2)],1.,interval=interval)
    # Local phase derivative near the finite record's beginning; not a support proof.
    start_f=float(np.median(np.abs(np.gradient(np.unwrap(np.angle(a[:101])),t[:101])))/(2*np.pi))
    rows=[]; full=[]
    def measure(a,b,dt,M,psd,fmin,taper,comparison,label):
        result=compare_waveforms(a,b,dt*M*MTSUN_S,fmin_hz=fmin,fmax_hz=512.,psd=psd,taper_fraction=taper,max_time_shift_s=20*M*MTSUN_S)
        row={'id':f'{comparison}-{M}-{psd}-{fmin}-{taper}-{dt}','label':label,'comparison':comparison,'mass_msun':M,'psd':psd,'fmin_hz':fmin,'fmax_hz':512.,'taper_fraction':taper,'dt_M':dt,**{k:result[k] for k in ['match','mismatch','phase_shift_rad','at_time_boundary']},'time_shift_M':result['time_shift_s']/(M*MTSUN_S),'start_frequency_hz':start_f/(M*MTSUN_S),'start_band_supported':bool(start_f/(M*MTSUN_S)<fmin),'window_removed_power':result['config']['unweighted_power_removed_by_taper']}
        row['id'] += '-'+label.replace(' ','').replace('/','-').replace(',','-')
        row['peak_window_weight'] = float(np.interp(0, np.linspace(*interval,len(a)), tukey(len(a),alpha=taper)))
        rows.append(row);full.append({'id':row['id'],**result});return row
    for M in [40,65,100,150]:
        for psd in ['white','aligo_historical']:
            for fmin in [20,30]:
                for taper in [.02,.05,.1]:
                    measure(a,b,1.,M,psd,fmin,taper,'surrogate','NR N2 / NRSur7dq4')
    baseline=next(r for r in rows if r['mass_msun']==65 and r['psd']=='aligo_historical' and r['fmin_hz']==20 and r['taper_fraction']==.02)
    for order in [3,4]:
        tx,hx=read_nr(order); hx=rotate_complete_modes(hx,quat)[(2,2)]
        tx,ax,bx=common_grid(times,nr,tx-tp,hx,1.,interval=interval)
        for M in [40,65,100,150]:
            for psd in ['white','aligo_historical']:
                for fmin in [20,30]:
                    for taper in [.02,.05,.1]:
                        measure(ax,bx,1.,M,psd,fmin,taper,'extrapolation',f'NR N2 / N{order}')
    # Half-step comparison changes both interpolation/sampling and FFT quadrature.
    tx,ax,bx=common_grid(times,nr,ts,hs[(2,2)],.5,interval=interval)
    fine=measure(ax,bx,.5,65,'aligo_historical',20,.02,'sampling','NR / surrogate, dt=0.5M')
    # Display only points supported after the optimized shift, without amplitude fitting.
    display_t=np.arange(-300.,80.01,.5)
    shift=baseline['time_shift_M']; phase=baseline['phase_shift_rad']
    assert (display_t-shift).min()>=ts.min() and (display_t-shift).max()<=ts.max()
    display_nr=np.interp(display_t,times,nr)
    display_sur=np.interp(display_t-shift,ts,hs[(2,2)])*np.exp(1j*phase)
    provenance={'nr':'SXS:BBH:0389 Lev3; N2/N3/N4 CoM-corrected','nr_url':'https://doi.org/10.5281/zenodo.3274661','nr_raw_sha256':source['source_sha256'],'model':'NRSur7dq4 / gwsurrogate 1.2.0','model_url':'https://doi.org/10.5281/zenodo.3348115','model_sha256':sha(model),'training_membership':'unverified','independent_validation':False,'parameters':params}
    protocol={'nr_peak_time_M':tp,'peak_method':'cubic spline of sum |h_lm|^2, all ell=2..4; local maximum refinement','normalization':'NR already t/M, rh/M with ChMass=Mref; no second division','reference_horizon_time_source':reference_time,'clock':'t_horizon=(t_peak_NR/M - 4300)*Mref; gauge-dependent numerical clock identification','frame':'all 21 NR modes passively rotated into coordinate orbital triad at reference epoch; not a calibrated radiation frame','surrogate_reference':'f_ref=0: initial model epoch -4300M','interval_M':interval,'sampling':'linear interpolation on common support; no extrapolation','time_search_M':[-20,20],'signal_kind':'complex h22, two-sided PSD(|f|), global complex phase/time optimization','mass':'detector-frame total mass; no redshift or distance inference','psd':'white or historical LALAdvLIGOPsd analytic shape; not measured detector PSD','taper':'symmetric Tukey; alpha scanned .02/.05/.1; long taper can suppress merger','display':'baseline65Msun/20..512Hz/historicalPSD/alpha.02; unwindowed complex modes; optimized time and phase only, no amplitude scaling','start_frequency_diagnostic':'median phase derivative over first100M; f_start<fmin necessary heuristic, not complete spectral support guarantee'}
    limits=['外部归档对照，训练成员关系尚未核实，不能声称训练集外验证。','复 h22 模式失配不是探测器网络失配，也不能推出形成渠道概率。','坐标轨道框架近似与辐射框架误差尚未隔离；此处为近无自旋案例。','N2/N3/N4 是外推阶敏感性，不是分辨率收敛结果或统计误差条。','有限长度、渐消窗和频段会改变失配；起始频率诊断不保证完整低频支持。']
    data={'schema':'blackhole-waveform-study/1','baseline_id':baseline['id'],'rows':rows,'waveform':{'time_M':display_t.tolist(),'nr_re':display_nr.real.tolist(),'surrogate_re':display_sur.real.tolist(),'residual_abs':np.abs(display_nr-display_sur).tolist()},'provenance':provenance,'protocol':protocol,'limits':limits,'sampling_check':{'baseline_mismatch':baseline['mismatch'],'half_step_mismatch':fine['mismatch'],'absolute_change':abs(fine['mismatch']-baseline['mismatch'])}}
    for p in [OUT/'waveform-comparison.json',ROOT.parent/'simulator/public/datasets/waveform-comparison.json']:
        p.write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
    (OUT/'metrics-full.json').write_text(json.dumps(full,indent=2),encoding='utf-8')
    print(json.dumps({'baseline':baseline,'sampling_check':data['sampling_check'],'rows':len(rows),'params':params},indent=2),flush=True)

if __name__=='__main__':main()
