"""Reproducible NRSur7dq4 evaluation and independent runtime checks."""
from pathlib import Path
import sys,json,hashlib,time,subprocess
import numpy as np
import h5py
from scipy.spatial.transform import Rotation
root=Path(__file__).resolve().parent
runtime=root/'adapter-build/zig-runtime/gwsurrogate-1.2.0'
sys.path.insert(0,str(runtime))
import gwsurrogate as gw
from gwsurrogate.precessing_utils import _utils
from gwsurrogate.new.precessing_surrogate import transformTimeDependentVector
out=root/'surrogate-results';out.mkdir(exist_ok=True)
model=root/'surrogate-data/NRSur7dq4.h5'
modelhash=hashlib.sha256(model.read_bytes()).hexdigest()
assert modelhash=='3fc2f58f6ed7a4fff58856b6c01a7dc9155954d2b9fa5d58126a7807826d45d5'

# Compare every stored fit, independently evaluating the polynomial in NumPy.
rng=np.random.default_rng(20261023)
errors=[]; checked=0
with h5py.File(model) as f:
    pairs=[]
    def find(name,d):
        if isinstance(d,h5py.Dataset) and ('bfOrders' in name):
            coef=name.replace('bfOrders','coefs')
            if coef in f: pairs.append((d[()],f[coef][()]))
    f.visititems(find)
    for orders,coefs in pairs:
        assert orders.dtype==np.dtype('int64')
        x=rng.uniform(-.6,.6,7);x[0]=np.log(rng.uniform(1,4))
        xp=x.copy();xp[0]=-.9857019407834238+1.4298059216576398*xp[0]
        ref=np.sum(coefs*np.prod(xp[None,:]**orders,axis=1))
        actual=_utils.eval_fit(orders,coefs,x,-.9857019407834238,1.4298059216576398,3,2)
        errors.append(abs(actual-ref));checked+=1
assert checked>1000 and max(errors)<1e-12,(checked,max(errors))

params={'q':2.,'chiA0':[.4,.1,.2],'chiB0':[-.1,.2,-.1],'f_low':0.,'f_ref':0.,'dt':2.,'units':'dimensionless','ellMax':4}
sur=gw.LoadSurrogate(str(model))
started=time.perf_counter()
t,h,d=sur(**params,precessing_opts={'return_dynamics':True})
elapsed=time.perf_counter()-started
print('evaluated',len(t),list(d),elapsed,flush=True)
arrays=[t,*h.values(),*d.values()]
assert all(np.all(np.isfinite(a)) for a in arrays)
assert np.all(np.diff(t)>0)
q=d['q_copr'].T
spin_error=max(np.max(np.abs(np.linalg.norm(d[k],axis=1)-np.linalg.norm(params[k+'0']))) for k in ['chiA','chiB'])
quat_error=float(np.max(np.abs(np.linalg.norm(q,axis=1)-1)))
assert quat_error<1e-12 and spin_error<1e-12
rot=Rotation.from_quat(q[:,[1,2,3,0]])
rotation_error=max(float(np.max(np.abs(rot.apply(d[k+'_copr'])-d[k]))) for k in ['chiA','chiB'])
assert rotation_error<1e-12
z=np.tile([0.,0.,1.],(len(t),1))
axis=rot.apply(z)
initial_error=max(float(np.max(np.abs(d[k][0]-params[k+'0']))) for k in ['chiA','chiB'])
assert initial_error<1e-7
# The same common times evaluated directly must reproduce the dt-generated output.
tx,hx,dx=sur(params['q'],params['chiA0'],params['chiB0'],f_low=0.,f_ref=0.,times=t[::7],ellMax=4,precessing_opts={'return_dynamics':True})
repeat_error=max(float(np.max(np.abs(hx[k]-v[::7]))) for k,v in h.items())
assert repeat_error<1e-12
np.savez_compressed(out/'nrsur7dq4-full.npz',t=t,**{'h_%d_%d'%k:v for k,v in h.items()},**d)
metadata={'params':params,'model':'NRSur7dq4','version':gw.__version__,'model_sha256':modelhash,'source_archive_sha256':'316a640576d32aaa45aaf9a6d876b240af356c89131e6af8257039b3818baae8','source':'https://doi.org/10.5281/zenodo.3348115','model_license':'CC-BY-4.0','code_license':'MIT','training_domain':{'q':[1,4],'spin_magnitude_max':.8},'time_unit':'GM/c^3','waveform_unit':'r*h/M','waveform_kind':'inertial-frame spherical harmonic mode h_22, not detector strain','reference_frequency':{'value':0,'unit':'cycles/M','meaning':'special value: initial model epoch; not a physical zero-frequency orbit'},'reference_time_M':float(t[0]),'time_convention':'model geometric time; t=0 at model total multimode waveform amplitude peak; not necessarily the h22 peak and not horizon common-formation time','frame':'inertial: z along L and x from lighter to heavier BH at reference epoch; y completes right-handed triad','quaternion':'[w,x,y,z]; active coprecessing-to-inertial rotation v_I = q v_C q^-1','spin_unit':'dimensionless chi = S/m_i^2','dynamics_valid_until_M':-100,'dynamics_cutoff_reason':'At t/M=-100 the NR training spins are extended by PN equations; later spins are unphysical fit parameters (Varma et al. 2019, Sec. III C). Direction display is conservatively cut at this boundary too.','dynamics_cutoff_source':'https://arxiv.org/abs/1905.09300','position_available':False,'separation_available':False,'physical_mass_msun':None,'distance_mpc':None,'notes':['No mapping to SXS coordinate clocks is assumed.','This demonstration is inside the training domain, not an accuracy assessment against independent NR waveforms.']}
dataset={'schema':'blackhole-dynamics/1','title':'NRSur7dq4 进动演示 q=2','metadata':metadata,'dynamics_valid_until_M':-100,'t':t.tolist(),'q_copr':q.tolist(),'chiA':d['chiA'].tolist(),'chiB':d['chiB'].tolist(),'orbphase':d['orbphase'].tolist(),'h22_real':h[(2,2)].real.tolist(),'h22_imag':h[(2,2)].imag.tolist(),'coprecessing_z_axis':axis.tolist()}
(out/'dynamics.json').write_text(json.dumps(dataset,ensure_ascii=False,separators=(',',':')),encoding='utf-8')
report={'status':'evaluated','samples':len(t),'time_range_M':[float(t[0]),float(t[-1])],'mode_count':len(h),'seconds':elapsed,'fit_polynomials_checked':checked,'max_fit_abs_error':max(errors),'max_quaternion_norm_error':quat_error,'max_spin_norm_error':float(spin_error),'max_rotation_vs_scipy_error':rotation_error,'initial_spin_error':initial_error,'repeat_grid_max_mode_error':repeat_error,'output_sha256':hashlib.sha256((out/'dynamics.json').read_bytes()).hexdigest(),'npz_sha256':hashlib.sha256((out/'nrsur7dq4-full.npz').read_bytes()).hexdigest(),'metadata':metadata}
(out/'validation.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k!='metadata'},indent=2))

