"""Verify restored channel density maps; no GW channel-fraction inference."""
from pathlib import Path
import hashlib
import json
from importlib.metadata import version
import h5py
import numpy as np
from scipy.stats import ks_2samp
import torch
from channel_flow import ChannelFlow,CHANNELS,PRODUCTS

BASE=Path(__file__).resolve().parent


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    torch.set_num_threads(4)
    design=json.loads((BASE/'runtime-design.json').read_text());gates=design['gates']
    report=dict(schema='blackhole-physical-channel-flow-runtime/1',status='pass',checks={},samples={},
                runtime={k:version(k) for k in ['numpy','scipy','torch','glasflow']},
                input_sha256={p.name:sha(p) for p in [BASE/'runtime-design.json',BASE/'channel_flow.py',Path(__file__),BASE/'product-inspection.json']},
                meaning='Numerical restoration of authors physical-channel emulators, without fitting any observed events or claiming calibrated learned tails.')
    saved={};count=design['sampling']['samples_per_channel'];n_test=design['sampling']['latent_test_points']
    reference=np.array([[8,.25,-.3,.1],[20,.6,.01,.5],[50,.9,.5,2],[90,.1,-.8,7.]])
    for i,channel in enumerate(CHANNELS):
        f=ChannelFlow(channel);ctx=[0.,2.] if channel=='CE' else [0.]
        mapped=f.to_network(reference);restored=f.to_physical(mapped)
        map_error=float(np.max(abs(restored/reference-1)));assert map_error<gates['coordinate_roundtrip_relative'],(channel,map_error)
        numerical=[]
        for point in reference:
            jac=np.zeros((4,4))
            for col in range(4):
                h=1e-5*max(1.,abs(point[col]));upper=point.copy();lower=point.copy();upper[col]+=h;lower[col]-=h
                jac[:,col]=(f.to_network(upper[None,:])[0]-f.to_network(lower[None,:])[0])/(2*h)
            numerical.append(np.linalg.slogdet(jac)[1])
        coordinate_error=float(np.max(abs(np.array(numerical)-f.log_coordinate_jacobian(reference))))
        assert coordinate_error<gates['coordinate_log_jacobian_absolute'],(channel,coordinate_error)
        rng=np.random.default_rng(design['sampling']['seed']+i)
        latent=rng.normal(size=(count,4));conditional=f.context(ctx,n_test)
        with torch.no_grad():
            network_points,inverse_det=f.network.inverse(torch.as_tensor(latent[:n_test]),conditional=conditional)
            returned,forward_det=f.network.forward(network_points,conditional=conditional)
        inverse_error=float(np.max(abs(returned.numpy()-latent[:n_test])))
        det_error=float(torch.max(abs(inverse_det+forward_det)))
        assert inverse_error<gates['network_latent_roundtrip_absolute'],(channel,inverse_error)
        assert det_error<gates['network_logdet_cancellation_absolute'],(channel,det_error)
        theta=f.to_physical(network_points.numpy())
        normal=-.5*np.sum(latent[:n_test]**2,axis=1)-2*np.log(2*np.pi)
        expected_original=normal+forward_det.numpy()+f.log_coordinate_jacobian(theta)
        actual=f.log_density(theta,ctx)
        original_density_error=float(np.max(abs(actual-expected_original)))
        # Physical coordinates are stored float64 numbers. Near q=1, inverse
        # logistic rounding means these are not exactly the original inverse
        # network points. Verify density at the represented point, while keeping
        # the original sampled-latent comparison and its failed gate visible.
        with torch.no_grad():
            represented_z,represented_det=f.network.forward(torch.as_tensor(f.to_network(theta)),conditional=conditional)
        expected=-.5*np.sum(represented_z.numpy()**2,axis=1)-2*np.log(2*np.pi)+represented_det.numpy()+f.log_coordinate_jacobian(theta)
        density_error=float(np.max(abs(actual-expected)))
        assert density_error<gates['density_change_of_variables_absolute'],(channel,density_error)
        # Autodiff of all4 composite outputs, independently taking the
        # determinant instead of using the layer's reported log determinant.
        determinant_errors=[]
        scale=torch.as_tensor(f.scale);upper=torch.as_tensor(f.upper)
        for point in theta[:design['sampling']['autodiff_points']]:
            def composite(raw):
                components=[]
                for col in range(4):
                    components.append(torch.atanh(raw[col]) if col==2 else torch.logit(raw[col]/upper[col])/scale[col])
                mapped=torch.stack(components)[None,:]
                result,_=f.network.forward(mapped,conditional=f.context(ctx,1))
                return result[0]
            jac=torch.autograd.functional.jacobian(composite,torch.as_tensor(point))
            numerical_det=float(torch.linalg.slogdet(jac)[1])
            with torch.no_grad():_,fd=f.network.forward(torch.as_tensor(f.to_network(point[None,:])),conditional=f.context(ctx,1))
            determinant_errors.append(abs(numerical_det-float(fd[0])-f.log_coordinate_jacobian(point[None,:])[0]))
        composite_error=max(determinant_errors)
        assert composite_error<gates['autodiff_composite_logdet_absolute'],(channel,composite_error)
        samples=f.from_latent(latent,ctx);saved[channel]=samples
        physical=np.isfinite(samples).all(axis=1)&(samples[:,0]>0)&(samples[:,1]>0)&(samples[:,1]<=1)&(abs(samples[:,2])<=1)&(samples[:,3]>=0)
        report['checks'][channel]=dict(state_keys=f.state_keys,raw_network_context=ctx,
            coordinate_roundtrip_relative=map_error,coordinate_log_jacobian_absolute=coordinate_error,
            latent_roundtrip_absolute=inverse_error,logdet_cancellation_absolute=det_error,
            density_change_of_variables_absolute=density_error,autodiff_composite_logdet_absolute=composite_error,
            original_sampled_latent_density_difference=original_density_error,
            original_sampled_latent_density_gate_passed=bool(original_density_error<gates['density_change_of_variables_absolute']),
            sample_rows=count,rows_outside_true_physical_domain=int((~physical).sum()),
            learned_coordinate_support=dict(chirp_mass_upper=float(f.upper[0]),mass_ratio_upper=float(f.upper[1]),redshift_upper=float(f.upper[3])))
        report['samples'][channel]=dict(quantiles_01_50_99=np.quantile(samples,[.01,.5,.99],axis=0).tolist(),
                                       negative_effective_spin_sample_fraction=float((samples[:,2]<0).mean()),
                                       scope='Conditional unregularized emulator draws at raw context; not a fitted channel fraction or true synthesis probability certification.')
        print(f'{channel}: strict tensors and coordinate/network Jacobians verified; {count} emulator draws',flush=True)
    with h5py.File(PRODUCTS/'emulation_samps.hdf5','r') as h:
        keys=[x.decode() for x in h['flow_samps/block0_items'][:]]
        assert keys==['mchirp','q','chieff','z']
        reference_data=h['flow_samps/block0_values'][:]
    f=ChannelFlow('CE');latent=np.random.default_rng(design['sampling']['seed']).normal(size=(count,4))
    candidates={}
    for name,ctx in [('raw_alpha',[0.,2.]),('log_alpha',[0.,float(np.log(2.))])]:
        draw=f.from_latent(latent,ctx)
        distances=[float(ks_2samp(draw[:,i],reference_data[:,i]).statistic) for i in range(4)]
        candidates[name]=dict(raw_context=ctx,marginal_ks_distances=distances,maximum_marginal_ks=max(distances),
                              declared_diagnostic_within_003=bool(max(distances)<=.03))
    report['published_example_check']=dict(reference_rows=len(reference_data),candidates=candidates,
        scope='Exploratory association of the published emulation plotting file with the paper Fig1 CE example. An empirical match is not recovery of training execution or an independent unseen-simulation test.')
    np.savez_compressed(BASE/'runtime-samples.npz',**saved)
    report['sample_cache_sha256']=sha(BASE/'runtime-samples.npz')
    report['all_original_gates_passed']=all(c['original_sampled_latent_density_gate_passed'] for c in report['checks'].values())
    if not report['all_original_gates_passed']:report['status']='qualified_represented_point_checks_pass_original_draw_gate_failed'
    report['numerical_diagnosis']='The original sampled-latent density identity fails its fixed1e-8 gate for NSC near q=0.9999999959553232. Float64 inverse-logistic rounding moves the remapped coordinate by4.93e-10 and latent values by5.83e-8. The represented-point change-of-variables check passes at the same gate. Original failure preserved; no gate was widened and no sample removed.'
    report['limits']='Coordinate/flow consistency only, with the explicit original-draw boundary failure above. GW event prior reduction, spin-conditionals for selection, historical thresholds, density regularization and population-tail validation must precede channel-fraction claims. No observed-event fit here.'
    (BASE/'runtime-results.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=report['status'],published_example_check=report['published_example_check'],physical_domain_rows={k:v['rows_outside_true_physical_domain'] for k,v in report['checks'].items()}),indent=2))


if __name__=='__main__':
    try:main()
    except Exception as error:
        path=BASE/'runtime-failures.json';history=json.loads(path.read_text()) if path.exists() else []
        history.append(dict(tool_sha256=sha(Path(__file__)),error_type=type(error).__name__,message=str(error)))
        path.write_text(json.dumps(history,indent=2)+'\n',encoding='utf-8')
        raise
