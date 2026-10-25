"""Paired event-bank and joint event/O3 perturbations, fixed69 observations."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
from scipy.integrate import cumulative_trapezoid

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
SEL=ROOT/'research/population-selection-audit'
CAT=ROOT/'research/population-catalog'
sys.path.insert(0,str(SEL))
from run_selection_study import init_model,density_pair,weights as native_weights
from run_control import weights as effective_weights


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def curve(grid,fs,means,selection):
    mix=fs[None,:,:]*means[:,:,0,None]+(1-fs[None,:,:])*means[:,:,1,None]
    with np.errstate(divide='ignore'):numerator=np.log(mix).sum(axis=0)
    efficiency=fs*selection[:,0,None]+(1-fs)*selection[:,1,None]
    logs=numerator-69*np.log(efficiency)
    p=np.exp(logs-logs.max(axis=1)[:,None]);p/=np.trapezoid(p,grid,axis=1)[:,None]
    cdf=cumulative_trapezoid(p,grid,axis=1,initial=0);cdf/=cdf[:,-1,None]
    q=np.array([np.interp([.05,.5,.95],c,grid) for c in cdf])
    return q,p,cdf,mix


def main():
    init_model();frozen=json.loads((ROOT/'results/model.json').read_text());sources=json.loads((CAT/'event-provenance.json').read_text())
    parent=json.loads((BASE/'paired-o3-results.json').read_text());validation=json.loads((BASE/'paired-o3-validation.json').read_text())
    assert validation['result_sha256']==sha(BASE/'paired-o3-results.json')
    with np.load(BASE/'paired-o3-checkpoints.npz',allow_pickle=False) as cp:
        a={k:cp[k] for k in cp.files if k in ['native_grid','effective_grid','native_fs','native_means','effective_means','full','selections','quantiles','full_quantiles','contrasts','analytic_contrast_covariance']}
    full_means=np.stack([a['native_means'],a['effective_means']],axis=1);assert full_means.shape==(69,2,6,2)
    draws=np.empty((512,69,2,6,2));covariances=[];details=[]
    for pos,event in enumerate(frozen['events']):
        index=event['index'];source=next(r for r in sources if r['index']==index);path=CAT/f'event-{index}.npz'
        assert sha(path)==source['cache_sha256']
        with np.load(path,allow_pickle=False) as bank:d={k:bank[k] for k in bank.files}
        _,prior,_=density_pair(d);d['prior']=prior
        w=np.stack([native_weights(d),effective_weights(d)],axis=1).reshape(len(prior),24)
        mu=w.mean(axis=0).reshape(2,6,2);wanted=full_means[pos]
        assert np.array_equal(mu==0,wanted==0)
        mask=wanted!=0;error=float(np.max(abs(mu[mask]/wanted[mask]-1)));assert error<1e-10
        m=len(w);n=6913;assert n<m
        covariance=(1-n/m)/n*np.cov(w,rowvar=False,ddof=1);covariances.append(covariance)
        rng=np.random.default_rng(np.random.SeedSequence([20261026,index]));fingerprint=hashlib.sha256()
        for rep in range(512):
            take=rng.choice(m,n,replace=False);fingerprint.update(take.astype('<i8').tobytes())
            draws[rep,pos]=w[take].mean(axis=0).reshape(2,6,2)
        details.append(dict(index=index,event=source['event'],full_rows=m,retained_rows=n,mean_regression_relative_error=error,
            all512_index_stream_sha256=fingerprint.hexdigest(),zero_weight_coordinates=np.flatnonzero(w.sum(axis=0)==0).tolist(),
            source_cache_sha256=source['cache_sha256']))
        print(f'Paired event bank {index}: {m} -> {n},512 shared draws',flush=True)
    covariances=np.array(covariances)
    grids=[a['native_grid'],a['effective_grid']];fractions=[a['native_fs'],np.broadcast_to(a['effective_grid'],(6,len(a['effective_grid'])))]
    # Derivatives in the actually sampled directions. If a whole component is
    # zero, its sample covariance is zero; its counterpart changes only a
    # per-event normalization. Both contribute no finite-bank shape variation.
    J=np.zeros((69,24,24));zero_directions=[]
    for p in range(2):
        for t in range(2):
            grid=grids[t];fs=fractions[t];q,density,cdf,mix=curve(grid,fs,full_means[:,t],a['full'][p,t])
            assert np.max(abs(q-a['full_quantiles'][p,t]))<1e-11
            for k in range(6):
                location=q[k,1];j=int(np.clip(np.searchsorted(grid,location)-1,0,len(grid)-2));slope=(cdf[k,j+1]-cdf[k,j])/(grid[j+1]-grid[j])
                out=np.ravel_multi_index((p,t,k),(2,2,6))
                for i in range(69):
                    if np.any(full_means[i,t,k]==0):
                        zero_directions.append(dict(policy=p,metric=t,model=k,event_position=i));continue
                    for c,basis in enumerate([fs[k],1-fs[k]]):
                        score=basis/mix[i,k]
                        integ=cumulative_trapezoid(density[k]*score,grid,initial=0)
                        value=-(np.interp(location,grid,integ)-.5*integ[-1])/slope
                        inp=np.ravel_multi_index((t,k,c),(2,6,2));J[i,out,inp]=value
                    pair=J[i,out,np.ravel_multi_index((t,k,0),(2,6,2)):np.ravel_multi_index((t,k,1),(2,6,2))+1]
                    assert abs(pair@full_means[i,t,k])<1e-10
    C=a['contrasts'];analytic_event=np.zeros((24,24))
    for i in range(69):analytic_event+=J[i]@covariances[i]@J[i].T
    analytic_event=C@analytic_event@C.T
    qs_event=[];qs_joint=[]
    for rep in range(512):
        event=[];joint=[]
        for p in range(2):
            ep=[];jp=[]
            for t in range(2):
                ep.append(curve(grids[t],fractions[t],draws[rep,:,t],a['full'][p,t])[0])
                jp.append(curve(grids[t],fractions[t],draws[rep,:,t],a['selections'][rep,p,t])[0])
            event.append(ep);joint.append(jp)
        qs_event.append(event);qs_joint.append(joint)
        if (rep+1)%64==0:print(f'Event/joint posterior propagation: {rep+1}/512',flush=True)
    qs_event=np.array(qs_event);qs_joint=np.array(qs_joint)
    nominal=a['full_quantiles'][...,1].reshape(24)@C.T
    ev=qs_event[...,1].reshape(512,24)@C.T
    jo=qs_joint[...,1].reshape(512,24)@C.T
    se=a['quantiles'][...,1].reshape(512,24)@C.T
    interaction=jo-ev-se+nominal
    linear_event=np.zeros((512,24))
    delta=draws.reshape(512,69,24)-full_means.reshape(1,69,24)
    for i in range(69):linear_event+=delta[:,i]@J[i].T
    linear_event=linear_event@C.T
    summaries=[]
    for k,row in enumerate(parent['contrasts']):
        summaries.append(dict(name=row['name'],full_bank_value_pp=float(nominal[k]),
            event_only_5_50_95_pp=np.quantile(ev[:,k],[.05,.5,.95]).tolist(),
            selection_only_5_50_95_pp=np.quantile(se[:,k],[.05,.5,.95]).tolist(),
            joint_5_50_95_pp=np.quantile(jo[:,k],[.05,.5,.95]).tolist(),
            event_only_std_pp=float(np.std(ev[:,k],ddof=1)),selection_only_std_pp=float(np.std(se[:,k],ddof=1)),joint_std_pp=float(np.std(jo[:,k],ddof=1)),
            linear_event_std_pp=float(np.sqrt(max(analytic_event[k,k],0))),
            linear_joint_std_pp=float(np.sqrt(max(analytic_event[k,k]+a['analytic_contrast_covariance'][k,k],0))),
            event_nonpositive_fraction=float(np.mean(ev[:,k]<=0)),joint_nonpositive_fraction=float(np.mean(jo[:,k]<=0)),
            interaction_5_50_95_pp=np.quantile(interaction[:,k],[.05,.5,.95]).tolist(),
            interaction_rms_pp=float(np.sqrt(np.mean(interaction[:,k]**2)))))
    report=dict(schema='blackhole-full69-joint-finite-bank-control/1',replicates=512,events=69,retained_rows_per_event=6913,
        event_seed=20261026,selection_seed=20261025,contrasts=summaries,event_banks=details,
        zero_component_sampling_directions=zero_directions,
        input_hashes=dict(paired_o3=sha(BASE/'paired-o3-results.json'),paired_o3_validation=sha(BASE/'paired-o3-validation.json'),design=sha(BASE/'paired-pe-design.json')),
        scope='Independent SRSWOR streams on each released event bank, shared across metrics/models/policies; paired with the specified independent O3 draw stream. Fixed69 observed events and O1/O2. Zero-component sampling covariance is not evidence of zero physical likelihood.',
        limits='Artificial finite-bank sensitivity only.6913-row draws do not reproduce historical cosmo rejection. PE-chain correlation, unobserved posterior/likelihood tails, real-campaign design and shared waveform systematics remain uncalibrated. Intervals are neither astrophysical credible intervals nor p-values.')
    (BASE/'paired-pe-results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    np.savez_compressed(BASE/'paired-pe-checkpoints.npz',event_means=draws,full_means=full_means,event_mean_covariances=covariances,
        quantiles_event=qs_event,quantiles_joint=qs_joint,contrasts_event=ev,contrasts_joint=jo,contrasts_selection=se,interaction=interaction,
        median_event_derivative=J,analytic_event_contrast_covariance=analytic_event,linear_event_contrast_perturbations=linear_event)
    print(json.dumps(dict(contrasts=summaries),indent=2))


if __name__=='__main__':main()
