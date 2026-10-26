"""Physical-channel preview at an explicitly data-supported context mapping.

No observed events are fitted. Raw and physical-domain-valid sample counts are
both recorded; no silent clipping, folding or renormalization of q is applied.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import torch
from channel_flow import ChannelFlow,CHANNELS

BASE=Path(__file__).resolve().parent


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    torch.set_num_threads(4)
    runtime=json.loads((BASE/'runtime-results.json').read_text())
    candidates=runtime['published_example_check']['candidates']
    assert candidates['log_alpha']['declared_diagnostic_within_003']
    assert not candidates['raw_alpha']['declared_diagnostic_within_003']
    # Chosen after the mapping diagnostic but before generating these previews.
    # Training natal-spin grid endpoints and the published alpha=2 example;
    # labels are a data-supported interpretation, not a recovered invocation.
    design=dict(seed=20261028,samples_per_case=32768,birth_spins=[0.,.1,.5],ce_efficiency=2.,
                ce_context='[chi_birth,ln(alpha_CE)]',other_context='[chi_birth]',
                scope='Interpolation-grid preview, not GW inference or a new synthesis training-validation set.')
    (BASE/'preview-design.json').write_text(json.dumps(design,indent=2)+'\n',encoding='utf-8')
    rows=[];saved={}
    for ci,channel in enumerate(CHANNELS):
        f=ChannelFlow(channel)
        z=np.random.default_rng(design['seed']+ci).normal(size=(design['samples_per_case'],4))
        for spin in design['birth_spins']:
            ctx=[spin,np.log(design['ce_efficiency'])] if channel=='CE' else [spin]
            samples=f.from_latent(z,ctx)
            good=np.isfinite(samples).all(axis=1)&(samples[:,0]>0)&(samples[:,1]>0)&(samples[:,1]<=1)&(abs(samples[:,2])<=1)&(samples[:,3]>=0)
            key=f'{channel}_birth{spin:g}';saved[key]=samples
            m=samples[:,0]*(1+samples[:,1])**.2/samples[:,1]**.6
            domain=good&(m<=100)&(m*samples[:,1]>=3)&(samples[:,3]<=1.9)&(abs(samples[:,2])<.99)
            rows.append(dict(channel=channel,birth_spin=spin,raw_context=list(map(float,ctx)),rows=len(samples),
                             raw_negative_effective_spin_fraction=float((samples[:,2]<0).mean()),
                             physical_domain_valid_rows=int(good.sum()),physical_domain_invalid_fraction=float(1-good.mean()),
                             restricted_example_domain_rows=int(domain.sum()),
                             restricted_example_domain_fraction=float(domain.mean()),
                             quantiles_01_50_99=np.quantile(samples,[.01,.5,.99],axis=0).tolist()))
        print(f'{channel}: three natal-spin contexts sampled',flush=True)
    np.savez_compressed(BASE/'channel-preview-samples.npz',**saved)
    report=dict(schema='blackhole-physical-channel-context-preview/1',design=design,cases=rows,
                input_sha256={p.name:sha(p) for p in [BASE/'runtime-results.json',BASE/'channel_flow.py',Path(__file__)]},
                sampling_cache_sha256=sha(BASE/'channel-preview-samples.npz'),
                limits='Raw emulator draws. Intrinsic prediction only with the stated context interpretation; no new channel branching inference, selection-corrected estimate, or trusted missing-tail bound. Invalid q is retained and counted. Domain fractions are emulator MC summaries, not certified physical support probabilities.')
    (BASE/'channel-preview.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps([dict(channel=r['channel'],birth_spin=r['birth_spin'],negative_fraction=r['raw_negative_effective_spin_fraction'],invalid_fraction=r['physical_domain_invalid_fraction']) for r in rows],indent=2))


if __name__=='__main__':main()
