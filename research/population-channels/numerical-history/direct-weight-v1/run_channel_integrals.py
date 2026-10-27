"""Full69 fixed-context channel integrals, with resumable hashed checkpoints.

Never runs root model.py and never changes earlier frozen evidence. All rows
and events are retained. A terminated run resumes only completed event records.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import hashlib
import json
import sys
import time
import h5py
import numpy as np
import torch
from conditional_channel_model import ConditionalChannels,physical_domain,CHANNELS,FAMILIES

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
CAT=ROOT/'research/population-catalog'
SEL=ROOT/'research/population-selection-audit'
MIX=ROOT/'research/population-mixture'
sys.path.insert(0,str(SEL))
from run_selection_study import density_pair,physical,PARAMS


def sha(p):
    with p.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def save(p,d):
    tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8');tmp.replace(p)


def stats():
    return dict(rows=0,positive_rows=np.zeros((2,6),np.int64),total=np.zeros((2,6)),
                outer=np.zeros((2,6,6)),maximum=np.zeros((2,6)))


def add(s,w):
    s['rows']+=len(w);s['positive_rows']+=np.count_nonzero(w,axis=0)
    s['total']+=w.sum(axis=0);s['outer']+=np.einsum('nki,nkj->kij',w,w)
    s['maximum']=np.maximum(s['maximum'],w.max(axis=0))


def serial(s):
    return {k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in s.items()}


def normalizers(f,protocol):
    path=BASE/'channel-normalization.json'
    report=json.loads(path.read_text()) if path.exists() else dict(protocol=protocol,channels=[])
    assert report['protocol']==protocol,'Preserve/version outputs when protocol changes.'
    design=f.design['normalization'];n=design['draws_per_channel'];birth=f.design['context']['birth_spin']
    for ci,(name,flow) in enumerate(zip(CHANNELS,f.models)):
        if any(r['channel']==name for r in report['channels']):continue
        rng=np.random.default_rng(design['seed']+ci);counts=[0,0];invalid=0;rows=0
        context=[birth,np.log(f.design['context']['ce_efficiency'])] if name=='CE' else [birth]
        for start in range(0,n,8192):
            theta=flow.from_latent(rng.normal(size=(min(8192,n-start),4)),context)
            counts[int(start>=n//2)]+=int(physical_domain(theta).sum())
            invalid+=int((theta[:,1]>1).sum());rows+=len(theta)
        fraction=sum(counts)/n
        report['channels'].append(dict(channel=name,rows=rows,in_domain=sum(counts),fraction=fraction,
                                        half_fractions=[x/(n//2) for x in counts],raw_q_above_one=invalid,
                                        bernoulli_mc_standard_error=float(np.sqrt(fraction*(1-fraction)/n))))
        save(path,report);print(f'Normalizer {name}: {fraction:.8f}, {n} draws',flush=True)
    return report


def event_integrals(f,protocol):
    path=BASE/'channel-event-integrals.json'
    cache=json.loads(path.read_text()) if path.exists() else dict(protocol=protocol,events=[])
    assert cache['protocol']==protocol,'Preserve/version outputs when protocol changes.'
    frozen=json.loads((ROOT/'results/model.json').read_text())
    provenance=json.loads((CAT/'event-provenance.json').read_text())
    for pos,event in enumerate(frozen['events']):
        index=event['index'];info=next(r for r in provenance if r['index']==index)
        assert info['native_source_identity_resolved']
        source=CAT/f'event-{index}.npz';fingerprint=sha(source)
        assert fingerprint==info['cache_sha256']
        if any(r['index']==index and r['source_cache_sha256']==fingerprint for r in cache['events']):continue
        with np.load(source,allow_pickle=False) as bank:d={k:bank[k] for k in PARAMS}
        rows=len(d['mass_1']);summary=stats();started=time.monotonic()
        for start in range(0,rows,8192):
            part={k:v[start:start+8192] for k,v in d.items()};_,part['prior'],_=density_pair(part)
            add(summary,f.weights(part))
            if start and start%65536==0:print(f'Event {pos+1}/69 {info["event"]}: {min(start+8192,rows)}/{rows} rows',flush=True)
        record=dict(index=index,event=info['event'],source_cache_sha256=fingerprint,
                    **serial(summary),mean=(summary['total']/rows).tolist(),seconds=time.monotonic()-started)
        cache['events'].append(record);save(path,cache)
        print(f'Completed event {pos+1}/69 {info["event"]}: all {rows} rows retained',flush=True)
    ordered=[next(r for r in cache['events'] if r['index']==e['index']) for e in frozen['events']]
    assert len(ordered)==69 and sum(r['rows'] for r in ordered)==11180578
    return ordered


def select_stats(f,raw,masks,factor,label):
    union=np.logical_or.reduce(list(masks.values()));indices=np.flatnonzero(union)
    summaries={k:stats() for k in masks}
    for start in range(0,len(indices),8192):
        take=indices[start:start+8192];d=physical({k:v[take] for k,v in raw.items()})
        w=f.weights(d)*factor[take,None,None]
        for name,mask in masks.items():
            selected=w[mask[take]]
            if len(selected):add(summaries[name],selected)
        if start and start%65536==0:print(f'Selection {label}: {min(start+8192,len(indices))}/{len(indices)} detected rows',flush=True)
    return {k:serial(v) for k,v in summaries.items()}


def selections(f,protocol):
    path=BASE/'channel-selection-integrals.json'
    cache=json.loads(path.read_text()) if path.exists() else dict(protocol=protocol,banks={})
    assert cache['protocol']==protocol
    provenance=json.loads((SEL/'o3-provenance.json').read_text())
    assert sha(SEL/'o3-bank.npz')==provenance['cache_sha256']
    with np.load(SEL/'o3-bank.npz',allow_pickle=False) as bank:
        full={k[5:]:bank[k] for k in bank.files if k.startswith('full_')}
        old={k[9:]:bank[k] for k in bank.files if k.startswith('combined_') and k!='combined_o3_to_full'}
    search=sorted(k for k in full if k.startswith('ifar_'))
    if 'early' not in cache:
        early=(old['name']!=b'o3')&(np.logical_or.reduce([old[k]>1 for k in search])|(old['optimal_snr_net']>10))
        cache['early']=select_stats(f,old,{'early':early},old['mixture_weight']/provenance['combined_attributes']['total_generated'],'O1/O2')['early']
        save(path,cache)
    time_fraction=provenance['attributes']['analysis_time_s']/provenance['combined_attributes']['analysis_time_s']
    for name in f.design['selection']['banks']:
        if name in cache['banks']:continue
        if name=='pure_bbh':raw=full;attributes=provenance['attributes']
        else:
            source=json.loads((MIX/'download.json').read_text());hdf=MIX/'sources'/source['file']
            assert sha(hdf)==source['sha256']
            keys=['mass1_source','mass2_source','redshift','sampling_pdf','mixture_weight']+[f'spin{i}{axis}' for i in [1,2] for axis in 'xyz']+search
            with h5py.File(hdf,'r') as handle:
                raw={k:handle['injections'][k][:] for k in keys}
                attributes={k:float(handle.attrs[k]) for k in ['analysis_time_s','total_generated']}
        union=np.logical_or.reduce([raw[k]>1 for k in search])
        matched=np.logical_or.reduce([raw[k]>1 for k in search if k!='ifar_cwb'])
        factor=raw['mixture_weight']*time_fraction/attributes['total_generated']
        cache['banks'][name]=dict(accepted_rows=len(union),total_generated=attributes['total_generated'],
                                  policies=select_stats(f,raw,dict(matched_filter=matched,with_cwb=union),factor,name))
        cache['o3_time_fraction']=time_fraction;save(path,cache)
        print(f'Completed selection {name}',flush=True)
    return cache


def main():
    torch.set_num_threads(4)
    bridge=json.loads((BASE/'fast-measure-validation.json').read_text());assert bridge['status']=='pass'
    for file,digest in bridge['input_sha256'].items():assert sha(BASE/file)==digest
    for folder,result in [('population-catalog','native-catalog'),('population-mixture',None)]:
        directory=ROOT/'research'/folder
        validation=directory/(result+'-validation.json' if result else 'validation.json')
        product=directory/(result+'-results.json' if result else 'results.json')
        assert json.loads(validation.read_text())['result_sha256']==sha(product)
    pieces=[Path(__file__),BASE/'conditional_channel_model.py',BASE/'fast_measure.py',BASE/'channel_flow.py',
            BASE/'likelihood-design.json',BASE/'fast-measure-validation.json',BASE/'product-inspection.json',
            ROOT/'research/population-stage5/run_tail_study.py',ROOT/'research/population-stage4/recover_prior.py',
            SEL/'run_selection_study.py',CAT/'event-provenance.json']
    hashes={p.relative_to(ROOT).as_posix():sha(p) for p in pieces}
    protocol=hashlib.sha256(json.dumps(hashes,sort_keys=True).encode()).hexdigest()
    save(BASE/'channel-integral-protocol.json',dict(protocol=protocol,input_sha256=hashes))
    f=ConditionalChannels()
    norm=normalizers(f,protocol)
    used=event_integrals(f,protocol)
    selection=selections(f,protocol)
    save(BASE/'channel-integrals-complete.json',dict(status='integrals_computed_not_yet_independently_validated',protocol=protocol,
          events=len(used),rows=sum(r['rows'] for r in used),channels=CHANNELS,hidden_spin_conditionals=FAMILIES,
          input_sha256={p.name:sha(p) for p in [BASE/'channel-normalization.json',BASE/'channel-event-integrals.json',BASE/'channel-selection-integrals.json']},
          scope='Fixed-context restricted-domain conditional calculation. No completed branching-fraction claim.'))
    print('All69 event, normalization and selection integrals saved; independent validation and posterior analysis remain.',flush=True)


if __name__=='__main__':main()
