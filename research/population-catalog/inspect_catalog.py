"""Bounded metadata-only inspection of fixed full-catalog candidate HDF files."""
from pathlib import Path
import argparse
import hashlib
import json
import h5py
import numpy as np
from http_hdf import RangeReader

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
META=['meta_data/meta_data/f_ref','meta_data/meta_data/reference_frequency',
      'config_file/config/reference-frequency','config_file/config/reference_frequency',
      'meta_data/meta_data/backward_spin_evolution','meta_data/sampler/pe_algorithm',
      'meta_data/meta_data/approximant']
SPINS=['a_1','a_2','tilt_1','tilt_2','cos_tilt_1','cos_tilt_2','phi_12',
       'spin_1x','spin_1y','spin_1z','spin_2x','spin_2y','spin_2z',
       'cos_tilt_1_infinity_only_prec_avg','cos_tilt_2_infinity_only_prec_avg']
PRIOR_FIELDS=['mass_1','mass_2','chirp_mass','mass_ratio','luminosity_distance','a_1','a_2','tilt_1','tilt_2','phi_12','phi_jl','theta_jn']


def plain(value):
    if isinstance(value,bytes):return value.decode('utf-8',errors='replace')
    if isinstance(value,np.ndarray):return plain(value.tolist())
    if isinstance(value,np.generic):return plain(value.item())
    if isinstance(value,list):return [plain(v) for v in value]
    return value


def inspect(h):
    labels=[]
    for label,g in h.items():
        if not isinstance(g,h5py.Group) or 'posterior_samples' not in g:continue
        samples=g['posterior_samples']
        if isinstance(samples,h5py.Dataset):
            fields=list(samples.dtype.names or []);shape=list(samples.shape);chunks=plain(samples.chunks)
        else:
            fields=list(samples.keys());shape=None;chunks=None
        meta={}
        for key in META:
            if key in g:
                obj=g[key]
                assert isinstance(obj,h5py.Dataset) and obj.size<=1024,key
                meta[key]=plain(obj[()])
        priors={}
        for key in PRIOR_FIELDS:
            path='priors/analytic/'+key
            if path in g:
                obj=g[path]
                assert isinstance(obj,h5py.Dataset) and obj.size<=1024,path
                priors[key]=plain(obj[()])
        labels.append(dict(label=label,shape=shape,chunks=chunks,spin_fields_present=[k for k in SPINS if k in fields],
            all_fields=fields,metadata=meta,analytic_prior_strings=priors))
    return labels


def main():
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--indices',nargs='+',type=int)
    group.add_argument('--all',action='store_true')
    parser.add_argument('--local',action='store_true',help='Use checksum-verified local file when available; fail if missing')
    args=parser.parse_args();plan=json.loads((BASE/'source-plan.json').read_text())
    output=BASE/'metadata-results.json';results=json.loads(output.read_text()) if output.exists() else []
    for index in ([e['index'] for e in plan['events']] if args.all else args.indices):
        e=next(v for v in plan['events'] if v['index']==index);s=e['source']
        if any(r['index']==index and r['source_md5']==s['md5'] and r.get('inspection_schema',1)>=2 and r['mode']==('local' if args.local else 'remote') for r in results):
            print(f'Skip completed metadata index {index}',flush=True);continue
        if args.local:
            candidates=[folder/s['file'] for folder in [ROOT/'research/population-stage7/source-data',ROOT/'research/population-stage5/source-data',BASE/'source-data']]
            existing=[p for p in candidates if p.exists()];assert len(existing)==1,s['file']
            path=existing[0];assert 'md5:'+hashlib.md5(path.read_bytes()).hexdigest()==s['md5']
            with h5py.File(path,'r') as h:labels=inspect(h)
            access=dict(local_path=str(path.relative_to(ROOT)).replace('\\','/'),whole_file_md5_verified=True)
        else:
            reader=RangeReader(s['url'],s['bytes'],cache_key=s['md5'])
            try:
                with h5py.File(reader,'r') as h:labels=inspect(h)
                access=dict(requests=reader.requests,network_bytes=reader.network_bytes,object_identity_header=reader.identity,
                    cached_bytes=reader.cached_bytes,
                    ranges=reader.ranges,whole_file_md5_verified=False,
                    limits='Exact Content-Range/size/header checks and per-block hashes; no full-file checksum for remote metadata reads.')
            finally:reader.close()
        row=dict(inspection_schema=2,index=index,event=e['event'],source_file=s['file'],source_md5=s['md5'],mode='local' if args.local else 'remote',labels=labels,access=access)
        results=[v for v in results if not (v['index']==index and v['mode']==row['mode'])]+[row]
        output.write_text(json.dumps(results,indent=2),encoding='utf-8')
        print(json.dumps(dict(index=index,event=e['event'],labels=[dict(label=v['label'],metadata=v['metadata'],phi12='phi_12' in v['spin_fields_present']) for v in labels],access={k:v for k,v in access.items() if k!='ranges'})),flush=True)


if __name__=='__main__':main()
