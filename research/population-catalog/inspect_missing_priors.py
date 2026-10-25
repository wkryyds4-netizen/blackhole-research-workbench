"""Read relevant raw configuration declarations for the six analytic gaps."""
from pathlib import Path
import hashlib
import json
import re
import h5py
from inspect_catalog import plain
from http_hdf import RangeReader

BASE=Path(__file__).resolve().parent


def main():
    summary=json.loads((BASE/'catalog-summary.json').read_text());plan=json.loads((BASE/'source-plan.json').read_text())
    output=BASE/'missing-prior-configurations.json';rows=json.loads(output.read_text()) if output.exists() else []
    indices=[r['index'] for r in summary['events'] if not r['full_core_prior_strings_present']]
    for index in indices:
        if any(r['index']==index and r.get('configuration_schema',1)>=2 for r in rows):continue
        event=next(e for e in plan['events'] if e['index']==index);s=event['source']
        reader=RangeReader(s['url'],s['bytes'],cache_key=s['md5'],block_size=1024*1024,min_interval=3)
        try:
            with h5py.File(reader,'r') as h:
                group=h['C01:IMRPhenomXPHM'];paths=[];selected={};version={}
                if 'config_file' in group:
                    def visit(name,item):
                        if not isinstance(item,h5py.Dataset):return
                        paths.append(name)
                        if item.size<=1024 and re.search(r'prior|mass|chirp|mchirp|ratio|comp-|spin|tilt|dist|freq|fref|f-ref|f_ref|q-min|q-max|eta|approx',name,re.I):
                            selected[name]=plain(item[()])
                    group['config_file'].visititems(visit)
                if 'version' in group:
                    def versions(name,item):
                        if isinstance(item,h5py.Dataset) and item.size<=64:version[name]=plain(item[()])
                    if isinstance(group['version'],h5py.Dataset):version['value']=plain(group['version'][()])
                    else:group['version'].visititems(versions)
                metadata={}
                for key in ['f_ref','f_low','approximant']:
                    path='meta_data/meta_data/'+key
                    if path in group:metadata[key]=plain(group[path][()])
            row=dict(configuration_schema=2,index=index,event=event['event'],source=s,configuration_paths=paths,selected_prior_configuration=selected,
                analysis_version=version,metadata=metadata,
                access=dict(requests=reader.requests,network_bytes=reader.network_bytes,cached_bytes=reader.cached_bytes,ranges=reader.ranges),
                limits='Recorded configuration declarations only. Omitted flags may use historical sampler defaults; no inferred prior is asserted here.')
            rows=[r for r in rows if r['index']!=index]+[row];output.write_text(json.dumps(rows,indent=2),encoding='utf-8')
            print(json.dumps(dict(index=index,event=event['event'],selected=selected,version=version,metadata=metadata)),flush=True)
        finally:reader.close()


if __name__=='__main__':main()
