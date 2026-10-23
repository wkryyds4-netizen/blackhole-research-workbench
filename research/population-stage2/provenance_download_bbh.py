"""Download the single explicitly identified public BBH source candidate."""
from pathlib import Path
import hashlib,json,time
import requests
from concurrent.futures import ThreadPoolExecutor,as_completed
BASE=Path(__file__).resolve().parent
OUT=BASE/'source-data';OUT.mkdir(exist_ok=True)
name='o1+o2+o3_bbhpop_real+semianalytic-LIGO-T2100377-v2.hdf5'
meta=json.loads((BASE/'provenance-sources/zenodo-7890398.json').read_text())
record=next(f for f in meta['files'] if f['key']==name)
dest=OUT/name
url=record['links']['self'].replace('/7890398/','/5636816/')
if not dest.exists():
    def fetch(start):
        end=min(start+2*1024*1024,record['size'])-1
        for attempt in range(4):
            try:
                with requests.get(url,headers={'Range':f'bytes={start}-{end}'},stream=True,timeout=(20,45)) as r:
                    if r.status_code!=206:raise ValueError(f'Range request returned {r.status_code}')
                    assert r.headers.get('Content-Range')==f'bytes {start}-{end}/{record["size"]}'
                    data=r.content
                    assert len(data)==end-start+1
                    return start,data
            except (requests.RequestException,ValueError,AssertionError) as exc:
                if attempt==3:raise
                time.sleep(1)
    with dest.with_suffix('.part').open('wb') as output, ThreadPoolExecutor(max_workers=4) as pool:
        n=0
        for task in as_completed([pool.submit(fetch,start) for start in range(0,record['size'],2*1024*1024)]):
            start,data=task.result();output.seek(start);output.write(data);n+=len(data)
            print('downloaded',n,'/',record['size'],flush=True)
    assert hashlib.md5(dest.with_suffix('.part').read_bytes()).hexdigest()==record['checksum'].split(':')[1]
    dest.with_suffix('.part').replace(dest)
data=dest.read_bytes()
assert len(data)==record['size']
assert hashlib.md5(data).hexdigest()==record['checksum'].split(':')[1]
manifest=dict(file=str(dest.relative_to(BASE)),url=url,bytes=len(data),md5=hashlib.md5(data).hexdigest(),sha256=hashlib.sha256(data).hexdigest(),retrieved='2026-10-24')
(BASE/'provenance-download.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print(json.dumps(manifest,indent=2))
