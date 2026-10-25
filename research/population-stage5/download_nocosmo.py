"""Download one fixed 335 MB public PE file, with resumable bounded chunks."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import time
import requests

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
NAME='IGWN-GWTC2p1-v2-GW190805_211137_PEDataRelease_mixed_nocosmo.h5'

def main():
    meta=json.loads((ROOT/'research/population-stage3/sources/zenodo-6513631.json').read_text())
    rec=next(f for f in meta['files'] if f['key']==NAME)
    assert rec['size']==335235368
    dest=BASE/'source-data'/NAME
    chunks=BASE/'source-data/chunks';chunks.mkdir(parents=True,exist_ok=True)
    size=4*1024*1024
    url=rec['links']['self']
    if not dest.exists():
        def fetch(start):
            end=min(start+size,rec['size'])-1
            part=chunks/f'{start}.part'
            if part.exists() and part.stat().st_size==end-start+1:return start
            for attempt in range(4):
                try:
                    r=requests.get(url,headers={'Range':f'bytes={start}-{end}'},timeout=(20,45))
                    r.raise_for_status()
                    assert r.status_code==206
                    assert r.headers['Content-Range']==f'bytes {start}-{end}/{rec["size"]}'
                    assert len(r.content)==end-start+1
                    part.write_bytes(r.content)
                    return start
                except (requests.RequestException,AssertionError):
                    if attempt==3:raise
                    time.sleep(1)
        starts=list(range(0,rec['size'],size))
        with ThreadPoolExecutor(max_workers=4) as pool:
            tasks=[pool.submit(fetch,start) for start in starts]
            for n,job in enumerate(as_completed(tasks),1):
                job.result()
                if n%10==0 or n==len(starts):print(f'chunks {n}/{len(starts)}',flush=True)
        with dest.with_suffix('.part').open('wb') as output:
            for start in starts:output.write((chunks/f'{start}.part').read_bytes())
        data=dest.with_suffix('.part').read_bytes()
        assert len(data)==rec['size']
        assert 'md5:'+hashlib.md5(data).hexdigest()==rec['checksum']
        dest.with_suffix('.part').replace(dest)
    data=dest.read_bytes()
    assert len(data)==rec['size']
    assert 'md5:'+hashlib.md5(data).hexdigest()==rec['checksum']
    # Verified final file supersedes only this downloader's explicitly named chunks.
    for start in range(0,rec['size'],size):
        part=chunks/f'{start}.part'
        if part.exists():part.unlink()
    if not any(chunks.iterdir()):chunks.rmdir()
    result=dict(record=6513631,file=NAME,url=url,bytes=len(data),md5=hashlib.md5(data).hexdigest(),
                sha256=hashlib.sha256(data).hexdigest())
    (BASE/'download.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
