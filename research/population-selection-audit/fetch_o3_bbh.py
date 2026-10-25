"""Fetch one fixed O3 BBH summary, with bounded resumable ranges and checksums."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import time
import requests

BASE=Path(__file__).resolve().parent
NAME='endo3_bbhpop-LIGO-T2100113-v12.hdf5'


def main():
    folder=BASE/'sources';meta=json.loads((folder/'zenodo-7890437.json').read_text())
    assert meta['id']==7890437
    entry=next(x for x in meta['files'] if x['key']==NAME)
    assert entry['size']==88842240 and entry['checksum']=='md5:dd6ba479ae01470320a3bdc68dbb3efc'
    dest=folder/NAME;size=4*1024*1024;starts=list(range(0,entry['size'],size))
    chunks=folder/(NAME+'.chunks')
    if not dest.exists():
        chunks.mkdir(exist_ok=True)
        def fetch(start):
            end=min(start+size,entry['size'])-1;part=chunks/f'{start}.part'
            if part.exists() and part.stat().st_size==end-start+1:return
            for attempt in range(5):
                try:
                    r=requests.get(entry['links']['self'],headers={'Range':f'bytes={start}-{end}'},timeout=(20,45))
                    r.raise_for_status()
                    assert r.status_code==206 and r.headers.get('Content-Range')==f'bytes {start}-{end}/{entry["size"]}'
                    assert len(r.content)==end-start+1
                    part.write_bytes(r.content);return
                except (requests.RequestException,AssertionError):
                    if attempt==4:raise
                    time.sleep(min(2**attempt,8))
        with ThreadPoolExecutor(max_workers=4) as pool:
            for n,job in enumerate(as_completed([pool.submit(fetch,s) for s in starts]),1):
                job.result();print(f'O3 BBH: {n}/{len(starts)} chunks',flush=True)
        temporary=dest.with_suffix('.part')
        with temporary.open('wb') as out:
            for start in starts:out.write((chunks/f'{start}.part').read_bytes())
        assert temporary.stat().st_size==entry['size'] and 'md5:'+hashlib.md5(temporary.read_bytes()).hexdigest()==entry['checksum']
        temporary.replace(dest)
    content=dest.read_bytes()
    assert len(content)==entry['size'] and 'md5:'+hashlib.md5(content).hexdigest()==entry['checksum']
    record=dict(file=NAME,url=entry['links']['self'],bytes=len(content),official_md5=entry['checksum'],sha256=hashlib.sha256(content).hexdigest())
    target=folder/'manifest.json';previous=json.loads(target.read_text())
    target.write_text(json.dumps([r for r in previous if r['file']!=NAME]+[record],indent=2),encoding='utf-8')
    # Only remove the explicit range files after full-file verification.
    if chunks.exists():
        for start in starts:
            part=chunks/f'{start}.part'
            if part.exists():part.unlink()
        if not any(chunks.iterdir()):chunks.rmdir()
    print(json.dumps(record,indent=2),flush=True)


if __name__=='__main__':main()
