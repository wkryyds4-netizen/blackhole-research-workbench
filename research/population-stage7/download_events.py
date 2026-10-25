"""Bounded resumable downloads of the eight preselected official PE files."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import hashlib
import json
import time
import requests

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]


def hashes(path):
    md5=hashlib.md5();sha=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):
            md5.update(block);sha.update(block)
    return md5.hexdigest(),sha.hexdigest()


def download(entry,event,index):
    name=entry['file'];assert Path(name).name==name and name.endswith('.h5')
    folder=BASE/'source-data';folder.mkdir(parents=True,exist_ok=True)
    dest=folder/name;chunks=folder/(name+'.chunks')
    size=4*1024*1024;starts=list(range(0,entry['bytes'],size))
    if not dest.exists():
        chunks.mkdir(exist_ok=True)
        def fetch(start):
            end=min(start+size,entry['bytes'])-1;part=chunks/f'{start}.part'
            if part.exists() and part.stat().st_size==end-start+1:return
            for attempt in range(5):
                try:
                    r=requests.get(entry['url'],headers={'Range':f'bytes={start}-{end}'},timeout=(20,45))
                    r.raise_for_status()
                    assert r.status_code==206
                    assert r.headers.get('Content-Range')==f'bytes {start}-{end}/{entry["bytes"]}'
                    assert len(r.content)==end-start+1
                    part.write_bytes(r.content);return
                except (requests.RequestException,AssertionError):
                    if attempt==4:raise
                    time.sleep(min(2**attempt,8))
        print(f'Start {event} {name} ({entry["bytes"]} bytes)',flush=True)
        with ThreadPoolExecutor(max_workers=4) as pool:
            jobs=[pool.submit(fetch,start) for start in starts]
            for n,job in enumerate(as_completed(jobs),1):
                job.result()
                if n%10==0 or n==len(starts):print(f'{event} {"nocosmo" if "nocosmo" in name else "cosmo"}: {n}/{len(starts)} chunks',flush=True)
        temporary=dest.with_suffix('.part')
        with temporary.open('wb') as out:
            for start in starts:out.write((chunks/f'{start}.part').read_bytes())
        md5,sha=hashes(temporary)
        assert temporary.stat().st_size==entry['bytes'] and 'md5:'+md5==entry['md5']
        temporary.replace(dest)
    md5,sha=hashes(dest)
    assert dest.stat().st_size==entry['bytes'] and 'md5:'+md5==entry['md5']
    # Remove only this fixed file's explicit chunks after the full checksum passes.
    if chunks.exists():
        for start in starts:
            part=chunks/f'{start}.part'
            if part.exists():part.unlink()
        if not any(chunks.iterdir()):chunks.rmdir()
    result=dict(entry,event=event,index=index,sha256=sha)
    print(f'Verified {event} {name}: {sha}',flush=True)
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--indices',nargs='+',type=int,default=[19,20,25,50])
    args=parser.parse_args();assert set(args.indices)<=set([19,20,25,50])
    plan=json.loads((ROOT/'research/population-stage6/followup-source-plan.json').read_text())
    output=BASE/'downloads.json';records=json.loads(output.read_text()) if output.exists() else []
    for event in plan['events']:
        if event['index'] not in args.indices:continue
        # Finish the smaller cosmo file before its associated original-prior file.
        for entry in sorted(event['files'],key=lambda f:f['bytes']):
            result=download(entry,event['event'],event['index'])
            records=[r for r in records if r['file']!=result['file']]+[result]
            output.write_text(json.dumps(records,indent=2),encoding='utf-8')


if __name__=='__main__':main()
