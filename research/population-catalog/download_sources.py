"""Bounded, resumable full-file downloads for explicitly selected catalog events."""
from pathlib import Path
import argparse
import hashlib
import json
import time
import requests

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]


class CatalogRateLimit(IOError):
    pass


def md5(path):
    digest=hashlib.md5()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):digest.update(block)
    return 'md5:'+digest.hexdigest()


def download_one(event):
    index=event['index'];s=event['source'];name=s['file'];folder=BASE/'source-data';folder.mkdir(exist_ok=True)
    assert Path(name).name==name and name.endswith('.h5') and s['bytes']<1024**3
    candidates=[folder/name,ROOT/'research/population-stage7/source-data'/name,ROOT/'research/population-stage5/source-data'/name]
    path=next((p for p in candidates if p.exists()),candidates[0]);reused=path.exists();part=path.with_suffix('.partial')
    if not reused:
        for attempt in range(4):
            start=part.stat().st_size if part.exists() else 0
            assert start<=s['bytes']
            if start==s['bytes']:break
            try:
                headers={'Accept-Encoding':'identity'}
                if start:headers['Range']=f'bytes={start}-{s["bytes"]-1}'
                print(f'Download index {index}: {start}/{s["bytes"]}',flush=True)
                with requests.get(s['url'],headers=headers,stream=True,timeout=(10,45)) as response:
                    if response.status_code==429:raise CatalogRateLimit('HTTP429; Retry-After='+str(response.headers.get('Retry-After')))
                    response.raise_for_status()
                    if start:assert response.status_code==206 and response.headers.get('Content-Range')==f'bytes {start}-{s["bytes"]-1}/{s["bytes"]}'
                    else:assert response.status_code==200
                    count=start;last=count
                    with part.open('ab' if start else 'wb') as out:
                        for chunk in response.iter_content(1024*1024):
                            count+=len(chunk);assert count<=s['bytes'];out.write(chunk)
                            if count-last>=64*1024*1024:print(f'Index {index}: {count}/{s["bytes"]}',flush=True);last=count
                assert count==s['bytes'];break
            except requests.RequestException:
                if attempt==3:raise
                time.sleep(min(2**attempt,8))
        assert part.stat().st_size==s['bytes']
        checked=part
    else:checked=path
    assert checked.stat().st_size==s['bytes']
    md=hashlib.md5();sha=hashlib.sha256()
    with checked.open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):md.update(block);sha.update(block)
    assert 'md5:'+md.hexdigest()==s['md5'],'Source checksum mismatch; retained file for diagnosis'
    if not reused:part.replace(path)
    return dict(index=index,event=event['event'],source=s,sha256=sha.hexdigest(),whole_file_md5_verified=True,
        local_file=str(path.relative_to(ROOT)).replace('\\','/'),reused_local_file=reused)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--indices',nargs='+',type=int,required=True);args=parser.parse_args()
    plan=json.loads((BASE/'source-plan.json').read_text());output=BASE/'download-results.json';rows=json.loads(output.read_text()) if output.exists() else []
    for index in args.indices:
        row=download_one(next(e for e in plan['events'] if e['index']==index));rows=[r for r in rows if r['index']!=index]+[row]
        temporary=output.with_suffix('.tmp');temporary.write_text(json.dumps(rows,indent=2),encoding='utf-8');temporary.replace(output)
        print(json.dumps(dict(index=index,sha256=row['sha256'],verified=True)),flush=True)


if __name__=='__main__':main()
