"""Pinned small channel products, checksummed before safe archive inspection."""
from pathlib import Path,PurePosixPath
from concurrent.futures import ThreadPoolExecutor,as_completed
import argparse
import hashlib
import json
import tarfile
import requests

BASE=Path(__file__).resolve().parent
FOLDER=BASE/'sources'


def sha(path,method='sha256'):
    with path.open('rb') as stream:return hashlib.file_digest(stream,method).hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--record',type=int,default=14967688)
    parser.add_argument('--files',nargs='+',default=['flow_models.tar.gz','plot_data.tar.gz'])
    args=parser.parse_args()
    meta=json.loads((FOLDER/f'zenodo-{args.record}.json').read_text())
    assert meta['id']==args.record
    outputs=[]
    for name in args.files:
        assert '/' not in name and '\\' not in name and name.endswith('.tar.gz')
        entry=next(x for x in meta['files'] if x['key']==name)
        assert entry['size']<=130000000,'This fetcher deliberately selects small products first.'
        path=FOLDER/f'{args.record}-{name}';size=entry['size'];step=4*1024*1024
        starts=list(range(0,size,step));chunks=FOLDER/f'{args.record}-{name}.chunks'
        if not path.exists():
            chunks.mkdir(exist_ok=True)
            def fetch(start):
                end=min(start+step,size)-1;part=chunks/f'{start}.part'
                if part.exists() and part.stat().st_size==end-start+1:return
                for attempt in range(3):
                    try:
                        with requests.get(entry['links']['self'],headers={'Range':f'bytes={start}-{end}'},stream=True,timeout=(15,40)) as response:
                            response.raise_for_status()
                            assert response.status_code==206 and response.headers.get('Content-Range')==f'bytes {start}-{end}/{size}'
                            raw=b''.join(response.iter_content(65536))
                        assert len(raw)==end-start+1
                        part.write_bytes(raw);return
                    except (requests.RequestException,AssertionError):
                        if attempt==2:raise
            with ThreadPoolExecutor(max_workers=4) as pool:
                for n,job in enumerate(as_completed([pool.submit(fetch,start) for start in starts]),1):
                    job.result();print(f'{name}: {n}/{len(starts)} chunks',flush=True)
            assembled=path.with_suffix('.assembled')
            with assembled.open('wb') as stream:
                for start in starts:stream.write((chunks/f'{start}.part').read_bytes())
            assert assembled.stat().st_size==size and 'md5:'+sha(assembled,'md5')==entry['checksum']
            assembled.replace(path)
        assert path.stat().st_size==size and 'md5:'+sha(path,'md5')==entry['checksum']
        record=dict(file=path.name,url=entry['links']['self'],record=args.record,bytes=size,
                    official_md5=entry['checksum'],sha256=sha(path),license=meta['metadata'].get('license'))
        inventory=[]
        with tarfile.open(path,'r:gz') as archive:
            for member in archive:
                p=PurePosixPath(member.name)
                assert not p.is_absolute() and '..' not in p.parts
                assert not member.issym() and not member.islnk()
                inventory.append(dict(name=member.name,bytes=member.size,regular_file=member.isfile()))
                assert len(inventory)<10000
        outputs.append(dict(source=record,members=inventory))
        manifest=FOLDER/'manifest.json';old=json.loads(manifest.read_text())
        manifest.write_text(json.dumps([x for x in old if x['file']!=path.name]+[record],indent=2)+'\n',encoding='utf-8')
        print(json.dumps(outputs[-1],indent=2),flush=True)
    (BASE/f'archive-inventory-{args.record}.json').write_text(json.dumps(outputs,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
