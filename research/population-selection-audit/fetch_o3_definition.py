"""Read fixed O3 documentation and notebook; never execute downloaded code."""
from pathlib import Path
import hashlib
import json
import requests

BASE=Path(__file__).resolve().parent


def main():
    folder=BASE/'sources';folder.mkdir(parents=True,exist_ok=True)
    url='https://zenodo.org/api/records/7890437'
    path=folder/'zenodo-7890437.json'
    if not path.exists():
        r=requests.get(url,timeout=(15,40));r.raise_for_status();assert r.json()['id']==7890437
        path.write_bytes(r.content)
    meta=json.loads(path.read_text());assert meta['id']==7890437
    records=[dict(file=path.name,url=url,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),bytes=path.stat().st_size)]
    for name in ['o3-sensitivity-estimates.md','o3b-vt-dr.ipynb']:
        entry=next(f for f in meta['files'] if f['key']==name);out=folder/name
        if not out.exists():
            r=requests.get(entry['links']['self'],timeout=(15,40));r.raise_for_status()
            assert len(r.content)==entry['size'] and 'md5:'+hashlib.md5(r.content).hexdigest()==entry['checksum']
            out.write_bytes(r.content)
        content=out.read_bytes()
        assert len(content)==entry['size'] and 'md5:'+hashlib.md5(content).hexdigest()==entry['checksum']
        records.append(dict(file=name,url=entry['links']['self'],sha256=hashlib.sha256(content).hexdigest(),
                            bytes=len(content),official_md5=entry['checksum']))
    target=folder/'manifest.json';previous=json.loads(target.read_text())
    previous=[x for x in previous if x['file'] not in {r['file'] for r in records}]+records
    target.write_text(json.dumps(previous,indent=2),encoding='utf-8')
    print(json.dumps(records,indent=2))


if __name__=='__main__':main()
