"""Fetch the fixed GWTC-4 O1/O2/O3-only sensitivity release, never O4 data.

This control tests a newly documented semianalytic selection representation.
It is not silently substituted for the historical frozen input.
"""
from pathlib import Path
import hashlib
import json
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BASE=Path(__file__).resolve().parent
NAME='mixture-semi_o1_o2-real_o3-cartesian_spins_20250503134659UTC.hdf'


def main():
    session=requests.Session()
    session.mount('https://',HTTPAdapter(max_retries=Retry(total=3,backoff_factor=.5,status_forcelist=[429,502,503,504],respect_retry_after_header=False)))
    folder=BASE/'sources';meta_path=folder/'zenodo-16740128.json'
    if not meta_path.exists():
        response=session.get('https://zenodo.org/api/records/16740128',timeout=(10,35));response.raise_for_status()
        assert response.json()['id']==16740128;meta_path.write_bytes(response.content)
    meta=json.loads(meta_path.read_text());assert meta['id']==16740128
    entry=next(x for x in meta['files'] if x['key']==NAME);assert entry['size']==257339480
    path=folder/NAME
    if not path.exists():
        temporary=path.with_suffix('.partial');count=0
        with session.get(entry['links']['self'],stream=True,timeout=(10,45)) as response:
            response.raise_for_status()
            with temporary.open('wb') as out:
                for chunk in response.iter_content(1024*1024):
                    count+=len(chunk);assert count<=entry['size'];out.write(chunk)
                    if count%(32*1024*1024)==0:print(f'Downloaded {count}/{entry["size"]} bytes',flush=True)
        assert count==entry['size']
        assert 'md5:'+hashlib.md5(temporary.read_bytes()).hexdigest()==entry['checksum']
        temporary.replace(path)
    content=path.read_bytes();assert len(content)==entry['size'] and 'md5:'+hashlib.md5(content).hexdigest()==entry['checksum']
    record=dict(file=NAME,url=entry['links']['self'],bytes=len(content),official_md5=entry['checksum'],sha256=hashlib.sha256(content).hexdigest())
    manifest=folder/'manifest.json';rows=json.loads(manifest.read_text());rows=[r for r in rows if r['file']!=NAME]+[record]
    manifest.write_text(json.dumps(rows,indent=2),encoding='utf-8');print(json.dumps(record,indent=2),flush=True)


if __name__=='__main__':main()
