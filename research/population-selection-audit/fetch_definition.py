"""Fetch the fixed official injection definition, validating its release checksum."""
from pathlib import Path
import hashlib
import json
import requests

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]


def main():
    meta=json.loads((ROOT/'research/population-stage2/provenance-sources/zenodo-7890398.json').read_text())
    entry=next(x for x in meta['files'] if x['key']=='o1+o2+o3-sensitivity-estimates.md')
    folder=BASE/'sources';folder.mkdir(parents=True,exist_ok=True)
    path=folder/entry['key']
    if not path.exists():
        r=requests.get(entry['links']['self'],timeout=(15,40));r.raise_for_status()
        assert len(r.content)==entry['size']
        assert 'md5:'+hashlib.md5(r.content).hexdigest()==entry['checksum']
        path.write_bytes(r.content)
    content=path.read_bytes()
    assert len(content)==entry['size'] and 'md5:'+hashlib.md5(content).hexdigest()==entry['checksum']
    record=dict(file=entry['key'],url=entry['links']['self'],bytes=len(content),
                official_md5=entry['checksum'],sha256=hashlib.sha256(content).hexdigest())
    (folder/'manifest.json').write_text(json.dumps([record],indent=2),encoding='utf-8')
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
