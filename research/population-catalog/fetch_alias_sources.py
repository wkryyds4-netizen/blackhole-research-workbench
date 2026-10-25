"""Restore fixed historical conversion sources using the recorded checksums."""
from pathlib import Path
import hashlib
import json
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BASE=Path(__file__).resolve().parent


def main():
    folder=BASE/'sources';manifest=json.loads((folder/'manifest.json').read_text())
    session=requests.Session();session.mount('https://',HTTPAdapter(max_retries=Retry(total=2,backoff_factor=.5,status_forcelist=[429,502,503,504],respect_retry_after_header=True)))
    for name in ['pesummary-f4e1c37-commit.json','pesummary-f4e1c37-conversions.py']:
        record=next(r for r in manifest if r['file']==name);path=folder/name
        if not path.exists():
            with session.get(record['url'],stream=True,timeout=(10,30)) as response:
                response.raise_for_status();chunks=[];count=0
                for block in response.iter_content(65536):
                    count+=len(block);assert count<=2000000;chunks.append(block)
                data=b''.join(chunks)
            assert len(data)==record['bytes'] and hashlib.sha256(data).hexdigest()==record['sha256']
            path.write_bytes(data)
        assert hashlib.sha256(path.read_bytes()).hexdigest()==record['sha256']
        print('Verified '+name,flush=True)


if __name__=='__main__':main()
