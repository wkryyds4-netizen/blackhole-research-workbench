"""Download two bounded PE pilots, verifying published size and MD5."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import time
import requests

BASE = Path(__file__).resolve().parent
PILOTS = [(6513631, 'GW190805_211137'), (5546663, 'GW191103_012549')]

def download(record, event):
    meta = json.loads((BASE / f'sources/zenodo-{record}.json').read_text())
    candidates = [f for f in meta['files'] if event in f['key'] and f['key'].endswith('_mixed_cosmo.h5')]
    assert len(candidates) == 1
    entry = candidates[0]
    assert entry['size'] < 60_000_000
    path = BASE / 'source-data' / entry['key']
    url = entry['links']['self']
    if not path.exists():
        def fetch(start):
            end = min(start + 2 * 1024 * 1024, entry['size']) - 1
            for attempt in range(3):
                try:
                    r = requests.get(url, headers={'Range': f'bytes={start}-{end}'}, timeout=(15, 40))
                    r.raise_for_status()
                    assert r.status_code == 206
                    assert r.headers.get('Content-Range') == f'bytes {start}-{end}/{entry["size"]}'
                    assert len(r.content) == end - start + 1
                    return start, r.content
                except (requests.RequestException, AssertionError):
                    if attempt == 2:
                        raise
                    time.sleep(1)
        path.parent.mkdir(exist_ok=True)
        with path.with_suffix('.part').open('wb') as output, ThreadPoolExecutor(max_workers=4) as pool:
            for start, data in pool.map(fetch, range(0, entry['size'], 2 * 1024 * 1024)):
                output.seek(start)
                output.write(data)
        data = path.with_suffix('.part').read_bytes()
        assert len(data) == entry['size']
        assert 'md5:' + hashlib.md5(data).hexdigest() == entry['checksum']
        path.with_suffix('.part').replace(path)
    data = path.read_bytes()
    assert len(data) == entry['size']
    assert 'md5:' + hashlib.md5(data).hexdigest() == entry['checksum']
    result = dict(record=record, event=event, file=path.name, url=url, bytes=len(data),
                  checksum=entry['checksum'], sha256=hashlib.sha256(data).hexdigest())
    print(json.dumps(result), flush=True)
    return result

if __name__ == '__main__':
    result = [download(record, event) for record, event in PILOTS]
    (BASE / 'pilot-downloads.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
