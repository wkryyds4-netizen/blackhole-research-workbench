"""Fetch only pinned public SXS horizon files; verify Zenodo checksums."""
from pathlib import Path
import hashlib, json, requests
from concurrent.futures import ThreadPoolExecutor
ROOT = Path(__file__).resolve().parent
TARGETS = [('0001', 'Lev5'), ('0389', 'Lev3'), ('0389', 'Lev2'), ('0003', 'Lev5')]

def fetch(pair):
    simulation, level = pair
    record = json.loads((ROOT/'adapter-sources'/f'sxs-{simulation}-record.json').read_text())
    out = ROOT/'sxs-data'/simulation/level
    out.mkdir(parents=True, exist_ok=True)
    manifest = {'simulation': f'SXS:BBH:{simulation}', 'record_id': record['id'],
                'doi': record['doi'], 'license': record['metadata']['license'], 'level': level, 'files': []}
    for f in record['files']:
        key = f['key'].replace(':', '/')
        if f'/{level}/' not in '/'+key or not key.endswith(('Horizons.h5','metadata.json')):
            continue
        url = f['links']['self']; dest = out/key.split('/')[-1]
        if not dest.exists():
            r = requests.get(url, timeout=120); r.raise_for_status(); dest.write_bytes(r.content)
        data = dest.read_bytes()
        assert len(data) == f['size']
        assert 'md5:'+hashlib.md5(data).hexdigest() == f['checksum'], dest
        manifest['files'].append({'file':dest.name, 'url':url, 'bytes':len(data),
            'md5':f['checksum'], 'sha256':hashlib.sha256(data).hexdigest()})
    assert len(manifest['files']) == 2
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps(manifest), flush=True)

if __name__ == '__main__':
    with ThreadPoolExecutor(max_workers=3) as pool: list(pool.map(fetch, TARGETS))
