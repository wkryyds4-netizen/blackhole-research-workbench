"""Download pinned SXS strain products and verify the archive's MD5."""
from pathlib import Path
import concurrent.futures, hashlib, json, requests

ROOT = Path(__file__).resolve().parent
RECORD = json.loads((ROOT / 'adapter-sources/sxs-2327-record.json').read_text())

def fetch(key):
    entry = RECORD['files']['entries'][key]
    path = ROOT / 'waveform-data/2327' / key.replace(':', '/')
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        response = requests.get(entry['links']['content'], timeout=180)
        response.raise_for_status()
        path.write_bytes(response.content)
    data = path.read_bytes()
    assert len(data) == entry['size']
    assert hashlib.md5(data).hexdigest() == entry['checksum'].split(':')[1]
    return dict(key=key, bytes=len(data), url=entry['links']['content'],
                md5=hashlib.md5(data).hexdigest(), sha256=hashlib.sha256(data).hexdigest())

if __name__ == '__main__':
    import sys
    if '--0389' in sys.argv:
        record=json.loads((ROOT/'adapter-sources/sxs-0389-record.json').read_text())
        entries=[]
        for entry in record['files']:
            if entry['key']=='Lev3/rhOverM_Asymptotic_GeometricUnits_CoM.h5':
                path=ROOT/'waveform-data/0389'/entry['key']; path.parent.mkdir(parents=True,exist_ok=True)
                if not path.exists():
                    with requests.get(entry['links']['self'], timeout=180, stream=True) as response:
                        response.raise_for_status()
                        with path.open('wb') as f:
                            for block in response.iter_content(1024*1024): f.write(block)
                data=path.read_bytes()
                assert len(data)==entry['size']
                assert hashlib.md5(data).hexdigest()==entry['checksum'].split(':')[1]
                entries.append(dict(key=entry['key'],bytes=len(data),md5=hashlib.md5(data).hexdigest(),sha256=hashlib.sha256(data).hexdigest(),url=entry['links']['self']))
        out=ROOT/'wave-source-results';out.mkdir(exist_ok=True)
        (out/'0389-download-manifest.json').write_text(json.dumps(dict(doi='10.5281/zenodo.3274661',license='CC-BY-4.0',files=entries),indent=2))
        print(json.dumps(entries,indent=2));sys.exit(0)
    keys = [f'Lev{lev}:Strain_N2.{ext}' for lev in (2,3) for ext in ('h5','json')]
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        entries = list(pool.map(fetch, keys))
    out=ROOT/'wave-source-results'; out.mkdir(exist_ok=True)
    (out/'download-manifest.json').write_text(json.dumps(dict(
        doi='10.26138/SXS:BBH:2327v3.0', record='f70bw-gcd50',
        license='Record does not declare license; local research snapshot', files=entries),indent=2))
    print(json.dumps(entries, indent=2))
