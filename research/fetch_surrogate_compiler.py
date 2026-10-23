import hashlib,json,zipfile
from pathlib import Path
import requests
root=Path(__file__).resolve().parent/'adapter-build'
index=requests.get('https://ziglang.org/download/index.json',timeout=60).json()
version='0.14.1'
entry=index[version]['x86_64-windows']
archive=root/Path(entry['tarball']).name
if not archive.exists():
    with requests.get(entry['tarball'],stream=True,timeout=90) as r:
        r.raise_for_status()
        with archive.open('wb') as f:
            for chunk in r.iter_content(1024*1024): f.write(chunk)
assert hashlib.sha256(archive.read_bytes()).hexdigest()==entry['shasum']
with zipfile.ZipFile(archive) as z: z.extractall(root)
(root/'zig-manifest.json').write_text(json.dumps({'version':version,**entry},indent=2))
print('Verified and extracted',archive)
