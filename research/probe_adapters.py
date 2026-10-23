from pathlib import Path
import json,requests
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).resolve().parent/'adapter-sources';ROOT.mkdir(exist_ok=True)
TARGETS={
 'sxs-0001-record.json':'https://zenodo.org/api/records/2621619',
 'sxs-0389-record.json':'https://zenodo.org/api/records/3274661',
 'sxs-0003-record.json':'https://zenodo.org/api/records/13167015',
 'gwsurrogate-pypi.json':'https://pypi.org/pypi/gwsurrogate/json',
 'surrogate-catalog.py':'https://raw.githubusercontent.com/sxs-collaboration/gwsurrogate/13cd6a5a4667e0dcae626ddd11677a0d140b54eb/gwsurrogate/catalog.py',
 'NRSur7dq4.md':'https://raw.githubusercontent.com/sxs-collaboration/gwsurrogate/13cd6a5a4667e0dcae626ddd11677a0d140b54eb/docs/models/NRSur7dq4.md',
 'spec_horizons_h5.py':'https://raw.githubusercontent.com/sxs-collaboration/sxs/b7d25d5aadc2c15ffa4e989fe1c2f661040660b6/sxs/horizons/spec_horizons_h5.py',
 'precessing_surrogate.py':'https://raw.githubusercontent.com/sxs-collaboration/gwsurrogate/13cd6a5a4667e0dcae626ddd11677a0d140b54eb/gwsurrogate/new/precessing_surrogate.py',
}
def fetch(pair):
 name,url=pair
 try:
  r=requests.get(url,timeout=50);r.raise_for_status();(ROOT/name).write_bytes(r.content)
  if name.endswith('record.json'):
   d=r.json();print(name,json.dumps({'title':d['metadata']['title'],'license':d['metadata'].get('license'),'version':d['metadata'].get('version'),'files':[{k:f.get(k) for k in ['key','size','checksum','links']} for f in d['files'] if f['key'].endswith(('Horizons.h5','metadata.json'))]},ensure_ascii=False),flush=True)
  elif name.endswith('pypi.json'):
   d=r.json();print('PyPI',d['info']['version'],[(x['filename'],x['size']) for x in d['urls']],flush=True)
  else:print(name,len(r.content),flush=True)
 except Exception as e:print(name,repr(e),flush=True)
with ThreadPoolExecutor(max_workers=4) as p:list(p.map(fetch,TARGETS.items()))
