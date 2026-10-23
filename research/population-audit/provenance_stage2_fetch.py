"""Fetch explicitly enumerated small provenance metadata, never execute source."""
from pathlib import Path
import hashlib, json
import requests
OUT=Path(__file__).resolve().parents[1]/'population-stage2'/'provenance-sources'
OUT.mkdir(parents=True,exist_ok=True)
URLS={
 'notebook-history.json':'https://api.github.com/repos/ColmTalbot/gwpopulation/commits?path=examples/GWTC3.ipynb&per_page=100',
 'notebook-current.ipynb':'https://raw.githubusercontent.com/ColmTalbot/gwpopulation/main/examples/GWTC3.ipynb',
 'zenodo-7890398.json':'https://zenodo.org/api/records/7890398',
 'zenodo-5636816.json':'https://zenodo.org/api/records/5636816',
 'notebook-original.ipynb':'https://raw.githubusercontent.com/ColmTalbot/gwpopulation/2daa71e8c56362f2f64e2166a0055478557f239a/examples/GWTC3.ipynb',
 'repository-tree-2024.json':'https://api.github.com/repos/ColmTalbot/gwpopulation/git/trees/2daa71e8c56362f2f64e2166a0055478557f239a?recursive=1',
}
manifest=[]
for name,url in URLS.items():
    try:
        r=requests.get(url,timeout=40);r.raise_for_status()
        if len(r.content)>5_000_000:raise ValueError('Metadata exceeds size bound')
        (OUT/name).write_bytes(r.content)
        manifest.append(dict(file=name,url=url,final_url=r.url,bytes=len(r.content),sha256=hashlib.sha256(r.content).hexdigest()))
        print(name,len(r.content),flush=True)
    except Exception as e:
        manifest.append(dict(file=name,url=url,error=type(e).__name__))
        print(name,type(e).__name__,flush=True)
(OUT/'fetch-manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
