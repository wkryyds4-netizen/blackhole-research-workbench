"""Read public repository metadata and small source files; never execute them."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import requests
import json

ROOT=Path(__file__).resolve().parent/'snapshots'
REPOS=['vijayvarma392/binaryBHexp','MdKasif0/Gravitational-Wave-Simulation','algometrix/blackhole-sim','cardiffgravity/waveform-fitter','gwastro/pycbc','sxs-collaboration/gwsurrogate','dalyagergely/SC_interactive','TianQinSYSU/GWSpace']

def fetch(repo):
    out=ROOT/repo.replace('/','__');out.mkdir(parents=True,exist_ok=True)
    s=requests.Session();s.headers['User-Agent']='GW-visualization-research'
    def get(url):
        r=s.get(url,timeout=35);r.raise_for_status();return r
    try:
        meta=get('https://api.github.com/repos/'+repo).json()
        branch=meta['default_branch']
        commit=get(f'https://api.github.com/repos/{repo}/commits/{branch}').json()
        sha=commit['sha']
        tree=get(f'https://api.github.com/repos/{repo}/git/trees/{sha}?recursive=1').json()
        paths=[e['path'] for e in tree.get('tree',[]) if e['type']=='blob']
        record={k:meta.get(k) for k in ['full_name','html_url','homepage','default_branch','pushed_at','archived','license']}
        record.update(commit=sha,commit_date=commit['commit']['committer']['date'],retrieved='2026-10-23',paths=paths)
        (out/'metadata.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
        wanted=[p for p in paths if ('/' not in p and (p.lower().startswith(('readme','license','copying')) or p in ['package.json','setup.py','pyproject.toml']))]
        if repo.startswith('MdKasif0/'):
            wanted += [p for p in paths if p.endswith(('.ts','.tsx','.glsl')) and any(x in p.lower() for x in ['physic','simulat','waveform','audio','store'])]
        if repo.startswith('algometrix/'):
            wanted += [p for p in paths if p.endswith(('.ts','.md')) and any(x in p.lower() for x in ['binary','theory','inspiral','physics.test'])]
        if repo.startswith('cardiffgravity/'):
            wanted += [p for p in paths if p.endswith('.js') and not any(x in p for x in ['min.js','jquery','bootstrap','d3.'])]
        if repo.startswith('vijayvarma392/'):
            wanted += ['binaryBHexp.py','interactive.py','docs/index.md']
        saved=[]
        for p in dict.fromkeys(wanted):
            if p not in paths:continue
            file=get(f'https://raw.githubusercontent.com/{repo}/{sha}/{p}')
            if len(file.content)>500000:continue
            target=out/p;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(file.content);saved.append(p)
        print(json.dumps({'repo':repo,'sha':sha[:10],'date':record['commit_date'],'license':(meta.get('license') or {}).get('spdx_id'),'saved':saved},ensure_ascii=False),flush=True)
    except Exception as exc:print(repo,type(exc).__name__,str(exc)[:160],flush=True)

if __name__=='__main__':
    with ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(fetch,REPOS))
