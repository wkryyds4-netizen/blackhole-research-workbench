"""Bounded read-only primary-source fetches; no upstream code execution."""
from pathlib import Path
import argparse
import hashlib
import json
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BASE=Path(__file__).resolve().parent
FOLDER=BASE/'sources'
SESSION=requests.Session()
SESSION.mount('https://',HTTPAdapter(max_retries=Retry(total=3,backoff_factor=.5,status_forcelist=[429,502,503,504],respect_retry_after_header=False)))


def fetch(url,name,limit=4000000):
    FOLDER.mkdir(exist_ok=True)
    with SESSION.get(url,stream=True,timeout=(15,40)) as response:
        chunks=[];count=0
        for block in response.iter_content(65536):
            count+=len(block)
            if count>limit:raise ValueError('Response exceeds recorded source-size bound')
            chunks.append(block)
        content=b''.join(chunks)
        record=dict(file=name,url=url,final_url=response.url,http_status=response.status_code,
                    bytes=len(content),sha256=hashlib.sha256(content).hexdigest())
    (FOLDER/name).write_bytes(content)
    manifest=FOLDER/'manifest.json'
    previous=json.loads(manifest.read_text()) if manifest.exists() else []
    manifest.write_text(json.dumps([x for x in previous if x['file']!=name]+[record],indent=2)+'\n',encoding='utf-8')
    return content,record


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--repo',default='michaelzevin/AMAZE')
    parser.add_argument('--record',type=int,default=4448170)
    parser.add_argument('--ref')
    parser.add_argument('--files',nargs='*')
    args=parser.parse_args()
    assert args.repo.count('/')==1 and all(c.isalnum() or c in '_-./' for c in args.repo)
    if args.files:
        assert args.ref and len(args.ref)==40 and all(c in '0123456789abcdef' for c in args.ref)
        for path in args.files:
            assert not path.startswith('/') and '..' not in Path(path).parts
            _,rec=fetch(f'https://raw.githubusercontent.com/{args.repo}/{args.ref}/{path}',
                        'github-'+args.repo.replace('/','--')+'-'+args.ref[:12]+'-'+path.replace('/','--'))
            print(json.dumps(rec),flush=True)
        return
    raw,rec=fetch(f'https://zenodo.org/api/records/{args.record}',f'zenodo-{args.record}.json')
    if rec['http_status']==200:
        d=json.loads(raw)
        print(json.dumps(dict(record=d['id'],title=d['metadata']['title'],license=d['metadata'].get('license'),
                              files=[{k:f[k] for k in ['key','size','checksum']} for f in d['files']]),indent=2),flush=True)
    else:print(json.dumps(rec),flush=True)
    raw,rec=fetch('https://api.github.com/repos/'+args.repo,'github-'+args.repo.replace('/','--')+'-repo.json')
    assert rec['http_status']==200
    repo=json.loads(raw);branch=repo['default_branch']
    raw,rec=fetch('https://api.github.com/repos/'+repo['full_name']+'/commits/'+branch,
                  'github-'+args.repo.replace('/','--')+'-commit.json')
    assert rec['http_status']==200
    ref=json.loads(raw)['sha']
    raw,rec=fetch('https://api.github.com/repos/'+repo['full_name']+'/git/trees/'+ref+'?recursive=1',
                  'github-'+args.repo.replace('/','--')+'-'+ref[:12]+'-tree.json')
    assert rec['http_status']==200
    tree=json.loads(raw);assert not tree.get('truncated')
    print(json.dumps(dict(repo=repo['full_name'],branch=branch,ref=ref,license=repo.get('license'),
                          files=[r['path'] for r in tree['tree'] if r['type']=='blob']),indent=2))


if __name__=='__main__':main()
