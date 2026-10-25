"""Bounded public GitLab API snapshots, pinned to an explicitly recorded commit."""
from pathlib import Path
import argparse
import hashlib
import json
from urllib.parse import quote
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BASE=Path(__file__).resolve().parent
FOLDER=BASE/'sources'


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('project',type=int,choices=[4440,9769,8759,1438])
    parser.add_argument('--ref',default='master')
    parser.add_argument('--before',help='Latest commit no later than this ISO timestamp')
    parser.add_argument('--files',nargs='*')
    args=parser.parse_args()
    session=requests.Session()
    session.mount('https://',HTTPAdapter(max_retries=Retry(total=2,backoff_factor=.4,status_forcelist=[429,502,503,504],respect_retry_after_header=False)))
    api=f'https://git.ligo.org/api/v4/projects/{args.project}'
    manifest_path=FOLDER/'manifest.json';manifest=json.loads(manifest_path.read_text())
    def fetch(url,name,params=None):
        with session.get(url,params=params,stream=True,timeout=(10,30)) as response:
            response.raise_for_status();chunks=[];n=0
            for block in response.iter_content(65536):
                n+=len(block)
                if n>2000000:raise ValueError('Source response exceeds 2 MB bound')
                chunks.append(block)
            content=b''.join(chunks)
            (FOLDER/name).write_bytes(content)
            record=dict(file=name,url=response.url,bytes=n,sha256=hashlib.sha256(content).hexdigest())
            manifest[:]=[v for v in manifest if v['file']!=name]+[record]
            manifest_path.write_text(json.dumps(manifest,indent=2),encoding='utf-8')
            return content,response.headers
    prefix=f'gitlab-{args.project}'
    if args.files:
        assert len(args.ref)==40 and all(c in '0123456789abcdef' for c in args.ref)
        for source in args.files:
            assert not source.startswith('/') and '..' not in Path(source).parts
            name=f'{prefix}-{args.ref[:12]}-'+source.replace('/','--')
            fetch(api+'/repository/files/'+quote(source,safe='')+'/raw',name,{'ref':args.ref})
            print(name,flush=True)
    else:
        if args.before:
            stamp=args.before.replace(':','').replace('-','')
            candidates=json.loads(fetch(api+'/repository/commits',prefix+'-before-'+stamp+'.json',dict(ref_name=args.ref,until=args.before,per_page=1))[0])
            assert len(candidates)==1
            commit=candidates[0]
        else:
            commit=json.loads(fetch(api+'/repository/commits/'+quote(args.ref,safe=''),prefix+'-commit.json')[0])
        ref=commit['id'];rows=[]
        for page in range(1,5):
            content,headers=fetch(api+'/repository/tree',f'{prefix}-{ref[:12]}-tree-{page}.json',dict(ref=ref,recursive='true',per_page=100,page=page))
            rows+=json.loads(content)
            if not headers.get('X-Next-Page'):break
        else:raise ValueError('Repository tree exceeds four-page bound')
        print(json.dumps(dict(project=args.project,commit=ref,committed_date=commit['committed_date'],paths=[r['path'] for r in rows if r['type']=='blob']),indent=2))


if __name__=='__main__':main()
