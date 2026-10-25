"""Restore remaining full-catalog files, reusing verified local originals."""
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
from pathlib import Path
import argparse
import json
import shutil
from download_sources import BASE,ROOT,download_one


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--remaining',action='store_true',required=True);parser.add_argument('--workers',type=int,default=2)
    args=parser.parse_args();assert 1<=args.workers<=2
    plan=json.loads((BASE/'source-plan.json').read_text());output=BASE/'download-results.json';rows=json.loads(output.read_text()) if output.exists() else []
    done={r['index'] for r in rows if r['whole_file_md5_verified'] and (ROOT/r.get('local_file','research/population-catalog/source-data/'+r['source']['file'])).exists()}
    todo=[e for e in plan['events'] if e['index'] not in done]
    new_bytes=sum(e['source']['bytes'] for e in todo if not any((folder/e['source']['file']).exists() for folder in [BASE/'source-data',ROOT/'research/population-stage7/source-data',ROOT/'research/population-stage5/source-data']))
    assert shutil.disk_usage(BASE).free>new_bytes+5*1024**3,'Insufficient free disk for the fixed bounded plan'
    print(json.dumps(dict(already_verified=len(done),pending=len(todo),maximum_new_bytes=new_bytes)),flush=True)
    iterator=iter(todo);errors=[];limited=False
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures={}
        def submit():
            e=next(iterator,None)
            if e is not None:futures[pool.submit(download_one,e)]=e
        for _ in range(args.workers):submit()
        while futures:
            ready,_=wait(futures,return_when=FIRST_COMPLETED)
            for future in ready:
                event=futures.pop(future)
                try:
                    row=future.result();rows=[r for r in rows if r['index']!=row['index']]+[row]
                    temp=output.with_suffix('.tmp');temp.write_text(json.dumps(rows,indent=2),encoding='utf-8');temp.replace(output);done.add(row['index'])
                    print(json.dumps(dict(verified=len(done),index=row['index'],sha256=row['sha256'],reused=row['reused_local_file'])),flush=True)
                except Exception as error:
                    failure=dict(index=event['index'],error=type(error).__name__+': '+str(error));errors.append(failure);print(json.dumps(failure),flush=True)
                    limited |= '429' in str(error)
                if not limited:submit()
    report=dict(verified=len(done),planned=69,missing=sorted({e['index'] for e in plan['events']}-done),errors=errors)
    (BASE/'download-batch-status.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report),flush=True)
    if report['missing']:raise SystemExit(1)


if __name__=='__main__':main()
