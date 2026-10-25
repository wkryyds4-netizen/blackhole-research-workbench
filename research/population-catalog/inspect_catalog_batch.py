"""Process-isolated HDF metadata readers; only the parent commits checkpoints."""
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
import argparse
import json
from inspect_catalog import BASE,inspect,h5py,RangeReader


def worker(e):
    s=e['source'];reader=RangeReader(s['url'],s['bytes'],cache_key=s['md5'],block_size=1024*1024,min_interval=3)
    try:
        with h5py.File(reader,'r') as h:labels=inspect(h)
        return dict(inspection_schema=2,index=e['index'],event=e['event'],source_file=s['file'],source_md5=s['md5'],mode='remote',labels=labels,
            access=dict(requests=reader.requests,network_bytes=reader.network_bytes,cached_bytes=reader.cached_bytes,
                object_identity_header=reader.identity,ranges=reader.ranges,whole_file_md5_verified=False,
                limits='Exact Content-Range/size/header checks and per-block hashes; no full-file checksum for remote metadata reads.'))
    finally:reader.close()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=2);args=parser.parse_args();assert 1<=args.workers<=2
    plan=json.loads((BASE/'source-plan.json').read_text());output=BASE/'metadata-results.json'
    rows=json.loads(output.read_text()) if output.exists() else []
    done={r['index'] for r in rows if r['mode']=='remote' and r.get('inspection_schema',1)>=2}
    todo=[e for e in plan['events'] if e['index'] not in done];errors=[];iterator=iter(todo);limited=False
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures={}
        def submit():
            e=next(iterator,None)
            if e is not None:futures[pool.submit(worker,e)]=e
        for _ in range(args.workers):submit()
        while futures:
            ready,_=wait(futures,return_when=FIRST_COMPLETED)
            for future in ready:
                e=futures.pop(future)
                try:
                    row=future.result()
                    rows=[v for v in rows if not (v['index']==row['index'] and v['mode']=='remote')]+[row]
                    temporary=output.with_suffix('.tmp');temporary.write_text(json.dumps(rows,indent=2),encoding='utf-8');temporary.replace(output)
                    done.add(row['index'])
                    print(json.dumps(dict(completed=len(done),index=row['index'],event=row['event'],network_bytes=row['access']['network_bytes'],cached_bytes=row['access']['cached_bytes'])),flush=True)
                except Exception as error:
                    failure=dict(index=e['index'],event=e['event'],error=type(error).__name__+': '+str(error));errors.append(failure)
                    (BASE/'metadata-errors.json').write_text(json.dumps(errors,indent=2),encoding='utf-8');print(json.dumps(failure),flush=True)
                    limited |= '429' in str(error)
                if not limited:submit()
    print(json.dumps(dict(completed=len(done),failed=len(errors),missing=sorted({e['index'] for e in plan['events']}-done))),flush=True)
    if errors:raise SystemExit(1)


if __name__=='__main__':main()
