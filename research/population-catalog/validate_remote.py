"""Validate range-reader output and accessed bytes against a verified local HDF."""
from pathlib import Path
import hashlib
import json

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]


def main():
    rows=json.loads((BASE/'metadata-results.json').read_text());checks=[]
    for index in sorted({r['index'] for r in rows}):
        local=next((r for r in rows if r['index']==index and r['mode']=='local'),None)
        remote=next((r for r in rows if r['index']==index and r['mode']=='remote'),None)
        if local is None or remote is None:continue
        assert local['source_md5']==remote['source_md5'] and local['labels']==remote['labels']
        path=ROOT/local['access']['local_path']
        assert 'md5:'+hashlib.md5(path.read_bytes()).hexdigest()==remote['source_md5']
        total=0
        with path.open('rb') as stream:
            for block in remote['access']['ranges']:
                stream.seek(block['start']);data=stream.read(block['end']-block['start']+1)
                assert hashlib.sha256(data).hexdigest()==block['sha256'];total+=len(data)
        assert total<=remote['access']['network_bytes']+remote['access'].get('cached_bytes',0)
        checks.append(dict(index=index,all_metadata_equal=True,ranges_exact=len(remote['access']['ranges']),bytes_exact=total,whole_local_md5_verified=True))
    assert checks,'Need at least one local-versus-remote overlap'
    report=dict(schema='blackhole-remote-hdf-validation/1',checks=checks,
        scope='Exact requested-byte and parsed-metadata agreement on overlapping local verified files only. Other remote files do not gain a full-file checksum certification.')
    (BASE/'remote-validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
