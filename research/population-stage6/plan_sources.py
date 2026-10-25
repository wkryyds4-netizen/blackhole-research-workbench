"""Resolve the predeclared follow-up ranking using already retrieved release metadata."""
from pathlib import Path
import json
import hashlib

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]


def main():
    screen=json.loads((BASE/'catalog-screen.json').read_text())
    manifest=json.loads((ROOT/'research/population-stage3/sources/manifest.json').read_text())
    metadata={};sources=[]
    for record in [6513631,5546663]:
        path=ROOT/f'research/population-stage3/sources/zenodo-{record}.json'
        sha=hashlib.sha256(path.read_bytes()).hexdigest()
        assert sha==next(x['sha256'] for x in manifest if x['file']==path.name)
        metadata[record]=json.loads(path.read_text())
        sources.append(dict(record=record,file=path.name,sha256=sha))
    rows=[]
    for ix in screen['followup_indices']:
        event=next(e for e in screen['events'] if e['index']==ix)
        files=[]
        for record,meta in metadata.items():
            for file in meta['files']:
                if event['event'] in file['key'] and file['key'].endswith(('_mixed_cosmo.h5','_mixed_nocosmo.h5')):
                    files.append(dict(record=record,file=file['key'],bytes=file['size'],
                        md5=file['checksum'],url=file['links']['self']))
        assert len(files)==2
        rows.append(dict(index=ix,event=event['event'],identity_status=event['identity_status'],
            selection_reasons=[key for key,order in screen['rankings'].items() if ix in order[:3]],
            files=files,analysis_label_status='not yet inspected or matched for these new candidates'))
    result=dict(schema='blackhole-population-followup-source-plan/1',metadata_snapshots=sources,
        screen_sha256=hashlib.sha256((BASE/'catalog-screen.json').read_bytes()).hexdigest(),events=rows,
        total_bytes=sum(f['bytes'] for e in rows for f in e['files']),
        first_pair=rows[0]['event'],downloads_performed=False,
        caution='Filename and metadata match only. No confirmed tutorial row identities, prior support, or sampling counts for the new files.')
    (BASE/'followup-source-plan.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(dict(candidates=[e['event'] for e in rows],total_bytes=result['total_bytes'],first_pair=result['first_pair']),indent=2))


if __name__=='__main__':main()
