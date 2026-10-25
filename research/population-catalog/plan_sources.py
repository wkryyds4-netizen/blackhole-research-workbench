"""Resolve fixed official full-catalog candidates without asserting PE row identity."""
from pathlib import Path
import hashlib
import json
import re
from astropy.time import Time

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]


def main():
    model_path=ROOT/'results/model.json';model=json.loads(model_path.read_text())
    paths=[ROOT/f'research/population-stage3/sources/zenodo-{n}.json' for n in [6513631,5546663]]
    records=[json.loads(p.read_text()) for p in paths];events=[]
    for e in model['events']:
        candidates=[]
        for record in records:
            for f in record['files']:
                name=f['key']
                if not (e['name'] in name and name.endswith('_nocosmo.h5')):continue
                match=re.search(r'GW(\d{2})(\d{2})(\d{2})_(\d{2})(\d{2})(\d{2})_',name)
                assert match,name
                yy,mm,dd,hh,minute,ss=match.groups()
                gps=float(Time(f'20{yy}-{mm}-{dd}T{hh}:{minute}:{ss}',format='isot',scale='utc').gps)
                candidates.append(dict(record=record['id'],file=name,bytes=f['size'],md5=f['checksum'],
                    url=f['links']['self'],filename_gps=gps,filename_to_catalog_gps_abs_seconds=abs(gps-e['gps'])))
        candidates.sort(key=lambda row:row['filename_to_catalog_gps_abs_seconds'])
        assert candidates and candidates[0]['filename_to_catalog_gps_abs_seconds']<1.5,e['name']
        if len(candidates)>1:assert candidates[1]['filename_to_catalog_gps_abs_seconds']>60
        events.append(dict(index=e['index'],event=e['name'],source=candidates[0],
            status='Filename/time-resolved candidate; tutorial PE label and exact row identity must be checked separately.',
            rejected_alternatives=candidates[1:]))
    assert len(events)==69 and len({e['source']['file'] for e in events})==69
    plan=dict(schema='blackhole-full-catalog-source-plan/1',model_sha256=hashlib.sha256(model_path.read_bytes()).hexdigest(),
        metadata=[dict(path=str(p.relative_to(ROOT)).replace('\\','/'),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths],
        events=events,total_full_file_bytes=sum(e['source']['bytes'] for e in events),
        scope='Only a fixed source plan. Metadata presence, prior support and tutorial source identity are not established by filenames.')
    (BASE/'source-plan.json').write_text(json.dumps(plan,indent=2),encoding='utf-8')
    print(json.dumps(dict(events=len(events),total_full_file_bytes=plan['total_full_file_bytes'],ambiguous_names_resolved=[dict(event=e['event'],chosen=e['source']['file']) for e in events if e['rejected_alternatives']]),indent=2))


if __name__=='__main__':main()
