"""Summarize complete metadata availability, preserving unverified source identity."""
from pathlib import Path
import collections
import hashlib
import json
import re

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]


def summarize_prior(values):
    if not isinstance(values,list) or len(values)!=1 or not isinstance(values[0],str):return dict(available=False)
    text=values[0]
    return dict(available=True,prior_class=text.split('(')[0],
        bounds={key:float(value) for key,value in re.findall(r'(minimum|maximum)=([\d.eE+\-]+)',text)})


def main():
    plan=json.loads((BASE/'source-plan.json').read_text());raw=json.loads((BASE/'metadata-results.json').read_text())
    metadata={r['index']:r for r in raw if r['mode']=='remote' and r.get('inspection_schema',1)>=2}
    expected={e['index'] for e in plan['events']};assert set(metadata)==expected,sorted(expected-set(metadata))
    verified={}
    for row in json.loads((ROOT/'research/population-stage7/event-provenance.json').read_text()):
        if row['all_tutorial_rows_exact_in_both']:verified[row['index']]=dict(label_unique=True,evidence='research/population-stage7/event-provenance.json')
    old=json.loads((ROOT/'research/population-stage5/representation-match.json').read_text())
    pilot=json.loads((ROOT/'research/population-stage3/pe-pilot-match.json').read_text())
    matched38=next(v for v in pilot['events'] if v['original_index']==38)
    if all(old['exact_matches'].values()) and any(v['full_match'] for v in matched38['labels']):verified[38]=dict(label_unique=True,evidence='research/population-stage3/pe-pilot-match.json + research/population-stage5/representation-match.json')
    for row in json.loads((BASE/'event-provenance.json').read_text()):
        if row.get('native_source_identity_resolved',False):verified[row['index']]=dict(label_unique=row['uniquely_verified_tutorial_source'],evidence='research/population-catalog/event-provenance.json')
    frequencies=collections.Counter();classes=collections.defaultdict(collections.Counter);rows=[]
    for pos,e in enumerate(plan['events']):
        m=metadata[e['index']];assert m['source_file']==e['source']['file'] and m['source_md5']==e['source']['md5']
        analyses=[v for v in m['labels'] if v['label']=='C01:IMRPhenomXPHM'];assert len(analyses)==1,e['event']
        a=analyses[0];freq=a['metadata'].get('meta_data/meta_data/f_ref');config=a['metadata'].get('config_file/config/reference-frequency')
        assert isinstance(freq,list) and len(freq)==1 and float(freq[0])>0,e['event']
        fref=float(freq[0]);frequencies[str(fref)]+=1
        if config:assert len(config)==1 and float(config[0])==fref
        prior={k:summarize_prior(v) for k,v in a.get('analytic_prior_strings',{}).items()}
        for k,v in prior.items():
            if v['available']:classes[k][v['prior_class']]+=1
        required=['a_1','a_2','cos_tilt_1','cos_tilt_2','phi_12','spin_1x','spin_1y','spin_1z','spin_2x','spin_2y','spin_2z']
        rows.append(dict(position=pos,index=e['index'],event=e['event'],source_file=m['source_file'],analysis_label=a['label'],
            native_reference_hz=fref,config_reference_present=bool(config),rows=a['shape'][0],
            native_spin_fields_complete=all(k in a['spin_fields_present'] for k in required),
            infinity_tilts_present=all(k in a['spin_fields_present'] for k in ['cos_tilt_1_infinity_only_prec_avg','cos_tilt_2_infinity_only_prec_avg']),
            prior_declarations=prior,full_core_prior_strings_present=all(k in prior and prior[k]['available'] for k in ['chirp_mass','mass_ratio','luminosity_distance','a_1','a_2','tilt_1','tilt_2']),
            tutorial_label_and_rows_independently_verified=e['index'] in verified and verified[e['index']]['label_unique'],
            native_posterior_source_identity_verified=e['index'] in verified,source_identity_evidence=verified.get(e['index'])))
    counts=dict(events=len(rows),native_spin_fields_complete=sum(r['native_spin_fields_complete'] for r in rows),
        infinity_tilts_present=sum(r['infinity_tilts_present'] for r in rows),
        core_prior_strings_present=sum(r['full_core_prior_strings_present'] for r in rows),
        config_reference_present=sum(r['config_reference_present'] for r in rows),reference_hz=dict(sorted(frequencies.items())),
        native_posterior_source_identity_verified=sum(r['native_posterior_source_identity_verified'] for r in rows))
    report=dict(schema='blackhole-full-catalog-metadata/1',counts=counts,
        prior_classes={k:dict(v) for k,v in classes.items()},events=rows,
        metadata_results_sha256=hashlib.sha256((BASE/'metadata-results.json').read_bytes()).hexdigest(),
        source_plan_sha256=hashlib.sha256((BASE/'source-plan.json').read_bytes()).hexdigest(),
        successful_access_totals=dict(requests=sum(m['access']['requests'] for m in metadata.values()),network_bytes=sum(m['access']['network_bytes'] for m in metadata.values()),cached_bytes=sum(m['access'].get('cached_bytes',0) for m in metadata.values()),
            scope='Access for successful stored inspections, excluding discarded passes and failed attempts.'),
        all_candidate_metadata_present=True,all_tutorial_source_rows_verified=len(verified)==len(rows),common_reference_population_ready=False,
        limits='Complete candidate-file metadata inventory. File naming/GPS matching is not PE-label or row identity. Range inspections do not verify whole-file MD5. Presence of analytic prior strings does not establish omitted likelihood tails or historical configuration completeness. Injection reference remains unresolved.')
    (BASE/'catalog-summary.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(dict(counts=counts,prior_classes=report['prior_classes'],missing_prior_events=[r['event'] for r in rows if not r['full_core_prior_strings_present']]),indent=2))


if __name__=='__main__':main()
