"""Reproduce a read-only provenance audit from safe JSON/NPZ and pickle opcodes.

No pickle object is deserialized. No inference input is changed.
"""
from pathlib import Path
import hashlib
import json
import pickletools
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
def read(path): return json.loads(path.read_text(encoding='utf-8'))
def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()

audit=read(ROOT/'data/sample_audit.json')
catalog_paths=sorted((ROOT/'data').glob('GWTC*.json'))
catalogs={p.stem:read(p)['events'] for p in catalog_paths}
all_by_name={key:{e['commonName']:e for e in vals.values()} for key,vals in catalogs.items()}
merged={}
for path in catalog_paths:
    for event in catalogs[path.stem].values():
        if event['mass_1_source'] is not None: merged[event['commonName']]=event
recovered=[]
for item in audit:
    if not item['included'] or item['far'] is not None: continue
    old=all_by_name['GWTC-1-confident'][item['name']]
    overwritten=merged[item['name']]
    assert overwritten['catalog.shortName']=='GWTC-2.1-confident'
    assert overwritten['far'] is None and old['far'] is not None
    assert old['far_unit']=='yr^-1'
    recovered.append(dict(index=item['index'],name=item['name'],
        posterior_identity_status='summary-matched; not original ID recovered',
        original_audit_far=None,
        recovered_far={k:old[k] for k in ['far','far_lower','far_upper','far_unit','version','GPS','jsonurl','catalog.shortName']},
        overwritten_record={k:overwritten[k] for k in ['far','far_unit','version','GPS','jsonurl','catalog.shortName']},
        cross_version=True,
        satisfies_far_less_than_1_per_year=bool(old['far']<1)))
assert len(recovered)==10
assert all(e['satisfies_far_less_than_1_per_year'] for e in recovered)

indices=list(range(len(audit)))
del indices[15]
del indices[38]
assert set(range(len(audit)))-set(indices)=={15,39}
extra=[dict(index=e['index'],name=e['name'],far=e['far'],
            far_unit=merged[e['name']]['far_unit'],source_url=merged[e['name']]['jsonurl'],
            reason=e['exclusion'],identity_status='summary-matched')
       for e in audit if e['index'] in indices and not e['included']]
assert len(extra)==2

string_ops={'SHORT_BINUNICODE','BINUNICODE','UNICODE','STRING','BINUNICODE8'}
scans={}
for name in ['gwtc-3-samples.pkl','gwtc-3-injections.pkl']:
    strings=sorted(set(arg for op,arg,pos in pickletools.genops((ROOT/'data'/name).read_bytes()) if op.name in string_ops))
    scans[name]=dict(method='pickletools.genops; no deserialization or callable execution',strings=strings,
                     explicit_GW_event_strings=[s for s in strings if s.startswith('GW')])
inj=np.load(ROOT/'data/injections.npz',allow_pickle=False)
live=OUT/'provenance-sources/GWTC-1-live.json'
source_files=['model.py','prepare_data.py','data/sample_audit.json','data/manifest.json',
              'data/posteriors.npz','data/injections.npz','data/gwtc-3-samples.pkl','data/gwtc-3-injections.pkl']
source_files += [str(p.relative_to(ROOT)).replace('\\','/') for p in catalog_paths]
output=dict(schema='blackhole-population-provenance-followup/1',date='2026-10-24',
    status='partial provenance repair; original inference unchanged; selection consistency not established',
    finding_far_overwrite=dict(cause='model.identify_events dictionary keyed by commonName; sorted catalogs overwrite GWTC-1 entries with GWTC-2.1 entries, including null FAR',
       file_order=[p.name for p in catalog_paths],missing_selected_count=10,recovered_catalog_far_count=10,
       selection_count_if_only_far_recovery_applied=69,
       live_GWTC1_byte_identical_to_local=live.exists() and digest(live)==digest(ROOT/'data/GWTC-1-confident.json'),
       caveat='GWOSC catalog FAR numeric values are preserved verbatim. Null upper/lower fields do not prove the numbers are exact measured rates rather than reported limits. Catalog FAR provenance is separate from PE release provenance.',events=recovered),
    tutorial_comparison=dict(source='https://colmtalbot.github.io/gwpopulation/examples/GWTC3.html',input_count=73,
       sequential_deletions=[15,38],original_deleted_indices=[15,39],inferred_deleted_names=[audit[i]['name'] for i in [15,39]],
       tutorial_remaining_count=len(indices),current_count=sum(e['included'] for e in audit),additional_current_exclusions=extra,
       caveat='Names remain conditional on local summary matching. The tutorial itself warns its event list is only qualitatively chosen.'),
    injection_provenance=dict(tutorial_points_to_release='https://zenodo.org/records/7890398',
       fields=list(inj.files),n_retained=len(inj['prior']),total_generated=int(inj['total_generated']),analysis_time=float(inj['analysis_time']),
       analysis_time_unit='years according to current pipeline convention; original preprocessing invocation absent',
       original_file_variant=None,original_pipeline_commit=None,original_ifar_threshold=None,original_snr_threshold=None,per_run_labels=None,
       retained_search_statistics=False,retained_original_row_ids=False,
       current_pipeline_reference='https://ratesandpopulations.docs.ligo.org/gwpopulation_pipe/_modules/gwpopulation_pipe/vt_helper.html',
       current_pipeline_defaults_only=dict(ifar_threshold_years=1,snr_threshold=10,
         selection='OR of IFAR above threshold, plus O1/O2 optimal network SNR above threshold',
         prior_transformation='sampling_pdf * mass_1_source * (2*pi*a_1**2) * (2*pi*a_2**2)',
         warning='Current documentation is not the historical tutorial invocation. Defaults do not establish the actual cuts or Jacobian applied to the local file.'),
       official_release_method='O1/O2 semi-analytic SNR; O3 real search statistics; thresholds must be specified separately',
       unresolved=['exact source HDF variant and checksum','original threshold and pipeline set','per-run detection masks','original Jacobian and mixture weighting configuration']),
    pe_provenance=dict(tutorial_release_links=['https://zenodo.org/records/5546663','https://zenodo.org/records/6513631'],
       release_options='GWTC-2.1 release provides cosmo and nocosmo; multiple waveform labels and Mixed are available',
       local_waveform_labels=None,local_prior_configs=None,local_event_identifiers=None,
       warning='Release options and pipeline formulas are reference documentation, not proof which options produced these tutorial arrays.'),
    pickle_opcode_inspection=scans,
    unresolved=['Original per-event IDs and PE waveform/prior metadata','Complete event-selection policy including mass classification','Matched injection selection and measures'],
    input_hashes={name:digest(ROOT/name) for name in source_files},
    source_snapshot_manifest='provenance-sources/manifest.json',
    changes_to_analysis='none')
(OUT/'provenance-followup.json').write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(dict(recovered=len(recovered),tutorial_count=len(indices),current_count=output['tutorial_comparison']['current_count'],extra=extra,live_catalog_identical=output['finding_far_overwrite']['live_GWTC1_byte_identical_to_local']),indent=2))
