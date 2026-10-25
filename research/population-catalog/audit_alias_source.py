"""Connect HDF version metadata, historical conversion source, and controls."""
from pathlib import Path
import ast
import hashlib
import json
import h5py

BASE=Path(__file__).resolve().parent


def main():
    record=next(r for r in json.loads((BASE/'download-results.json').read_text()) if r['index']==14)
    path=BASE/'source-data'/record['source']['file'];assert hashlib.sha256(path.read_bytes()).hexdigest()==record['sha256']
    with h5py.File(path,'r') as h:
        version=[v.decode() for v in h['version/pesummary'][()]]
        canonical=h['C01:IMRPhenomXPHM'];mixed=h['C01:Mixed']
        c_ref=canonical['meta_data/meta_data/f_ref'][()].tolist()
        m_ref='f_ref' in mixed['meta_data/meta_data']
        ref_column={label:'reference_frequency' in h[label+'/posterior_samples'].dtype.names for label in ['C01:IMRPhenomXPHM','C01:Mixed']}
    commit=json.loads((BASE/'sources/pesummary-f4e1c37-commit.json').read_text())['sha']
    assert commit=='f4e1c37161f66419724706eaa5d81b0c9c5097a6' and any('.gf4e1c37' in v for v in version)
    source=BASE/'sources/pesummary-f4e1c37-conversions.py';raw=source.read_bytes()
    manifest=json.loads((BASE/'sources/manifest.json').read_text());entry=next(v for v in manifest if v['file']==source.name)
    assert hashlib.sha256(raw).hexdigest()==entry['sha256']
    tree=ast.parse(raw.decode());methods={node.name:node for node in ast.walk(tree) if isinstance(node,ast.FunctionDef)}
    reference=methods['_reference_frequency'];text=ast.get_source_segment(raw.decode(),reference)
    assert 'extra_kwargs["f_ref"]' in text and 'self.append_data("reference_frequency", [20.] * nsamples)' in text
    evolution=methods['_evolve_spins'];etext=ast.get_source_segment(raw.decode(),evolution)
    assert 'samples[5], samples[6], samples[7][0]' in etext and 'evolve_limit="infinite_separation"' in etext
    queries=[]
    for name in ['pesummary-issues-mixed-fref.json','pesummary-issues-infinity-reference.json']:
        observed=next(v for v in manifest if v['file']==name);queries.append(dict(file=name,total_count=observed['total_count'],scope='Recorded search observation in the provenance manifest; not a new live search.'))
    report=dict(schema='blackhole-historical-label-reference-source/1',source_sha256=record['sha256'],pesummary_version=version,historical_commit=commit,
        canonical_reference_hz=c_ref,mixed_metadata_has_f_ref=m_ref,posterior_reference_frequency_column_present=ref_column,
        historical_conversion_source=dict(sha256=entry['sha256'],reference_method_line=reference.lineno,evolution_method_line=evolution.lineno,
            url=f'https://github.com/pesummary/pesummary/blob/{commit}/pesummary/gw/conversions/__init__.py'),
        source_behavior='Absent f_ref metadata and absent reference_frequency column lead to a 20 Hz default; backward evolution receives that reference.',
        targeted_prior_issue_searches=queries,
        causal_status='Version metadata, native posterior equality, missing Mixed f_ref, historical fallback source, and numerical 16/20 Hz controls support a reference-default explanation. Actual historical invocation/logs are not recovered; a small approximately 1.2e-6 release residual remains on one held-out row.',
        novelty='Not established. Zero hits in two targeted issue searches do not establish absence of prior reports. Common-reference population inference is existing literature.')
    (BASE/'alias-source-audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
