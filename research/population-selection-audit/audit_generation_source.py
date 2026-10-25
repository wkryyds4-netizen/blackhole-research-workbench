"""Read-only provenance checks of pinned upstream implementation snapshots."""
from pathlib import Path
import ast
import hashlib
import json

BASE=Path(__file__).resolve().parent
SOURCES=BASE/'sources'
REF='ebbc51a29a9cf1785c0f1bbdf6c3cef2e2cc3bfd'


def main():
    manifest={r['file']:r for r in json.loads((SOURCES/'manifest.json').read_text())}
    tree={r['path']:r for p in SOURCES.glob(f'gitlab-8759-{REF[:12]}-tree-*.json') for r in json.loads(p.read_text())}
    verified=[]
    def source(path):
        name=f'gitlab-8759-{REF[:12]}-'+path.replace('/','--');raw=(SOURCES/name).read_bytes()
        sha=hashlib.sha256(raw).hexdigest();assert sha==manifest[name]['sha256']
        blob=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest();assert blob==tree[path]['id']
        text=raw.decode();module=ast.parse(text)
        verified.append(dict(path=path,sha256=sha,git_blob_sha1=blob,
            url=f'https://git.ligo.org/reed.essick/gw-distributions/-/blob/{REF}/{path}'))
        return text,module
    detection,dtree=source('gwdistributions/transforms/detection/detection.py')
    waveform,wtree=source('gwdistributions/transforms/detection/waveform.py')
    generator,gtree=source('gwdistributions/generators.py')
    def method(module,cls,name):
        c=next(v for v in module.body if isinstance(v,ast.ClassDef) and v.name==cls)
        return next(v for v in c.body if isinstance(v,(ast.FunctionDef,ast.AsyncFunctionDef)) and v.name==name)
    def contains(text,node,needle):
        assert needle in ast.get_source_segment(text,node)
        return node.lineno
    observed=method(dtree,'ObservedSNRGivenOptimalSNR','transform')
    contains(detection,observed,'real, imag = be.normal(size=2)')
    contains(detection,observed,'obs_phs2 = (obs_snr2 + imag**2)')
    contains(detection,observed,'obs_phs_net += obs_phs2')
    assert 'num_degrees_of_freedom=2*num_detectors' in detection
    frequency=method(wtree,'Waveform','transform')
    contains(waveform,frequency,'if fref is None:')
    contains(waveform,frequency,'fref = wave_flow')
    # Identify the generator class without guessing its public class name.
    gc=next(v for v in gtree.body if isinstance(v,ast.ClassDef) and any(isinstance(w,ast.FunctionDef) and w.name=='logprob' for w in v.body))
    logprob=method(gtree,gc.name,'logprob');logtext=ast.get_source_segment(generator,logprob)
    assert 'self.distributions' in logtext and 'self.conditionals' not in logtext
    assert 'self._total_generated += 1' in generator and 'for func in self.conditionals:' in generator and 'if check:' in generator
    report=dict(schema='blackhole-generation-source-provenance/1',project=8759,commit=REF,
        cutoff='2025-05-03T13:46:59Z',commit_selection='Latest public master commit no later than the cumulative filename timestamp, not a recovered execution commit.',
        verified_sources=verified,
        findings=dict(independent_gaussian_quadratures_implementation=dict(confirmed=True,line=observed.lineno),
            phase_maximized_network_dof='2*N_IFO; source agrees with the declared ncx2 noise model',
            unspecified_reference_frequency=dict(behavior='If unset, fref defaults to waveform low frequency; that low frequency may depend on total detector-frame mass.',line=frequency.lineno),
            stored_generator_logprob=dict(behavior='Sum of base-distribution log probabilities; conditionals applied separately. All attempts counted before acceptance.',line=logprob.lineno)),
        release_execution_commit_verified=False,retention_threshold_configuration_verified=False,injection_spin_reference_recovered=False,
        dcc_probe=dict(http_status=200,scientific_document_available=False,message='Requested version of document does not exist.'),
        limits='Source-level availability and implementation verification strengthen the conditional model, but do not identify the configuration or runtime that produced the fixed LVK cumulative file or the old O3 campaign.')
    (BASE/'generation-source-audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
