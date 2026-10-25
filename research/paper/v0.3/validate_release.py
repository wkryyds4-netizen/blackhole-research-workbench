"""Offline integrity gate for the versioned bilingual manuscript artifacts.

This verifies the evidence chain and bounded review records. It is not a
scientific reanalysis, a fresh PDF compilation, or independent peer review.
"""
from pathlib import Path
import hashlib
import json
import re

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def read(name):return json.loads((BASE/name).read_text(encoding='utf-8'))


def main():
    evidence=read('evidence-manifest.json');build=read('build-manifest.json');inspection=read('inspection.json')
    qa=read('visual-qa.json');translation=read('translation-review.json')
    assert build['status']=='compiled'
    assert evidence['producer_sha256']==sha(BASE/'build_evidence.py')
    assert build['evidence_manifest_sha256']==sha(BASE/'evidence-manifest.json')
    assert inspection['build_manifest_sha256']==qa['build_manifest_sha256']==sha(BASE/'build-manifest.json')
    assert qa['inspection_sha256']==sha(BASE/'inspection.json')
    assert qa['status']=='pass' and translation['status']=='self_review_completed'
    assert qa['independent_peer_review'] is False and translation['independent_peer_review'] is False
    for name,digest in evidence['source_sha256'].items():assert sha(ROOT/name)==digest,name
    for name,digest in evidence['generated_sha256'].items():assert sha(BASE/name)==digest,name
    for name,digest in qa['review_tool_sha256'].items():assert sha(BASE/name)==digest,name
    structures=[]
    for entry in build['editions']:
        lang=entry['language'];project=BASE/f'latex-{lang}';pdf=ROOT/entry['output']
        assert entry['exit_code']==0 and entry['clean_diagnostics'] is True and not entry['diagnostics']
        assert pdf.read_bytes().startswith(b'%PDF-') and sha(pdf)==entry['pdf_sha256']
        for name,digest in entry['source_sha256'].items():assert sha(BASE/name)==digest,name
        for name,digest in entry['figures_sha256'].items():
            assert sha(project/'figures'/name)==sha(BASE/'figures'/name)==digest==evidence['figure_sha256'][name]
        q=next(x for x in qa['editions'] if x['language']==lang)
        inspected=next(x for x in inspection['editions'] if x['language']==lang)
        assert q['pdf_sha256']==inspected['pdf_sha256']==entry['pdf_sha256']
        assert q['pages']==len(inspected['pages'])==entry['pages']==14
        assert q['all_pages_checked_in_overviews'] is True
        assert q['native_resolution_pages']==({'en':[4,8,9,11],'zh':[3,7,9,11]}[lang])
        assert set(q['native_resolution_pages'])<=set(range(1,entry['pages']+1))
        for name,digest in translation['source_sha256'][lang].items():assert sha(BASE/name)==digest,name
        body=(project/'manuscript.tex').read_text(encoding='utf-8')
        structures.append((re.findall(r'\\label\{([^}]+)\}',body),re.findall(r'\\cite\{([^}]+)\}',body),
            body.count(r'\begin{equation}'),body.count(r'\begin{figure}'),body.count(r'\begin{table}')))
    assert structures[0]==structures[1] and structures[0][2:]==(14,5,4)
    assert (BASE/'latex-en/facts.tex').read_bytes()==(BASE/'latex-zh/facts.tex').read_bytes()
    assert (BASE/'latex-en/bibliography.tex').read_bytes()==(BASE/'latex-zh/bibliography.tex').read_bytes()
    paired=json.loads((ROOT/'research/population-invariant/paired-pe-validation.json').read_text())
    assert paired['original_step_derivative_gate_passed'] is False
    print('PASS: v0.3 bilingual source, validated-result links, figures, PDFs and bounded self-review records match. No publication or physical-robustness certification is implied.')


if __name__=='__main__':main()
