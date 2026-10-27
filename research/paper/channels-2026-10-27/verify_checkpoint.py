"""Check saved evidence, TeX/PDF and workbench hashes without running research."""
from pathlib import Path
import hashlib
import json
import re

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
CHANNEL = ROOT / 'research/population-channels'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def main():
    evidence = read(BASE / 'evidence-manifest.json')
    assert sha(CHANNEL / 'branching-synthesis.json') == evidence['synthesis_sha256']
    for name, digest in evidence['source_sha256'].items():
        assert sha(CHANNEL / name) == digest, name
    for name, digest in evidence['generated_sha256'].items():
        assert sha(BASE / name) == digest, name
    build = read(BASE / 'build-manifest.json')
    assert build['status'] == 'compiled'
    assert sha(BASE / 'evidence-manifest.json') == build['evidence_sha256']
    structures = []
    for edition in build['editions']:
        project = BASE / ('latex-' + edition['language'])
        assert edition['clean_diagnostics'] and edition['pages'] == 7
        for name, digest in edition['source_sha256'].items():
            assert sha(project / name) == digest, str(project / name)
        for name, digest in edition['figure_sha256'].items():
            assert sha(project / 'figures' / name) == digest, name
            assert sha(BASE / 'figures' / name) == digest, name
        output = ROOT / edition['output']
        assert sha(output) == edition['pdf_sha256'], str(output)
        assert sha(ROOT / 'simulator/public/research-docs' / output.name) == edition['pdf_sha256']
        for name in ['facts.tex', 'preamble.tex', 'references.tex']:
            assert sha(project / name) == sha(BASE / name), name
        assert sha(project / 'table-main.tex') == sha(BASE / ('table-main-' + edition['language'] + '.tex'))
        text = (project / 'manuscript.tex').read_text(encoding='utf-8')
        labels = re.findall(r'\\label\{([^}]+)\}', text)
        citations = sorted(set(re.findall(r'\\cite\w*\{([^}]+)\}', text)))
        structures.append((labels, citations))
    assert structures[0] == structures[1], 'Bilingual labels/citations differ'
    export = read(CHANNEL / 'workbench-export.json')
    assert sha(ROOT / 'simulator/public/datasets/channel-study.json') == export['sha256']
    assert sha(CHANNEL / 'export_workbench.py') == export['exporter_sha256']
    for source in export['source_reports']:
        assert sha(ROOT / source['repo_path']) == source['sha256'], source['repo_path']
    for paper in export['papers']:
        assert sha(ROOT / 'simulator/public' / paper['url']) == paper['sha256']
    print('PASS: saved evidence, bilingual TeX/PDF, shared figures and workbench hashes; no research rerun.')


if __name__ == '__main__':
    main()
