"""Offline release integrity check; does not claim to recompile TeX in CI."""
from pathlib import Path
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest = json.loads((HERE/'latex-build-manifest.json').read_text())
    sources = json.loads((HERE/'latex-source-manifest.json').read_text())
    qa = json.loads((HERE/'latex-visual-qa.json').read_text())
    assert qa['status'] == 'pass'
    assert len(manifest['editions']) == len(sources['editions']) == len(qa['editions']) == 2
    for item, source in zip(manifest['editions'], sources['editions']):
        lang = item['language']
        assert lang == source['language']
        project = HERE/f'latex-{lang}'
        text = (project/'main.tex').read_text(encoding='utf-8')
        assert sha(project/'main.tex') == item['main_tex_sha256'] == source['main_tex_sha256']
        assert sha(HERE/f'manuscript.{lang}.md') == item['manuscript_sha256'] == source['manuscript_sha256']
        assert sha(ROOT/item['output']) == item['pdf_sha256']
        inspected = next(e for e in qa['editions'] if e['language'] == lang)
        assert inspected['pdf_sha256'] == item['pdf_sha256'] and inspected['pages'] == item['pages']
        assert (ROOT/item['output']).read_bytes().startswith(b'%PDF-')
        for env, count in [('equation', 14), ('figure', 3), ('table', 6)]:
            assert text.count('\\begin{'+env+'}') == count
        bib = set(re.findall(r'\\bibitem\{([^}]+)\}', text))
        assert bib == {'ref'+str(i) for i in range(1, 10)}
        for refs in re.findall(r'\\cite\{([^}]+)\}', text):
            assert set(refs.split(',')) <= bib, f'Invalid citation {refs}'
        labels = set(re.findall(r'\\label\{([^}]+)\}', text))
        assert set(re.findall(r'\\(?:eqref|ref)\{([^}]+)\}', text)) <= labels
        assert r'\([0,1]\)' in text, 'The prior interval must remain an interval'
        for name, expected in item['figures_sha256'].items():
            assert sha(project/'figures'/name) == expected == sha(HERE/'figures'/name)
    print('PASS: both LaTeX sources, figures, reference targets and compiled PDF hashes match the release.')


if __name__ == '__main__':
    main()
