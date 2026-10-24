"""Check bilingual equations, table numbers, figures, references and hashes.

These deterministic checks complement human/agent translation review; they do
not establish semantic equivalence of prose or validate the underlying science.
Only the Python standard library is required.
"""
from pathlib import Path
import hashlib
import json
import re
import unicodedata

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def tables(source):
    groups = []
    current = []
    for line in source.splitlines() + ['']:
        if line.startswith('|'):
            if not re.fullmatch(r'[| :\-]+', line):
                current.append([cell.strip() for cell in line.strip('|').split('|')])
        elif current:
            groups.append(current)
            current = []
    return groups


def numbers(cell):
    cell = unicodedata.normalize('NFKC', cell).replace('−', '-')
    return re.findall(r'-?\d+(?:\.\d+)?', cell)


def main():
    sources = {lang: (HERE/f'manuscript.{lang}.md').read_text(encoding='utf-8') for lang in ('zh', 'en')}
    zh, en = sources['zh'], sources['en']
    equations = lambda source: [line for line in source.splitlines() if line.startswith('$$')]
    assert len(equations(zh)) == 14 and equations(zh) == equations(en), 'Equation mismatch'
    sections = lambda source: re.findall(r'^#{2,3} (\d+(?:\.\d+)?) ', source, re.M)
    assert sections(zh) == sections(en), 'Numbered section mismatch'
    a, b = tables(zh), tables(en)
    assert len(a) == len(b) == 6, 'Table count mismatch'
    compared_cells = 0
    for index, (left, right) in enumerate(zip(a, b), 1):
        assert len(left) == len(right), f'Table {index}: row mismatch'
        for row, (lrow, rrow) in enumerate(zip(left, right)):
            assert len(lrow) == len(rrow), f'Table {index}: column mismatch'
            if row == 0:
                continue
            for col, (lcell, rcell) in enumerate(zip(lrow, rrow), 1):
                assert numbers(lcell) == numbers(rcell), f'Table {index}, row {row}, column {col}: numeric mismatch'
                compared_cells += 1
    figures = lambda source: re.findall(r'^!\[.*?\]\((.*?)\)', source, re.M)
    assert figures(zh) == figures(en) and len(figures(en)) == 3, 'Figure mismatch'
    links = lambda source: re.findall(r'(?<!!)\[[^\]]+\]\((https?://[^)]+)\)', source)
    assert links(zh) == links(en), 'Reference URL mismatch'
    refs = lambda source: re.findall(r'^\[(\d+)\]', source, re.M)
    assert refs(zh) == refs(en) == list(map(str, range(1, 10))), 'Reference numbering mismatch'
    for lang, name in [('zh', 'build-manifest.json'), ('en', 'build-manifest.en.json')]:
        manifest = json.loads((HERE/name).read_text(encoding='utf-8'))
        assert hashlib.sha256((HERE/f'manuscript.{lang}.md').read_bytes()).hexdigest() == manifest['manuscript_sha256']
        assert hashlib.sha256((ROOT/manifest['output']).read_bytes()).hexdigest() == manifest['pdf_sha256']
    report = {'status': 'pass', 'date': '2026-10-24',
              'manuscript_sha256': {lang: hashlib.sha256((HERE/f'manuscript.{lang}.md').read_bytes()).hexdigest() for lang in sources},
              'identical_equations': 14, 'matched_tables': 6, 'compared_body_cells': compared_cells,
              'identical_figure_paths': 3, 'matched_numbered_references': 9,
              'numbered_sections_match': True, 'reference_urls_match': True,
              'pdf_and_source_hashes_match': True,
              'scope': 'Structural/numeric parity, not prose semantic validation. See translation-review.md for internal translation review.'}
    (HERE/'translation-validation.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
