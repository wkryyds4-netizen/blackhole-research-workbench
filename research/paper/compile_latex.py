"""Compile both reviewed LaTeX projects and record the exported PDF hashes."""
from pathlib import Path
import hashlib
import json
import re
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    source = json.loads((HERE/'latex-source-manifest.json').read_text())
    editions = []
    for item in source['editions']:
        lang = item['language']
        project = HERE/f'latex-{lang}'
        assert sha(project/'main.tex') == item['main_tex_sha256'], 'Regenerate source manifest first'
        engine = shutil.which(item['engine'])
        if not engine:
            raise RuntimeError(f"Install a TeX distribution and put {item['engine']} on PATH")
        version = subprocess.check_output([engine, '--version'], text=True).splitlines()[0]
        for _ in range(2):
            result = subprocess.run([engine, '-interaction=nonstopmode', '-halt-on-error',
                                     '-file-line-error', 'main.tex'], cwd=project,
                                    capture_output=True, text=True, errors='replace')
            if result.returncode:
                raise RuntimeError(f'{lang}: compilation failed; inspect {project / "main.log"}')
        log = (project/'main.log').read_text(errors='replace')
        forbidden = ('undefined references', 'undefined citations', 'Missing character:',
                     'Overfull', 'Rerun to get', 'Label(s) may have changed')
        assert not any(token in log for token in forbidden), f'{lang}: unresolved layout/reference diagnostics'
        pages = int(re.search(r'Output written on main.pdf \((\d+) pages', log)[1])
        output = ROOT/'output/pdf'/f'blackhole-population-robustness-latex-{lang}-v0.2.pdf'
        output.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(project/'main.pdf', output)
        editions.append({**item, 'engine_version': version, 'pages': pages,
                         'output': output.relative_to(ROOT).as_posix(), 'pdf_sha256': sha(output),
                         'log_sha256': sha(project/'main.log'),
                         'figures_sha256': {p.name: sha(p) for p in sorted((project/'figures').glob('*.pdf'))},
                         'passes': 2, 'references_resolved': True,
                         'missing_glyph_reports': 0, 'overfull_boxes': 0})
    report = {'date': '2026-10-24', 'status': 'compiled', 'editions': editions,
              'note': 'Actual TeX engine outputs. Byte identity across rebuilds is not promised. Reinspect rebuilt PDFs.'}
    (HERE/'latex-build-manifest.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
