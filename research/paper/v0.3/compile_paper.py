"""Compile the new bilingual revision with installed TeX engines.

No old manuscript, baseline result or v0.2 PDF is overwritten. Rendering is
for visual inspection; a successful compile does not assert scientific review.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
TMP=ROOT/'research/paper/tmp/v0.3'


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--render',action='store_true');args=parser.parse_args()
    evidence=json.loads((BASE/'evidence-manifest.json').read_text())
    for name,digest in evidence['source_sha256'].items():assert sha(ROOT/name)==digest,name
    for name,digest in evidence['generated_sha256'].items():assert sha(BASE/name)==digest,name
    editions=[]
    for lang,engine_name in [('en','pdflatex'),('zh','xelatex')]:
        project=BASE/f'latex-{lang}';engine=shutil.which(engine_name);assert engine,engine_name
        version=subprocess.check_output([engine,'--version'],text=True,errors='replace').splitlines()[0]
        for run in range(3):
            cp=subprocess.run([engine,'-interaction=nonstopmode','-halt-on-error','-file-line-error','main.tex'],
                cwd=project,capture_output=True,text=True,errors='replace')
            (project/'compile-console.txt').write_text(cp.stdout+'\n'+cp.stderr,encoding='utf-8')
            if cp.returncode:break
            log=(project/'main.log').read_text(errors='replace')
            if run>=1 and not any(x in log for x in ['Rerun to get','Label(s) may have changed','undefined references','undefined citations']):break
        log=(project/'main.log').read_text(errors='replace') if (project/'main.log').exists() else cp.stdout
        problems=[line for line in log.splitlines() if any(x in line for x in ['Overfull','Missing character:','undefined references','undefined citations','Rerun to get','Label(s) may have changed'])]
        clean=cp.returncode==0 and not problems
        match=re.search(r'Output written on main.pdf \((\d+) pages',log)
        entry=dict(language=lang,engine=engine_name,engine_version=version,passes=run+1,exit_code=cp.returncode,
            pages=int(match[1]) if match else None,clean_diagnostics=clean,diagnostics=problems,
            source_sha256={p.relative_to(BASE).as_posix():sha(p) for p in sorted(project.glob('*.tex'))},
            figures_sha256={p.name:sha(p) for p in sorted((project/'figures').glob('*.pdf'))},
            log_sha256=sha(project/'main.log') if (project/'main.log').exists() else None)
        if clean:
            output=ROOT/'output/pdf'/f'blackhole-population-robustness-latex-{lang}-v0.3.pdf'
            output.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(project/'main.pdf',output)
            entry.update(output=output.relative_to(ROOT).as_posix(),pdf_sha256=sha(output))
        if args.render and (project/'main.pdf').exists() and cp.returncode==0:
            renderer=shutil.which('pdftoppm');assert renderer,'pdftoppm'
            folder=TMP/lang;folder.mkdir(parents=True,exist_ok=True)
            subprocess.run([renderer,'-r','105','-png',str(project/'main.pdf'),str(folder/'page')],check=True,capture_output=True)
            entry['render_prefix']=folder.relative_to(ROOT).as_posix()+'/page'
        if cp.returncode:print(cp.stdout[-5000:])
        editions.append(entry)
    result=dict(schema='blackhole-bilingual-v03-build/1',date='2026-10-25',
        status='compiled' if all(e['clean_diagnostics'] for e in editions) else 'needs_repair',
        evidence_manifest_sha256=sha(BASE/'evidence-manifest.json'),editions=editions,
        note='Actual installed TeX engine outputs. The desktop editor compiled English but its Chinese fallback lacked a font; direct XeLaTeX avoids that fallback. No visual/semantic review is asserted by this build record.')
    (BASE/'build-manifest.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({**result,'editions':[{k:v for k,v in e.items() if k not in ['source_sha256','figures_sha256']} for e in editions]},indent=2))
    if result['status']!='compiled':sys.exit(1)


if __name__=='__main__':main()
