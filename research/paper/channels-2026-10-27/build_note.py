"""Build self-contained bilingual TeX projects and actual rendered PDFs."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys

BASE=Path(__file__).resolve().parent;ROOT=BASE.parents[2]
CHANNEL=ROOT/'research/population-channels'

def sha(p):
    with p.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--render',action='store_true');args=parser.parse_args()
    evidence=json.loads((BASE/'evidence-manifest.json').read_text())
    assert sha(CHANNEL/'branching-synthesis.json')==evidence['synthesis_sha256']
    for name,digest in evidence['source_sha256'].items():assert sha(CHANNEL/name)==digest
    for name,digest in evidence['generated_sha256'].items():assert sha(BASE/name)==digest
    editions=[]
    for lang,engine_name in [('en','pdflatex'),('zh','xelatex')]:
        project=BASE/f'latex-{lang}';assert (project/'main.tex').exists()
        for name in ['preamble.tex','references.tex','facts.tex']:shutil.copyfile(BASE/name,project/name)
        shutil.copyfile(BASE/f'table-main-{lang}.tex',project/'table-main.tex')
        (project/'figures').mkdir(exist_ok=True)
        for p in (BASE/'figures').glob('*.pdf'):shutil.copyfile(p,project/'figures'/p.name)
        engine=shutil.which(engine_name);assert engine,engine_name
        for attempt in range(3):
            result=subprocess.run([engine,'-interaction=nonstopmode','-halt-on-error','-file-line-error','main.tex'],cwd=project,capture_output=True,text=True,errors='replace')
            (project/'compile-console.txt').write_text(result.stdout+'\n'+result.stderr,encoding='utf-8')
            if result.returncode:break
            log=(project/'main.log').read_text(errors='replace')
            if attempt>=1 and not any(t in log for t in ['Rerun to get','Label(s) may have changed','undefined references','undefined citations']):break
        log=(project/'main.log').read_text(errors='replace') if (project/'main.log').exists() else result.stdout
        errors=[x for x in log.splitlines() if any(t in x for t in ['Overfull','Missing character:','undefined references','undefined citations','Rerun to get','Label(s) may have changed'])]
        entry=dict(language=lang,engine=engine_name,exit_code=result.returncode,passes=attempt+1,clean_diagnostics=result.returncode==0 and not errors,diagnostics=errors,
                   source_sha256={p.name:sha(p) for p in project.glob('*.tex')},figure_sha256={p.name:sha(p) for p in (project/'figures').glob('*.pdf')})
        if result.returncode==0:
            output=ROOT/'output/pdf'/f'blackhole-channel-completion-latex-{lang}-2026-10-27.pdf'
            output.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(project/'main.pdf',output)
            info=subprocess.check_output(['pdfinfo',str(output)],text=True,errors='replace')
            entry.update(output=output.relative_to(ROOT).as_posix(),pdf_sha256=sha(output),pages=int(re.search(r'Pages:\s+(\d+)',info)[1]))
            if args.render:
                folder=ROOT/'tmp/pdfs/channel-note'/lang;folder.mkdir(parents=True,exist_ok=True)
                subprocess.run(['pdftoppm','-r','100','-png',str(output),str(folder/'page')],check=True,capture_output=True)
                entry['render_prefix']=folder.relative_to(ROOT).as_posix()+'/page'
        else:print(result.stdout[-4000:],flush=True)
        editions.append(entry)
    report=dict(schema='blackhole-bilingual-channel-note-build/1',status='compiled' if all(e['clean_diagnostics'] for e in editions) else 'needs_repair',editions=editions,evidence_sha256=sha(BASE/'evidence-manifest.json'),scope='Actual installed TeX outputs; this record does not assert visual or scientific review.')
    (BASE/'build-manifest.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({**report,'editions':[{k:v for k,v in e.items() if k not in ['source_sha256','figure_sha256']} for e in editions]},indent=2),flush=True)
    if report['status']!='compiled':sys.exit(1)


if __name__=='__main__':main()
