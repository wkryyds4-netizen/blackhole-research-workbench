"""Inspect PDF text/geometry and prepare per-page contact sheets.

Run with the bundled PDF Python environment. This creates inspection evidence,
not an automatic claim of visual or scientific correctness.
"""
from pathlib import Path
import hashlib
import json
import re
from collections import Counter
from pypdf import PdfReader
import pdfplumber
from PIL import Image, ImageOps, ImageDraw

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
TMP=ROOT/'research/paper/tmp/v0.3'


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    manifest=json.loads((BASE/'build-manifest.json').read_text())
    evidence=json.loads((BASE/'evidence-manifest.json').read_text())
    assert manifest['status']=='compiled' and manifest['evidence_manifest_sha256']==sha(BASE/'evidence-manifest.json')
    sources=[];editions=[]
    for entry in manifest['editions']:
        lang=entry['language'];project=BASE/f'latex-{lang}';path=ROOT/entry['output']
        assert sha(path)==entry['pdf_sha256']
        for name,digest in entry['source_sha256'].items():assert sha(BASE/name)==digest,name
        for name,digest in entry['figures_sha256'].items():assert sha(project/'figures'/name)==digest==evidence['figure_sha256'][name]
        body=(project/'manuscript.tex').read_text(encoding='utf-8')
        refs=(project/'bibliography.tex').read_text(encoding='utf-8')
        labels=re.findall(r'\\label\{([^}]+)\}',body);assert len(labels)==len(set(labels))
        assert set(re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',body))<=set(labels)
        bibliography=re.findall(r'\\bibitem\{([^}]+)\}',refs);assert len(bibliography)==len(set(bibliography))==15
        citations=[ref for group in re.findall(r'\\cite\{([^}]+)\}',body) for ref in group.split(',')]
        assert set(citations)<=set(bibliography)
        fact_uses=Counter(re.findall(r'\\(Data[A-Za-z]+)',body));assert set(fact_uses)<=set(evidence['facts'])
        sources.append(dict(fact_uses=dict(fact_uses),labels=labels,citations=citations,
            equations=body.count(r'\begin{equation}'),figures=body.count(r'\begin{figure}'),tables=body.count(r'\begin{table}')))
        reader=PdfReader(path);assert len(reader.pages)==entry['pages']
        assert 'Chengheng' in (reader.metadata.author or '')
        pages=[]
        with pdfplumber.open(path) as pdf:
            for n,page in enumerate(pdf.pages,1):
                text=page.extract_text() or '';assert len(text)>70,(lang,n,'Unexpected sparse page')
                assert '\ufffd' not in text,(lang,n,'Replacement glyph')
                outside=[dict(text=w['text'],x0=w['x0'],x1=w['x1'],top=w['top'],bottom=w['bottom']) for w in page.extract_words()
                    if w['x0']<25 or w['x1']>page.width-25 or w['top']<20 or w['bottom']>page.height-20]
                assert not outside,(lang,n,outside)
                render=TMP/lang/f'page-{n:02d}.png';assert render.exists()
                pages.append(dict(page=n,characters=len(text),width=page.width,height=page.height,
                    text_inside_page_safety_margins=True,render_sha256=sha(render)))
        # Four-page sheets are an overview only; inspect relevant pages at
        # original rendered resolution as well, especially formulas/tables.
        sheets=[]
        for start in range(1,len(pages)+1,4):
            sheet=Image.new('RGB',(1000,1460),'#dedede');draw=ImageDraw.Draw(sheet)
            selected=[]
            for offset,n in enumerate(range(start,min(start+4,len(pages)+1))):
                image=Image.open(TMP/lang/f'page-{n:02d}.png').convert('RGB')
                image.thumbnail((480,690));x=10+(offset%2)*500;y=25+(offset//2)*730
                sheet.paste(image,(x,y));draw.text((x,y-17),f'{lang} page {n}',fill='black');selected.append(n)
            target=TMP/lang/f'contact-{start:02d}.png';sheet.save(target)
            sheets.append(dict(pages=selected,path=target.relative_to(ROOT).as_posix(),sha256=sha(target)))
        editions.append(dict(language=lang,pdf_sha256=entry['pdf_sha256'],pages=pages,contact_sheets=sheets))
    assert sources[0]==sources[1],'Bilingual fact usage, equations, figures, tables, labels or citations differ'
    for name in ['g-table.tex','finite-table.tex','support-table.tex','completion-table.tex']:
        numerical=[]
        for lang in ['en','zh']:
            # Compare actual data rows, excluding the translated header.
            content=(BASE/f'latex-{lang}'/name).read_text(encoding='utf-8').split(r'\midrule')[1]
            numerical.append(re.findall(r'-?\d+(?:\.\d+)?',content))
        assert numerical[0]==numerical[1],name
    result=dict(schema='blackhole-bilingual-v03-text-geometry-inspection/1',status='text_and_geometry_checked_visual_review_required',
        build_manifest_sha256=sha(BASE/'build-manifest.json'),source_structure=sources[0],editions=editions,
        scope='Source/translation structure, numeric tables, PDF text, broad page margins and rendered-file inventory. Not a peer review or an automatic visual-quality certificate.')
    (BASE/'inspection.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(dict(status=result['status'],editions=[dict(language=e['language'],pages=len(e['pages'])) for e in editions],source_structure=sources[0]),indent=2))


if __name__=='__main__':main()
