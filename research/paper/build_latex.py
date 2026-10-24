"""Generate complete, editable LaTeX projects from the reviewed bilingual text.

The manuscripts remain the content source; equations use native LaTeX, figures
use vector PDFs, and references are embedded to avoid a missing BibTeX step.
"""
from pathlib import Path
import hashlib
import json
import re
import shutil

HERE = Path(__file__).resolve().parent
SUPER = dict(zip('⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁿ', '0123456789+-=()n'))
SUB = dict(zip('₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ᵢⱼₖ', '0123456789+-=()ijk'))
GREEK = {'α': r'\alpha', 'β': r'\beta', 'θ': r'\theta', 'σ': r'\sigma',
         'λ': r'\lambda', 'μ': r'\mu', 'π': r'\pi', 'ξ': r'\xi',
         'ε': r'\epsilon', 'Φ': r'\Phi'}
SYMBOL = {'≤': r'\leq', '≥': r'\geq', '≈': r'\approx', '×': r'\times',
          '−': '-', '°': r'{}^\circ', '∞': r'\infty', '∝': r'\propto'}


def inline(value):
    protected = []
    def hold(tex):
        protected.append(tex)
        return f'ZZLATEX{len(protected)-1}ZZ'
    value = re.sub(r'\[([^\]]+)\]\((https?://[^)]+)\)',
                   lambda m: hold(r'\href{' + m[2].replace('%', r'\%').replace('#', r'\#') + '}{' + inline(m[1]) + '}'), value)
    value = re.sub(r'`([^`]+)`', lambda m: hold(r'\path{' + m[1] + '}'), value)
    value = re.sub(r'\*\*(.*?)\*\*', lambda m: hold(r'\textbf{' + inline(m[1]) + '}'), value)
    value = value.replace('[0,1]', hold(r'\([0,1]\)'))
    value = re.sub(r'\[(\d+(?:,\d+)*)\]', lambda m: hold(r'\cite{' + ','.join('ref'+v for v in m[1].split(',')) + '}'), value)
    value = re.sub(r'Equations \((\d+)\) and \((\d+)\)', lambda m: hold(r'Equations~\eqref{eq:'+m[1]+r'} and~\eqref{eq:'+m[2]+'}'), value)
    value = re.sub(r'(Equation |式)\((\d+)\)', lambda m: hold(m[1].rstrip()+r'~\eqref{eq:'+m[2]+'}'), value)
    value = re.sub(r'Section (\d+(?:\.\d+)?)', lambda m: hold(r'Section~\ref{sec:'+m[1]+'}'), value)
    value = re.sub(r'第(\d+(?:\.\d+)?)节', lambda m: hold(r'第~\ref{sec:'+m[1]+r'}~节'), value)
    value = re.sub(r'Table (\d+)', lambda m: hold(r'Table~\ref{tab:'+m[1]+'}'), value)
    value = re.sub(r'表(\d+)', lambda m: hold(r'表~\ref{tab:'+m[1]+'}'), value)
    value = re.sub(r'([AB])σ', lambda m: hold(r'\('+m[1]+r'_\sigma\)'), value)
    base = 'A-Za-z' + ''.join(GREEK)
    value = re.sub('(['+base+'])(['+''.join(SUB)+']+)',
        lambda m: hold(r'\('+GREEK.get(m[1],m[1])+'_{'+''.join(SUB[c] for c in m[2])+r'}\)'), value)
    value = re.sub('['+''.join(SUPER)+']+',
        lambda m: hold(r'\textsuperscript{'+''.join(SUPER[c] for c in m[0])+'}'), value)
    for char, macro in {**GREEK, **SYMBOL}.items():
        value = value.replace(char, hold(r'\('+macro+r'\)'))
    value = value.replace('·', hold(r'\textperiodcentered{}'))
    replacements = {'\\': r'\textbackslash{}', '&': r'\&', '%': r'\%', '$': r'\$',
                    '#': r'\#', '_': r'\_', '{': r'\{', '}': r'\}', '~': r'\textasciitilde{}', '^': r'\textasciicircum{}'}
    value = ''.join(replacements.get(c,c) for c in value)
    # Long repository paths outside Markdown code spans also need break points.
    value = re.sub(r'(research/[A-Za-z0-9/_.\\-]+\.py)', lambda m: hold(r'\path{'+m[0].replace(r'\_', '_')+'}'), value)
    for index in range(len(protected)-1,-1,-1):
        value = value.replace(f'ZZLATEX{index}ZZ', protected[index])
    return value


CAPTIONS = {
 'en': ['Fixed settings defining the 81 conditional analyses.',
        'Run membership and recovered injection mixture weights.',
        'Paired differences of posterior medians and empirical numerical perturbations (percentage points).',
        'Perturbation standard deviations of individual model medians (percentage points).',
        'Independent numerical validation. Errors are in proportion units unless specified otherwise.',
        'Reproduction entry points for the scientific calculations.'],
 'zh': ['构成81个条件分析的固定设置。', '注入运行归属、计数与恢复的混合权重。',
        '成对后验中位数差及经验数值扰动，单位为百分点。', '各模型后验中位数的扰动标准差，单位为百分点。',
        '独立数值验证。未特别注明的误差以比例为单位。', '科学计算复现入口。']}


def generate(lang):
    source = HERE/f'manuscript.{lang}.md'
    lines = source.read_text(encoding='utf-8').splitlines()
    root = HERE/f'latex-{lang}'
    assert root.is_dir(), 'Create the dedicated LaTeX project before generating sources'
    images = root/'figures'; images.mkdir(exist_ok=True)
    chinese = lang == 'zh'
    if chinese:
        preamble = r'''% !TeX program = xelatex
\documentclass[UTF8,a4paper,11pt,fontset=fandol]{ctexart}
'''
    else:
        preamble = r'''% !TeX program = pdflatex
\documentclass[a4paper,11pt]{article}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern}
'''
    preamble += r'''\usepackage[margin=25mm]{geometry}
\usepackage{amsmath,amssymb,bm}
\usepackage{graphicx,booktabs,tabularx,array}
\usepackage[numbers,sort&compress]{natbib}
\usepackage{microtype}
\usepackage{placeins}
\usepackage[colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue]{hyperref}
\usepackage{url}
\urlstyle{tt}
\setlength{\emergencystretch}{4em}
\setlength{\parskip}{0.25em}
\setlength{\parindent}{1em}
\renewcommand{\arraystretch}{1.18}
\newcolumntype{Y}{>{\raggedright\arraybackslash}X}
\graphicspath{{figures/}}
\hypersetup{pdfauthor={Blackhole Research Workbench project research draft}}
'''
    title = lines[0].removeprefix('# ')
    subtitle = next(l.removeprefix('## ') for l in lines if l.startswith('## '))
    parts = [preamble, r'\title{'+inline(title)+r'\\[0.4em]\large '+inline(subtitle)+'}',
             r'\author{Blackhole Research Workbench\\\small '+('项目研究稿' if chinese else 'Project research draft')+'}',
             r'\date{'+('2026年10月24日' if chinese else '24 October 2026')+r'\\\small v0.2 --- '+('未同行评审' if chinese else 'Not peer reviewed')+'}',
             r'\begin{document}', r'\maketitle']
    status = next(l for l in lines if l.startswith(('稿件状态：','Manuscript status:')))
    parts.append(r'{\small '+inline(status)+r'\par}')
    abstract_header = '### 摘要' if chinese else '### Abstract'
    ai = lines.index(abstract_header)
    abstract = next(l for l in lines[ai+1:] if l.strip())
    parts += [r'\begin{abstract}', inline(abstract), r'\end{abstract}']
    keywords = next(l for l in lines if l.startswith(('**关键词','**Keywords')))
    parts.append(r'\noindent '+inline(keywords)+r'\par')
    if chinese:
        ai = lines.index('### Abstract')
        abstract = next(l for l in lines[ai+1:] if l.strip())
        parts += [r'\section*{Abstract}', inline(abstract)]
    index = next(i for i,l in enumerate(lines) if re.match(r'## 1 ',l))
    equations = figures = tables = references = 0
    bibliography = False
    while index < len(lines):
        line = lines[index].strip(); index += 1
        if not line: continue
        if line.startswith('$$'):
            equations += 1
            formula = re.sub(r'\s*\\qquad\s*\(\d+\)\s*$', '', line.strip('$').strip())
            if equations == 4:
                left,right = formula.split(r',\quad ')
                formula = r'\begin{aligned}'+left.replace('=', '&=',1)+r',\\'+right.replace('=', '&=',1)+r'\end{aligned}'
            parts += [r'\begin{equation}\label{eq:'+str(equations)+'}',formula,r'\end{equation}']
        elif line.startswith('!['):
            figures += 1
            path = re.search(r'\((.*?)\)',line)[1]
            stem = Path(path).stem
            shutil.copyfile(HERE/'figures'/f'{stem}.pdf',images/f'{stem}.pdf')
            while index < len(lines) and not lines[index].strip(): index += 1
            caption = re.sub(r'^(?:图\d+：|Figure \d+\.\s*)','',lines[index]); index += 1
            width = '.96' if figures != 3 else '.86'
            parts += [r'\begin{figure}[tbp]',r'\centering',r'\includegraphics[width='+width+r'\linewidth]{'+stem+'.pdf}',
                      r'\caption{'+inline(caption)+r'}\label{fig:'+str(figures)+'}',r'\end{figure}']
        elif line.startswith('|'):
            tables += 1; rows = [line]
            while index < len(lines) and lines[index].startswith('|'):
                rows.append(lines[index]); index += 1
            cells = [[v.strip() for v in row.strip('|').split('|')] for row in rows if not re.fullmatch(r'[| :\-]+',row)]
            n = len(cells[0]); spec = {2:'YY',3:'Y Y c',4:'l r r Y',5:'Y c c c Y'}[n]
            parts += [r'\begin{table}[tbp]',r'\centering\small',r'\caption{'+CAPTIONS[lang][tables-1]+r'}\label{tab:'+str(tables)+'}',
                      r'\begin{tabularx}{\linewidth}{'+spec+'}',r'\toprule']
            for row,c in enumerate(cells):
                parts.append(' & '.join(inline(v) for v in c)+r' \\')
                if row == 0:parts.append(r'\midrule')
            parts += [r'\bottomrule',r'\end{tabularx}',r'\end{table}']
        elif re.match(r'^## \d+ ',line):
            m = re.match(r'^## (\d+) (.*)',line)
            parts += [r'\FloatBarrier',r'\section{'+inline(m[2])+r'}\label{sec:'+m[1]+'}']
        elif re.match(r'^### \d',line):
            m = re.match(r'^### (\d+\.\d+) (.*)',line)
            parts.append(r'\subsection{'+inline(m[2])+r'}\label{sec:'+m[1]+'}')
        elif line.startswith(('## 附录A','## Appendix A')):
            title = re.sub(r'^## (附录A|Appendix A)\s*','',line)
            parts += [r'\FloatBarrier',r'\appendix',r'\section{'+inline(title)+r'}\label{app:reproduction}']
        elif line in ('## 参考文献','## References'):
            bibliography = True
            parts += [r'\FloatBarrier',r'\begingroup\small',r'\setlength{\bibsep}{2pt}',r'\begin{thebibliography}{99}']
        elif line.startswith('## '):
            parts += [r'\FloatBarrier',r'\section*{'+inline(line[3:])+'}']
        elif bibliography:
            m = re.match(r'^\[(\d+)\] (.*)',line); assert m,line
            references += 1
            parts.append(r'\bibitem{ref'+m[1]+'} '+inline(m[2]))
        else:
            # Keep checksum labels with their values across page boundaries.
            if line.endswith((':', '：')) and 'SHA-256' in line:
                parts.append(inline(line)+r'\par\nopagebreak')
            else:
                parts.append(inline(line)+'\n')
    parts += [r'\end{thebibliography}',r'\endgroup',r'\end{document}']
    assert (equations,figures,tables,references)==(14,3,6,9)
    main = root/'main.tex'; main.write_text('\n'.join(parts)+'\n',encoding='utf-8',newline='\n')
    return {'language':lang,'manuscript_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'main_tex_sha256':hashlib.sha256(main.read_bytes()).hexdigest(),
            'equations':equations,'figures':figures,'tables':tables,'references':references,
            'engine':'xelatex' if chinese else 'pdflatex'}


if __name__ == '__main__':
    report = {'date':'2026-10-24','editions':[generate(lang) for lang in ('zh','en')],
              'format':'Standard article/ctexart preprint layout; not an arXiv endorsement or submission.',
              'guidance':'https://info.arxiv.org/help/submit_tex.html'}
    (HERE/'latex-source-manifest.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(report,indent=2))
