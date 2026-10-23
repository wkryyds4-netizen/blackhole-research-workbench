"""Render the checked Markdown paper with embedded Chinese text and equations.

Requires reportlab, matplotlib and Pillow. Supply a licensed CJK TrueType font
with --font/--bold-font on other machines; no font files are redistributed.
"""
from pathlib import Path
import argparse
import hashlib
import html
import json
import re

import matplotlib
matplotlib.use('Agg')
from matplotlib.mathtext import math_to_image
from matplotlib.font_manager import FontProperties
from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, KeepTogether

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


SUPER = dict(zip('⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁿ', '0123456789+-=()n'))
SUB = dict(zip('₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ᵢⱼₖ', '0123456789+-=()ijk'))


def inline(text):
    # CJK fonts can silently omit mathematical Unicode. Preserve its meaning
    # with explicit baseline shifts using ordinary glyphs, including negatives.
    text = re.sub('[\u2010-\u2015\u2212]', '-', text)
    text = html.escape(text)
    text = re.sub(r'\[([^\]]+)\]\(([^)]+)\)', r'<link href="\2" color="#14567a">\1</link>', text)
    text = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', text)
    for mapping, tag in ((SUPER, 'super'), (SUB, 'sub')):
        text = re.sub('[' + ''.join(mapping) + ']+',
                      lambda m: '<' + tag + '>' + ''.join(mapping[c] for c in m[0]) + '</' + tag + '>', text)
    return text.replace('`', '')


def verify_font_coverage(source):
    prose = '\n'.join(line for line in source.splitlines()
                      if not line.startswith(('$$', '![')))
    visible = html.unescape(re.sub('<[^>]+>', '', inline(prose)))
    for name in ('PaperCJK', 'PaperBold'):
        glyphs = pdfmetrics.getFont(name).face.charToGlyph
        missing = sorted({c for c in visible if not c.isspace() and ord(c) not in glyphs})
        if missing:
            raise ValueError(f'{name} lacks visible glyphs: {missing!r}; select another font.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--font', type=Path, default=Path('C:/Windows/Fonts/simsun.ttc'))
    parser.add_argument('--bold-font', type=Path, default=Path('C:/Windows/Fonts/simhei.ttf'))
    parser.add_argument('--output', type=Path, default=ROOT/'output/pdf/blackhole-population-robustness-v0.2.pdf')
    args = parser.parse_args()
    if not args.font.is_file():
        parser.error('Supply --font pointing to a licensed CJK font with TrueType outlines.')
    pdfmetrics.registerFont(TTFont('PaperCJK', str(args.font), subfontIndex=0))
    bold = args.bold_font if args.bold_font.is_file() else args.font
    pdfmetrics.registerFont(TTFont('PaperBold', str(bold), subfontIndex=0))
    pdfmetrics.registerFontFamily('PaperCJK', normal='PaperCJK', bold='PaperBold', italic='PaperCJK', boldItalic='PaperBold')
    base = dict(fontName='PaperCJK', fontSize=10.2, leading=17, wordWrap='CJK', textColor=colors.HexColor('#17232c'),
                spaceAfter=7, allowWidows=0, allowOrphans=0)
    styles = {'body': ParagraphStyle('Body', **base),
              'small': ParagraphStyle('Small', **{**base, 'fontSize': 8.4, 'leading': 12.4, 'spaceAfter': 5}),
              'caption': ParagraphStyle('Caption', **{**base, 'fontSize': 8.8, 'leading': 13.3, 'textColor': colors.HexColor('#394e5a')}),
              'title': ParagraphStyle('Title', **{**base, 'fontName': 'PaperBold', 'fontSize': 19, 'leading': 27,
                                                 'alignment': TA_CENTER, 'spaceAfter': 12, 'keepWithNext': True}),
              'h2': ParagraphStyle('H2', **{**base, 'fontName': 'PaperBold', 'fontSize': 13, 'leading': 20,
                                           'spaceBefore': 12, 'spaceAfter': 7, 'keepWithNext': True}),
              'h3': ParagraphStyle('H3', **{**base, 'fontName': 'PaperBold', 'fontSize': 11, 'leading': 18,
                                           'spaceBefore': 8, 'spaceAfter': 5, 'keepWithNext': True})}
    src = HERE/'manuscript.zh.md'
    text = src.read_text(encoding='utf-8')
    verify_font_coverage(text)
    lines = text.splitlines()
    tmp = HERE/'tmp'; tmp.mkdir(exist_ok=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    width = A4[0]-108
    story = []
    i = 0; equation_count = 0; figure_count = 0; table_count = 0
    while i < len(lines):
        line = lines[i].strip()
        i += 1
        if not line:
            continue
        if line.startswith('$$'):
            equation_count += 1
            formula = line.strip('$').strip().replace(r'\mathcal L', r'\mathcal{L}')
            path = tmp/f'equation-{equation_count:02d}.png'
            math_to_image('$'+formula+'$', str(path), prop=FontProperties(size=13), dpi=230, color='#17232c')
            with PILImage.open(path) as im:
                w, h = im.size
            rendered_width = min(width, w*72/230)
            img = Image(str(path), width=rendered_width, height=h*rendered_width/w)
            img.hAlign = 'CENTER'
            story += [Spacer(1, 4), img, Spacer(1, 11)]
        elif line.startswith('!['):
            match = re.match(r'!\[.*?\]\((.*?)\)', line)
            path = HERE/match[1]
            with PILImage.open(path) as im:
                w, h = im.size
            rendered_width = min(width, w*72/240)
            img = Image(str(path), width=rendered_width, height=h*rendered_width/w)
            img.hAlign = 'CENTER'
            figure_count += 1
            while i < len(lines) and not lines[i].strip(): i += 1
            caption = []
            if i < len(lines) and re.match(r'图\d', lines[i]):
                caption = [Spacer(1, 5), Paragraph(inline(lines[i]), styles['caption'])]
                i += 1
            story.append(KeepTogether([Spacer(1, 6), img]+caption+[Spacer(1, 7)]))
        elif line.startswith('|'):
            rows = [line]
            while i < len(lines) and lines[i].strip().startswith('|'):
                rows.append(lines[i].strip()); i += 1
            parsed = [[v.strip() for v in row.strip('|').split('|')] for row in rows if not re.match(r'^\|[\s:|\-]+\|$', row)]
            n = len(parsed[0]); table_count += 1
            if n == 2: widths = [width*.4, width*.6]
            elif n == 3: widths = [width*.3, width*.56, width*.14]
            elif n == 4: widths = [width*.13, width*.22, width*.22, width*.43]
            else: widths = [width*.27]+[width*.73/(n-1)]*(n-1)
            cells = [[Paragraph(inline(v), styles['small']) for v in row] for row in parsed]
            table = Table(cells, colWidths=widths, repeatRows=1, hAlign='LEFT')
            table.setStyle(TableStyle([('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#e8eff3')),
                                      ('LINEBELOW', (0, 0), (-1, 0), .6, colors.HexColor('#567080')),
                                      ('LINEBELOW', (0, -1), (-1, -1), .4, colors.HexColor('#8ca0ac')),
                                      ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f7f9fa')]),
                                      ('VALIGN', (0, 0), (-1, -1), 'TOP'), ('LEFTPADDING', (0, 0), (-1, -1), 6),
                                      ('RIGHTPADDING', (0, 0), (-1, -1), 6), ('TOPPADDING', (0, 0), (-1, -1), 5),
                                      ('BOTTOMPADDING', (0, 0), (-1, -1), 4)]))
            story += [table, Spacer(1, 10)]
        elif line.startswith('# '):
            story.append(Paragraph(inline(line[2:]), styles['title']))
        elif line.startswith('## '):
            story.append(Paragraph(inline(line[3:]), styles['h2']))
        elif line.startswith('### '):
            story.append(Paragraph(inline(line[4:]), styles['h3']))
        else:
            paragraph = [line]
            while i < len(lines) and lines[i].strip() and not lines[i].startswith(('#', '|', '![', '$$')):
                paragraph.append(lines[i].strip()); i += 1
            value = ' '.join(paragraph)
            style = styles['small'] if re.match(r'^\[\d+\]', value) or len(value.strip('`')) == 64 else styles['body']
            story.append(Paragraph(inline(value), style))
    def page(canvas, doc):
        canvas.saveState()
        canvas.setFont('Helvetica', 8)
        canvas.setFillColor(colors.HexColor('#657a86'))
        canvas.drawString(54, A4[1]-31, 'BLACKHOLE RESEARCH WORKBENCH  |  RESEARCH DRAFT v0.2')
        canvas.setStrokeColor(colors.HexColor('#d8e1e6'))
        canvas.line(54, A4[1]-39, A4[0]-54, A4[1]-39)
        canvas.drawString(54, 30, '2026-10-24  |  Conditional population diagnostics')
        canvas.drawRightString(A4[0]-54, 30, str(doc.page))
        canvas.restoreState()
    doc = SimpleDocTemplate(str(args.output), pagesize=A4, leftMargin=54, rightMargin=54, topMargin=54, bottomMargin=49,
                            title='公开双黑洞群体推断的条件稳健性', author='Blackhole Research Workbench project research draft',
                            subject='Conditional population robustness, prior support, injection provenance and paired Monte Carlo diagnostics')
    doc.build(story, onFirstPage=page, onLaterPages=page)
    metadata = dict(manuscript_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),
                    pdf_sha256=hashlib.sha256(args.output.read_bytes()).hexdigest(),
                    equations=equation_count, figures=figure_count, tables=table_count,
                    inline_font_coverage='passed for regular and bold fonts after semantic sub/superscript conversion',
                    output=str(args.output.relative_to(ROOT)).replace('\\', '/'),
                    note='CJK fonts embedded as subsets from the build host; no font file redistributed. Visual verification recorded separately.')
    (HERE/'build-manifest.json').write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding='utf-8')
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
