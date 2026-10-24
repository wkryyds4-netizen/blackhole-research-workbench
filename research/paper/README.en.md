# Conditional Robustness of Binary-Black-Hole Population Inference

**Research technical paper v0.2, 24 October 2026. Not peer reviewed.**

[中文说明](README.md)

| Edition | PDF | Editable manuscript |
|---|---|---|
| English | [14-page LaTeX paper](../../output/pdf/blackhole-population-robustness-latex-en-v0.2.pdf) | [Full English text](manuscript.en.md) |
| Chinese | [14-page LaTeX paper](../../output/pdf/blackhole-population-robustness-latex-zh-v0.2.pdf) | [Chinese text](manuscript.zh.md) |

Both editions report the same research: 69 tutorial-derived events, 81 fixed population settings, a shared-support prior control, row-level provenance recovery for 41,972 injections, and 2,048 paired empirical resamples. Each contains 14 displayed equations, six tables, three figures, and nine references. The English edition is a full translation; both LaTeX editions now contain 14 pages.

The baseline spin-amplitude assumption gives a common-prior posterior-median difference of 6.734 percentage points between two tilt-width settings, with a 5th–95th percentile empirical perturbation range of 5.557–7.600 points. The corresponding low-spin stress case gives 1.115 points, with range −0.293–3.601 points. These ranges describe numerical sensitivity, not astrophysical credible intervals. The phenomenological mixture components are not formation-channel labels. Event-level original PE settings and event–injection selection consistency remain unresolved.

## LaTeX projects and compilation

Editable projects: [English pdfLaTeX](latex-en/main.tex) and [Chinese XeLaTeX](latex-zh/main.tex). Each directory is self-contained, with `main.tex` and three vector PDF figures. The preprint uses A4, 11pt type, 25mm margins, native equations, booktabs tables, numbered captions, and linked references. English uses article/Latin Modern; Chinese uses ctexart/Fandol.

The packaging follows the [official arXiv TeX source guidance](https://info.arxiv.org/help/submit_tex.html). This is not a journal-specific template or an arXiv submission. References are embedded in `thebibliography`, so BibTeX is unnecessary. The date is fixed.

Install TeX Live or MiKTeX with ctex/Fandol, natbib, microtype, placeins, and the standard packages; put `xelatex` and `pdflatex` on PATH. From the repository root:

```sh
python -X utf8 research/paper/compile_latex.py
python -X utf8 research/paper/validate_latex.py
```

The script runs each engine twice and exports to `output/pdf/`. After changing the Markdown manuscripts, regenerate sources with `python -X utf8 research/paper/build_latex.py`; after changing figure generation, run `build_figures.py` first. For direct TeX editing, do not rerun the generator because it overwrites `main.tex`; run `pdflatex main.tex` or `xelatex main.tex` twice inside the corresponding project instead.

The [source manifest](latex-source-manifest.json), [build manifest](latex-build-manifest.json), and [visual QA](latex-visual-qa.json) record sources, engines, page counts, and hashes. CI verifies released source/figure/PDF integrity, rather than recompiling TeX. Rebuilds can change PDF bytes and require renewed layout review and QA records.

## Earlier ReportLab exports

The earlier [Chinese 13-page](../../output/pdf/blackhole-population-robustness-v0.2.pdf) and [English 16-page](../../output/pdf/blackhole-population-robustness-en-v0.2.pdf) exports remain available as historical versions. They were not compiled with LaTeX. The following commands reproduce those earlier exports.

### Build earlier exports from frozen results

From the repository root:

```sh
python -m pip install -r research/paper/requirements.txt
python -X utf8 research/paper/build_figures.py
python -X utf8 research/paper/build_pdf.py --language en
python -X utf8 research/paper/build_pdf.py --language zh
python -X utf8 research/paper/validate_bilingual.py
```

Figure and PDF generation use tracked derived results; no raw-data download or inference rerun is required. The observed build environment uses Python 3.14. English uses matplotlib's bundled DejaVu Serif. Chinese defaults to Windows SimSun/SimHei; other hosts must supply legally usable CJK fonts with TrueType outlines:

```sh
python -X utf8 research/paper/build_pdf.py --language zh --font /path/to/cjk-regular.ttf --bold-font /path/to/cjk-bold.ttf
```

The repository does not redistribute font files. PDFs embed font subsets. The builder checks regular and bold glyph coverage, preserves negative signs, and converts Unicode sub/superscripts to explicit typesetting. Displayed equations are rendered with matplotlib. Final files are under `output/pdf/`; intermediate equations and previews are under ignored `research/paper/tmp/`. Custom fonts or different dependencies can change layout and file hashes, so rebuilt PDFs require visual inspection.

The [Chinese manifest](build-manifest.json) and [English manifest](build-manifest.en.json) record manuscript/PDF hashes and equation/figure/table counts. [Chinese visual QA](visual-qa.json) and [English visual QA](visual-qa.en.json) record inspected outputs. PDF byte identity across rebuilds or platforms is not promised.

## Review and reproducibility

- [Bilingual validation](translation-validation.json): identical displayed equations, matching numerical content in 116 table-body cells, section numbering, shared figure paths, reference URLs, and source/PDF hashes. The standard-library validator also runs in CI.
- [Translation review](translation-review.md): internal review of completeness, scientific terminology, numerical values, and limits of interpretation.
- [Scientific review](review.md): the earlier internal review of the Chinese source.
- [Scientific reproduction guide](../REPRODUCING_POPULATION.md): offline checks, public input retrieval and hashes, and stage-two commands.
- [Stage one](../population-stage1/README.md), [stage two](../population-stage2/README.md), and [literature comparison](../population-audit/literature-overlap.md): detailed calculation and provenance records.
- [Notices and licensing](../../NOTICE.md): original code and third-party scientific inputs are managed separately; raw inputs remain separate downloads.

Deterministic translation checks do not establish semantic equivalence by themselves. Internal AI-assisted review is not external peer review or novelty certification. Individual authors, affiliations, and funding have not been invented. Formal submission still requires author responsibility, event-side provenance and selection-consistency checks, and model validation appropriate to the intended scientific claims.
