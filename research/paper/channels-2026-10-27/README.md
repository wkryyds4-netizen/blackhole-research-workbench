# Conditional five-channel research note / 五渠道条件研究稿

This interim note accompanies the [24-case results](../../population-channels/RESULTS.en.md) ([中文](../../population-channels/RESULTS.md)). It reports fixed-context, restricted-domain, saved-bank inference, not measured cosmic formation fractions. Numerical posterior agreement does not certify finite-bank errors or missing support. The overall research project is paused at the user's request.

- [English PDF](../../../output/pdf/blackhole-channel-completion-latex-en-2026-10-27.pdf), 7 pages, pdfLaTeX.
- [中文 PDF](../../../output/pdf/blackhole-channel-completion-latex-zh-2026-10-27.pdf), 7 页，XeLaTeX。
- `latex-en/main.tex` and `latex-zh/main.tex` are standalone projects with their own supporting files and figures.
- `build-manifest.json` records actual compilation diagnostics, page counts and file hashes; `evidence-manifest.json` connects shared numbers and figures to research outputs. These records do not assert independent peer review.

Run `python research/paper/channels-2026-10-27/verify_checkpoint.py` from the repository root for an offline integrity check. It does not rerun research or TeX compilation. To rebuild with installed TeX engines, use `build_note.py`; `--render` additionally requires Poppler. Rebuilding can change PDF bytes and requires renewed inspection. Original upstream products and local numerical caches are retrieved separately, as documented in [the reproduction guide](../../REPRODUCING_POPULATION.md).

两版共享数值、公式结构、表格和矢量图。署名沿用澄衡（Chengheng）；文稿保留研究笔名说明与未同行评审状态。既有 v0.3 论文保持原快照。
