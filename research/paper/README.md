# 公开双黑洞群体推断的条件稳健性

**研究技术稿 v0.2，2026-10-24，未同行评审。**

[English guide](README.en.md) · [内部科学审阅](review.md) · [双语审阅](translation-review.md)

| 版本 | 阅读 PDF | 可编辑全文 |
|---|---|---|
| 中文版 | [14页 LaTeX PDF](../../output/pdf/blackhole-population-robustness-latex-zh-v0.2.pdf) | [中文全文](manuscript.zh.md) |
| English edition | [14-page LaTeX PDF](../../output/pdf/blackhole-population-robustness-latex-en-v0.2.pdf) | [Full English manuscript](manuscript.en.md) |

英文版为全文对应译本，两版使用相同研究版本及图表；LaTeX 成品现均为14页。14个显示公式、6表共116个数据单元格、章节编号及9条参考文献的对应关系已核对，见[一致性检查](translation-validation.json)。

本文将两个科学阶段整理为详细论文：69个教程来源事件、81个固定假设设置、共同支持先验控制、41972条注入的逐行来源恢复，以及2048次共享样本的经验重采样。正文包含14个公式、6张表和3幅图，讨论与已有研究的联系、复现方法和未解决的问题。

核心结果是有限模型族中的条件稳健性诊断。基准自旋幅度假设下，两个倾角宽度的共同先验后验中位数相差6.734个百分点；相应经验数值扰动范围为5.557–7.600个百分点。低自旋压力案例的差为1.115个百分点，经验范围为−0.293–3.601个百分点。这些经验范围不是天体物理可信区间，现象学分量也不是形成渠道标签。事件原始PE配置和事件—注入选择一致性仍待核实。

## LaTeX 论文与编译

可编辑项目：[中文 XeLaTeX](latex-zh/main.tex) · [English pdfLaTeX](latex-en/main.tex)。每个目录含完整 `main.tex` 和三张矢量 PDF 图，可单独复制编译。采用 A4、11pt、25mm 页边距、原生公式、三线表、自动编号及超链接引用。中文使用 `ctexart` 和 TeX 发行版的 Fandol 字体；英文使用 `article` 和 Latin Modern。

参考 [arXiv 官方 TeX 源文件指南](https://info.arxiv.org/help/submit_tex.html)，不是某期刊专属模板，也没有提交 arXiv。文献嵌入 `thebibliography`，无需 BibTeX；日期固定。

安装含 ctex/Fandol、natbib、microtype、placeins 等宏包的 TeX Live 或 MiKTeX，将 `xelatex`、`pdflatex` 加入 PATH。在仓库根目录运行：

```sh
python -X utf8 research/paper/compile_latex.py
python -X utf8 research/paper/validate_latex.py
```

脚本各编译两轮，将 PDF 导出到 `output/pdf/`。若修改 Markdown 内容，先运行 `python -X utf8 research/paper/build_latex.py`；若修改冻结结果的绘图，先运行下面的 `build_figures.py`。手动编辑 `main.tex` 后不要再运行生成器，以免覆盖；可在相应目录直接运行 `xelatex main.tex` 或 `pdflatex main.tex` 两次。

[源码清单](latex-source-manifest.json)、[编译清单](latex-build-manifest.json)和[页面检查记录](latex-visual-qa.json)保存引擎、页数和哈希。CI 核查已发布源码、图件及 PDF 的完整性，不在 CI 中重编译 TeX。重编译可能改变 PDF 字节哈希，需重新检查版面并更新检查记录。

## 早期 ReportLab 导出（保留作历史版本）

此前的 [中文13页](../../output/pdf/blackhole-population-robustness-v0.2.pdf)和[英文16页](../../output/pdf/blackhole-population-robustness-en-v0.2.pdf)来自 ReportLab；它们不是 LaTeX 编译结果。以下命令仅复现这些早期导出。

### 从冻结结果构建早期导出

在仓库根目录运行；构图和排版不需要下载原始数据，也不会重新运行群体推断：

```sh
python -m pip install -r research/paper/requirements.txt
python -X utf8 research/paper/build_figures.py
python -X utf8 research/paper/build_pdf.py
python -X utf8 research/paper/build_pdf.py --language en
python -X utf8 research/paper/validate_bilingual.py
```

本次构建使用Python 3.14。默认字体为Windows的SimSun/SimHei；其他平台需指定可合法使用、含所需中文字形且有TrueType轮廓的字体，例如：

```sh
python -X utf8 research/paper/build_pdf.py --font /path/to/cjk-regular.ttf --bold-font /path/to/cjk-bold.ttf
```

英文版默认使用matplotlib附带的DejaVu Serif字体，无需Windows中文字体。支持相同的`--font`和`--bold-font`选项。不随仓库分发字体文件。PDF嵌入所用字体子集；构建器检查正文及粗体字形覆盖，将Unicode上下标转换为显式上下标排版，保留负号含义；显示公式由matplotlib渲染。输出位于`output/pdf/`，临时公式和页面预览位于忽略目录`research/paper/tmp/`。字体或依赖变化会改变版式和PDF哈希，应重新逐页检查。

[build-manifest.json](build-manifest.json)记录稿件与PDF哈希、公式/图表数量；[visual-qa.json](visual-qa.json)记录本次页面检查。PDF的构建时间和字体子集可能使文件字节变化，不承诺跨平台字节相同。

英文版对应记录为[build-manifest.en.json](build-manifest.en.json)及[visual-qa.en.json](visual-qa.en.json)。双语检查只依赖Python标准库，并已接入CI；它核查结构、数值和文件哈希，不能替代正文语义审阅。

## 科学计算复现与材料

- [完整复现指南](../REPRODUCING_POPULATION.md)：离线检查、输入下载与哈希、阶段二重算命令。
- [阶段一](../population-stage1/README.md)：物理积分量、信息不变性和共同支持。
- [阶段二](../population-stage2/README.md)：成对经验重采样、协方差、来源恢复及权重控制。
- [文献重合分析](../population-audit/literature-overlap.md)：既有研究与本次实现贡献的边界。
- [许可与来源](../../NOTICE.md)：原创实现和第三方数据分别管理；原始大文件单独下载。

文稿署名为项目研究稿，未代填个人作者、机构或基金。内部AI辅助审阅不等于外部同行评审，也不构成新颖性或发表资格认证。正式投稿前应补齐作者责任确认、事件侧来源和选择一致性核验，以及与拟投稿范围相称的模型检验。
