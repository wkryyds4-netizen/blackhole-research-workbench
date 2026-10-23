# 群体输入来源复核：恢复 FAR 元数据，保留尚未解决的选择问题

日期：2026-10-24。结果见 [provenance-followup.json](provenance-followup.json)。本次不改变旧模型、后验数组、注入数组或已发布推断。

## 已确认的发现

现有 69 个入选事件中，10 个 FAR 缺失全部来自 GWTC-1 的十个双黑洞。原因可以直接从 `model.py::identify_events` 复现：按文件名排序读目录，并以 `commonName` 为字典键；后读入的 GWTC-2.1 同名记录覆盖 GWTC-1 记录，连同其 `far=null` 一起覆盖。后验参数更新记录和搜索显著性记录没有分开保存。

从原本就已下载的 GWTC-1 目录恢复的数值如下。单位直接取 API 的 `yr^-1`，不做单位猜测。2026-10-24 再取的 [GWOSC 官方 JSON](https://gwosc.org/eventapi/json/GWTC-1-confident/) 与本地文件 SHA256 完全相同。

| 事件 | GWTC-1 API FAR / yr⁻¹ |
|---|---:|
| GW150914 | 1e-7 |
| GW151012 | 0.00792 |
| GW151226 | 1e-7 |
| GW170104 | 1e-7 |
| GW170608 | 1e-7 |
| GW170729 | 0.02 |
| GW170809 | 1e-7 |
| GW170814 | 1e-7 |
| GW170818 | 4.2e-5 |
| GW170823 | 1e-7 |

这些数字按目录原值保存，不将很小的 FAR 解读为精确测量；API 的上下界字段为 null，未提供足够信息区分显示值和搜索背景给出的限值。每条恢复记录保留目录版本、GPS、URL、原单位和被覆盖的记录。**这解决的是跨目录 FAR 元数据缺失；后验数组事件身份仍来自摘要匹配。** 十项都满足当前 `<1/yr` 条件，因此仅恢复 FAR 不会改变原 69 个成员。

原 [GWPopulation 教程](https://colmtalbot.github.io/gwpopulation/examples/GWTC3.html) 对 73 项依次执行 `del posteriors[15]`、`del posteriors[38]`。第二次删除作用于已缩短的列表，故对应原始索引 15 和 39，按当前匹配为 GW190425 和 GW190814，余 71 项。当前模型还排除：

| 原索引 | 当前匹配名 | 当前目录 FAR / yr⁻¹ |
|---|---|---:|
| 53 | GW191204_110529 | 3.3 |
| 64 | GW200208_222617 | 4.8 |

因此，“69 个事件”这个计数相同本身不能证明我们的事件集、波形与官方 69 个 BBH 群体样本相同。教程本身说明事件选择只可作定性使用。

## 注入与先验：恢复了接口文档，没有恢复原运行配置

[教程指向的注入发布](https://zenodo.org/records/7890398) 同时包括 O1/O2 半解析 SNR 和 O3 实际搜索统计量；发布说明要求分别指定检出阈值。它提供若干子群及 mixture 文件，只有一个 release 链接不足以确定原始 HDF 文件。

本地 NPZ 有 41,972 条保留注入，`total_generated=73,280,198`，`analysis_time=2.000455579638502`。其中没有运行标签、检出统计量、原行号、阈值、原 HDF 文件名或预处理版本。通过 `pickletools.genops` 读取两个原 pickle 的文本操作码，不执行反序列化；也没有发现隐藏的 GW 事件名、波形名或检出阈值字段。

本次归档的[当前官方 `vt_helper` 文档源码](https://ratesandpopulations.docs.ligo.org/gwpopulation_pipe/_modules/gwpopulation_pipe/vt_helper.html) 的函数默认值是 IFAR > 1 yr、O1/O2 网络 SNR > 10。该代码把原始 Cartesian-spin 抽样密度转为质量比和自旋幅度/倾角密度时，使用 `mass_1 * (2*pi*a_1²) * (2*pi*a_2²)` 因子。**这些是核验原文件时的参考公式与默认参数，不是原教程运行采用这些选项的证据。** 不应根据数组长度或当前默认值填写历史阈值。

[GWTC-2.1 PE 发布](https://zenodo.org/records/6513631) 有 cosmo/nocosmo 文件及多个波形标签、Mixed 标签；本地简化数组仅保留数字 `prior`，缺少逐事件先验配置、波形和版本映射。当前 [`data_collection` 文档](https://ratesandpopulations.docs.ligo.org/gwpopulation_pipe/_modules/gwpopulation_pipe/data_collection.html) 支持多种质量与距离先验，进一步说明仅凭数字列名无法确认其物理测度。

## 复现与下一步

在仓库根目录执行：

```powershell
.\simulator\.adapters-venv\Scripts\python.exe -X utf8 research/population-audit/provenance_followup.py
```

该步骤完全离线，使用 JSON、`np.load(..., allow_pickle=False)` 与静态 pickle 操作码。它断言 10 个覆盖案例、单位、成员数，以及两次顺序删除对应的原始索引。输入 SHA256 均写入结果。可选的 `collect_provenance_sources.py` 只获取小文档；`provenance-sources/manifest.json` 记录下载 URL、日期、字节数和校验值。Zenodo 说明文件的 shell 下载遇到 TLS EOF，已记录失败；其内容通过 web 工具读到，未伪造本地快照。

要建立可作为正式推断起点的基线，还必须从原始 PE 文件/预处理清单恢复逐事件身份、波形和先验，并从原始注入 HDF 重建同一检出规则。下一版应把 `pe_record` 与 `search_record` 分开存储，禁止空值覆盖非空搜索信息；FAR 缺失也应明确标记，不能自动通过筛选。上述修复没有证明注入选择与事件选择一致，不能据此升级形成渠道结论。
