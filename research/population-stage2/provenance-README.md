# 注入来源阶段成果：41,972 行、检出掩码与密度转换的精确复现

日期：2026-10-24。本次保持原始 `data/injections.npz`、后验数据及旧模型不变。

本地教程注入数据已经与一个官方原始文件建立**逐行精确对应**，不再仅依靠计数推测。结果为 [provenance-injection-match.json](provenance-injection-match.json)，下载校验见 [provenance-download.json](provenance-download.json)。

## 原始文件与匹配证据

来源：[LVK O1+O2+O3 Search Sensitivity Estimates](https://zenodo.org/records/7890398)，文件 `o1+o2+o3_bbhpop_real+semianalytic-LIGO-T2100377-v2.hdf5`，63,149,928 bytes。官方 MD5 `c21337a25b7318fdfd471cce549a547e` 已验证；下载文件 SHA256 为 `81883c571a9fe0b30747c7bc9daca6d358fcd68ef375c7ace0b93aa9ea2e74fe`。文件存于 `source-data/`，其 `.gitignore` 排除 HDF 和下载临时文件，未上传或提交。

2021 年 [5636816](https://zenodo.org/records/5636816) 与 2023 年 [7890398](https://zenodo.org/records/7890398) 两个发布的全部 5 个 HDF 官方校验值相同。实际下载使用前者的 API；因此不存在因这两个记录版本不同而换用不同 HDF 内容的问题。

| 核对项 | 结果 |
|---|---|
| 原始 `/injections` 行数 | 289,620 |
| 本地 41,972 行对应的原始行号 | 全部唯一恢复，原始顺序一致 |
| 主质量、质量比、红移 | 每列逐元素完全相等 |
| 两个自旋幅度、两个倾角余弦 | 从 Cartesian spins 转换后每列逐元素完全相等 |
| 指定检出规则的原始行集合 | 与本地严格同集，对称差为 0 |
| 指定 Jacobian 转换后的抽样密度 | 与本地 `prior` 逐元素完全相等 |

匹配使用 `(mass1_source, redshift)` 两元精确键。原文件的主质量单列存在重复，故不能单靠质量唯一定位；两元键的 289,620 个组合全部唯一。随后独立核对全部七个物理坐标，结果均为 `np.array_equal=True`，并非容差意义上的接近。

## 可以精确复现的检出规则

对该 HDF 的全部行，以下规则恰好得到本地 41,972 行：

```python
found = any(ifar_search > 1 year for five searches)
found |= (name in ['o1', 'o2']) & (optimal_snr_net > 10)
```

五个搜索为 `cwb`、`gstlal`、`mbta`、`pycbc_bbh`、`pycbc_hyperbank`；路径分别为 `/injections/ifar_<search>`。采用 `far_<search> < 1/yr` 的等价表达也恢复同一集合。

**这是该有限文件上数据与选择规则的精确重建，不等于找回了当年的命令、软件版本，亦不证明阈值在数学上唯一。** 多个落于相邻统计量间隙内的阈值可能给出同一有限集合。

## 抽样密度与运行权重

本地 `prior` 严格等于：

```python
sampling_pdf * mass1_source * (2*pi*a_1**2) * (2*pi*a_2**2)
```

这对应从两个源参考系质量转到主质量与质量比，并将两组 Cartesian spins 转成幅度/倾角余弦、积分掉方位角的密度转换。它**没有**除以 `/injections/mixture_weight`：使用该除法的候选密度不再与本地 `prior` 相等。

[官方说明](https://zenodo.org/records/7890398/files/o1%2Bo2%2Bo3-sensitivity-estimates.md?download=1) 分别定义 `sampling_pdf` 与需在 Monte Carlo 求和中使用的 `mixture_weight`；其构造依据是 [Essick (2021), DOI 10.3847/2515-5172/ac2ba7](https://doi.org/10.3847/2515-5172/ac2ba7)。因此，在本地 `p_pop/prior` 上另乘 `mixture_weight` 是一次明确的权重控制实验，不会重复应用已经吸收的权重。

| 运行 | 原文件保留行数 | 本地检出行数 | 每行 mixture_weight（该运行内为常数） |
|---|---:|---:|---:|
| O1 | 51,194 | 1,623 | 1.0000057241333624 |
| O2 | 109,972 | 3,881 | 0.9999982192691211 |
| O3 | 128,454 | 36,468 | 0.9999992049793384 |

这里的“原文件保留行数”是 `n_accepted` 子集，不是每个运行生成的全部模拟次数，不能拿来充当分层 bootstrap 的原始 trial counts。

对任意非负总体权重，补入运行权重前后选择效率之比必落在 `min(w)` 与 `max(w)` 之间。因此，69 事件条件似然中这项修正对**对数似然形状的最大变化范围**有解析上界：

`69 * log(max(w)/min(w)) = 0.000517834611630186`。

这能限制本次漏乘近单位运行权重的影响，不能限制错误的事件先验、不同事件选择或有限样本 Monte Carlo 误差。

后续已单独执行 [mixture-weight-control.json](mixture-weight-control.json) 控制：乘回 `mixture_weight` 后，六个模型在共同 r 先验下的后验中位数最大改变 **6.6502e-9**（比例单位，即 **6.6502e-7 个百分点**）；原模型诱导先验下最大改变约 **3.65e-8**。4001/8001 网格对这项修正量的差为约 **1.88e-11**。将整数归一化改为原始非整数归一化，归一化条件后验的 TV 差约 **1.52e-14**，处于浮点数精度水平。该控制支持“本次近单位运行权重遗漏对所测六模型影响极小”，不支持泛化为所有数值问题已解决。

## 非整数归一化与 IID 限制

HDF 的 `total_generated` 为 **73,280,198.73534885**，本地值 **73,280,198** 等于整数截断；相对差约 `1.0035e-8`。原始 `analysis_time_s=63,129,577` 除以 365.25 天/年后，正好等于本地 **2.000455579638502 年**。

非整数 `total_generated` 应作为官方提供的有效归一化元数据保存。文件头及已核查资料没有恢复其详细构造过程，不能将其宣称为已证实的字面 IID 整数试验次数。对当前条件于事件数 N 的模型，这个参数无关的归一化变化仅添加似然常数，归一化后验严格不变；率推断需要单独处理。

现在已恢复运行标签和权重，可以开展按运行固定保留数的分层经验扰动，与全体保留行的 IID 经验扰动对照。但这仍是**固定已保留样本的数值敏感性诊断**；没有恢复各运行的完整生成设计、原始抽样计数和未保留事件，不能称为完整注入 campaign 的重复抽样覆盖率校准。

## 教程历史与尚未完成的证据链

[2024-06-15 固定 commit 的原始 notebook](https://github.com/ColmTalbot/gwpopulation/blob/2daa71e8c56362f2f64e2166a0055478557f239a/examples/GWTC3.ipynb) 和当前 notebook 都引用同样两个 Google Drive 数据对象。已检索的路径历史共 3 次提交；两份 notebook 均没有原 HDF 文件名或检出阈值生成命令。此处只说明本次实际查过的公共路径，不能排除作者私有或其他仓库存在生成脚本。记录见 [provenance-history.json](provenance-history.json)。

此次突破覆盖注入文件、行、有限样本检出规则及密度转换。逐事件 PE 波形/先验/身份链和与观测事件一致的完整选择政策尚待恢复，不能据此宣布整个群体分析已成为官方复现。

## 文件接口与复现

`provenance-injection-row-map.npz` 按原 `data/injections.npz` 顺序存储：

- `original_row_indices`：`/injections` 中的整数行号。
- `run`：Unicode `U2`，值为 `o1`、`o2`、`o3`。
- `mixture_weight`：原 HDF 对应行的浮点权重。

```powershell
.\simulator\.adapters-venv\Scripts\python.exe -X utf8 research/population-stage2/provenance_match_injections.py
.\simulator\.adapters-venv\Scripts\python.exe -X utf8 research/population-stage2/provenance_history_audit.py
```

第一个脚本离线验证逐行、掩码、转换和输入校验，遇到不相等直接失败；第二个复核固定 notebook 与官方发布清单。下载脚本为 `provenance_download_bbh.py`，采用受检验的 HTTP byte ranges，并在完成后验证官方 MD5。`provenance_remote_header.py` 记录前期仅用小块字节范围排查候选的证据；最终结论以完整文件校验及逐行匹配为准。
