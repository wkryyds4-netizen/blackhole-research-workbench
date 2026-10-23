# SXS 真实数据接入与验证记录

更新：第二阶段已接入 SXS:BBH:2327 三分辨率、对齐敏感性面板，并成功运行 GWSurrogate。最新汇总见 [第二阶段记录](../stage2-results/README.md)。下文保留第一阶段历史过程；其中“尚未运行”的状态已被后续成功运行替代。

日期：2026-10-23。工作台入口：<http://127.0.0.1:8765/workbench.html?dataset=sxs-0003-lev5>。

本轮已完成 SXS 视界数据接入。GWSurrogate 只完成源码、模型下载和构建调查，尚未生成波形或动力学；两者的状态不能混同。

## 已运行的数据

| 模拟 | 固定来源 | 分辨率 | 显示帧数 | 用途 |
|---|---|---|---:|---|
| SXS:BBH:0001 | [Zenodo 2621619](https://doi.org/10.5281/zenodo.2621619) | Lev5 | 3774 | 近无自旋等质量旋近和共同视界 |
| SXS:BBH:0003 | [Zenodo 13167015](https://doi.org/10.5281/zenodo.13167015) | Lev5 | 4746 | 近等质量、第一颗黑洞横向自旋约 0.499 的进动演示 |
| SXS:BBH:0389 | [Zenodo 3274661](https://doi.org/10.5281/zenodo.3274661) | Lev3 / Lev2 | 2264 / 2275 | 同一近无自旋模拟的分辨率比较 |

0003 的元数据明确标注 deprecated / superseded_by SXS:BBH:2327。工作台保留此标记；本算例用于历史数据接口验证，不能作为新的精度基准。没有将 Zenodo 的记录版本当成 SXS 最新目录版本。

每份原始 HDF5 和元数据均先验证 Zenodo MD5，再记录 SHA-256。下载清单位于 `../sxs-data/<ID>/<Lev>/manifest.json`。完整记录保留作者、DOI、许可与文献键。数据按 CC-BY-4.0 使用，前端署名见 `../../simulator/public/licenses/SXS-DATA.txt`。

## 转换和界面

适配器直接用 h5py 读取官方文档规定的视界中心、Christodoulou 质量和无量纲自旋数组；没有安装完整 SXS Python 包。A/B 是两个分量的视界，C 是共同视界。界面支持本地数据选择、拖动时间、轨迹尾迹、自旋箭头、质量/自旋/分离度曲线、点击曲线定位、共同视界跳转和来源完整导出。

`Mref` 是元数据中参考时刻两颗黑洞质量之和。时间为 `(t_source−t_C_first)/Mref`，质量与位置除以 Mref。位置、自旋同时旋转 `(x,y,z)→(x,−z,y)`，这是行列式 +1 的旋转，不是镜像。没有去除坐标质心漂移。显示中的 65 M☉ 是任意尺度，不是 GW150914 的拟合质量。

每个量独立按原始时间网格线性插值，禁止外推。显示压缩以原始分段线性数据为参考，对每分量约束位置误差 ≤0.003 Mref、质量误差 ≤1e-5 Mref、自旋误差 ≤1e-4。`test_sxs_adapter.py` 对四份数据的所有保留帧反向变换，并对其覆盖区间内的所有原始采样点重新核对压缩误差。这里没有估计原始数值相对论模拟的真实误差。

0001、0003 的分量视界数据分别在共同视界首样本前约 0.086 Mref、0.118 Mref 结束。这些区间显示最后 A/B 帧并报告缺口。C 出现时离散切换，不把两个视界插值成一个。视界时刻没有与无穷远波形的迟滞时间同步，t=0 不是波形峰值。

原始微小自旋和方向分量全部保留。超出光学模型支持范围时自动使用坐标视图；支持单 Kerr 的数值帧也只是把瞬时参数代入定常参考图像，不是并合的时变数值光学。视界质量曲线不能直接解释为辐射能量曲线。

## 本轮研究结果

0389 的 Lev2/Lev3 比较使用相同原始坐标时间、统一 Lev3 参考质量，未拟合时间或相位偏移。区间为原始 t∈[640, 5003.0752]，取 2001 个公共时刻，终点在较早共同视界前 50 个原始时间单位。

- 坐标分离度差最大约 **0.12336 Mref**。
- 单个中心位置向量差最大约 **0.76570 Mref**。
- 自旋向量差最大约 **1.0037×10⁻⁴**。
- 共同视界首样本时刻差约 **4.62542** 个原始时间单位。

这些是两分辨率间的坐标诊断，不是收敛阶、严格误差上界、波形失配或观测误差。位置差与分离度差的不同提示比较对轨道相位和时间约定敏感；目前尚未把相位积累、坐标效应和离散误差分别量化。

对“通用型泛化”的具体启示是：需要同时保存原始时钟、显示时钟、旋转、缺口、模型有效范围和数据弃用状态。一个可以加载的 JSON 不足以保证不同模型可直接叠加比较。下一步应先接入未弃用的进动参考数据，再比较参考频率、旋转和时间对齐规则，并把这些选择作为敏感性分析的维度。

群体形成渠道研究没有从这些动画得到新的比例估计。现有 69 事件、81 个条件模型的研究仍在原工作台保留；下一步需要补齐事件后验波形版本及注入选择阈值，再比较波形假设如何改变联合后验和层次推断结果。

## GWSurrogate 运行调查

官方 1.2.0 源码包 SHA-256：`316a640576d32aaa45aaf9a6d876b240af356c89131e6af8257039b3818baae8`。

NRSur7dq4 数据来自 [Zenodo 3348115](https://doi.org/10.5281/zenodo.3348115)，17,479,528 字节，官方 MD5 `8e033ba4e4da1534b3738ae51549fb98` 已匹配，SHA-256 和 CC-BY-4.0 记录见 `../surrogate-data/manifest.json`。本轮没有使用目录中的其他模型，也没有把模型文件存在当作模型运行成功。

独立 Python 3.12 环境中的安装遇到三个可区分的问题：Windows GBK 读取 README 失败；GCC 参数及 `__restrict__` 与 MSVC 不兼容；进一步构建遇到 C99 变长数组与复数表达不受 MSVC 支持。UTF-8 与编译选项的最小补丁已保存在 `../adapter-sources/gwsurrogate-windows.patch`，但**该补丁不足以安装成功**。完整日志保存了编译诊断。未改写数值公式，也未用自制近似冒充代理模型输出。

后续采用兼容编译环境后，应先验证自然三次样条、四元数归一化、旋转方向和固定输入输出，再接动力学视图。NRSur7dq4 的 q≤4、自旋幅度≤0.8 是模型训练范围；参考频率、时长与外推需逐次检查。其返回的自旋、轨道相位和旋转四元数不包含 SXS 视界中心轨迹，不能据此声称得到了真实数值分离度。

## 复现与检查

在项目根目录运行（PowerShell；先建独立 Python 3.12 环境）：

```powershell
simulator/.adapters-venv/Scripts/python.exe -m pip install -r research/adapter-requirements.txt
simulator/.adapters-venv/Scripts/python.exe -X utf8 research/fetch_sxs.py
simulator/.adapters-venv/Scripts/python.exe -X utf8 research/sxs_adapter.py
simulator/.adapters-venv/Scripts/python.exe -X utf8 research/test_sxs_adapter.py
node simulator/test-generalization.mjs
node simulator/qa-datasets.cjs
node simulator/qa-workbench.cjs
```

下载脚本使用已保存的固定记录，不搜索“最新版”。后续转换和科学核对可离线运行。浏览器检查需要现有 8765 本地服务和脚本注明的 Playwright/Edge 环境。

`sxs-validation.json` 给出误差、来源哈希和分辨率比较；`browser-validation.json` 给出 13 项真实数据界面检查；原工作台的 27 项检查也已通过。截图保存在本目录。这里的检查不构成完整数值相对论、光学或代理模型的物理验证。
