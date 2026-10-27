# 五渠道条件推断：24组完整对照

2026-10-27。全部69事件、11,180,578后验行，完整事件与选择积分、两套后验算法及共同正则项拆分均已完成。结论限定于声明的限制域、固定出生自旋0.1/CE效率2及额外自旋补全；不是宇宙形成比例测量。

| 比较轴 | 最大中位数移动 / 百分点 | 渠道 |
|---|---:|---|
| 正则化概率 | 8.654 | CE |
| 隐藏自旋补全 | 3.415 | CE |
| O3银行替换 | 0.614 | CE |
| 搜索口径 | 0.427 | NSC |

最大值分别在各自网格点取得，不能相加或当作可信区间。两套算法的分位数最大差为0.169个百分点；四次重要性重复最大跨度为0.455个百分点。归一化抽样半部替换的最大移动为0.051个百分点。

官方混合银行、匹配滤波、均匀基准条件下，CE内禀权重中位数从83.85%（无正则化）到76.17%（1%共同项）。但共同项的计算效率约8.06%，远高于未正则化CE约0.139%；它在检出模型混合中的责任份额中位数达到12.14%。这说明“小密度正则项”不必是“小选择修正”，不表示测得12.14%的未知物理渠道。

`branching-results.json`的`detected_quantiles`是**正则化渠道标签**。共同均匀项按beta分配给标签；`detected-component-results.json`另给出模型项、共同项和去共同项后的条件检出比例。边缘中位数不必加和为1。

## 限制与下一步

最小事件贡献集中量只有5.24，最小选择集中量218.21。这些值不是PE链或活动ESS，也未用于删除事件。后验求值稳定不代表银行尾部已经校准。历史注入参考/筛选、PE未见支持、波形误差和真实渠道分量自旋仍未解决。两个条件补全只是压力模型，并非已验证的形成物理。

下一项实质工作是对本批渠道积分做成对有限银行传播，不能复用旧g/r模型的抽稀范围。随后再检验渠道特定自旋补全和超参数变化。

## 图表、论文与证据

- [Posterior figure](channel-posteriors.svg)
- [Sensitivity and regularizer](channel-sensitivity.svg)
- [中文 LaTeX PDF](../../output/pdf/blackhole-channel-completion-latex-zh-2026-10-27.pdf) / [English LaTeX PDF](../../output/pdf/blackhole-channel-completion-latex-en-2026-10-27.pdf)
- [LaTeX projects and shared facts](../paper/channels-2026-10-27/)
- [Method and literature comparison](measure-method.md)
- [Complete synthesis](branching-synthesis.json)
- [Independent posterior validation](observed-mcmc-validation.json)


## 复现

先恢复来源模型及此前已验证的69事件和注入缓存；环境与来源见本目录README。在仓库根目录运行以下命令。缓存协议与源码哈希绑定，重跑校准/改动设计后不要混用旧输出；先保留当前版本。不运行根目录model.py。原直接乘积失败与短链失败均保留。

```text
python research/population-channels/validate_fast_measure.py
python research/population-channels/validate_channel_weights.py
python research/population-channels/run_channel_integrals.py
python research/population-channels/validate_channel_products.py
python research/population-channels/validate_branching_sampler.py
python research/population-channels/fit_branching.py
python research/population-channels/validate_independent_mcmc.py --refined
python research/population-channels/validate_observed_mcmc.py
python research/population-channels/derive_channel_components.py
python research/population-channels/summarize_branching.py
python research/paper/channels-2026-10-27/build_note.py --render
```
