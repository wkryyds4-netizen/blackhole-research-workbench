"""Bilingual study notes with identical result-derived tables."""
from pathlib import Path
import hashlib
import json
import numpy as np

BASE=Path(__file__).resolve().parent


def table(headers,rows):
    esc=lambda v:str(v).replace('|',r'\|')
    return '| '+' | '.join(map(esc,headers))+' |\n|'+'|'.join(['---']*len(headers))+'|\n'+''.join('| '+' | '.join(map(esc,row))+' |\n' for row in rows)


def main():
    d=json.loads((BASE/'results.json').read_text());v=json.loads((BASE/'validation.json').read_text())
    assert v['status']=='pass' and v['result_sha256']==hashlib.sha256((BASE/'results.json').read_bytes()).hexdigest()
    summaries=[]
    for name,case in d['cases'].items():
        q=np.array(case['quantiles'])
        summaries.append([name,f"{case['smallest_magnitude_contrast_pp']:.6f}",f"{case['largest_magnitude_contrast_pp']:.6f}",
                          f'{100*q[:,1].min():.4f}–{100*q[:,1].max():.4f}%',f"{case['median_sensitivity_pp']['joint_mass_redshift_maximum_span']:.6f}"])
    summary=table(['Bank:policy','Min amplitude contrast (pp)','Max amplitude contrast (pp)','All-model median range','Max mass/redshift median span (pp)'],summaries)
    case=d['cases']['official_mixture:matched_filter'];contrasts=case['magnitude_contrasts']
    rows=[];support=[]
    for mi,mass in enumerate(['alpha2_peak10','alpha3_peak10','alpha2_nopeak']):
        for zi,kappa in enumerate([0,2,4]):
            values=[next(c['value_pp'] for c in contrasts if c['mass_id']==mass and c['kappa']==kappa and c['shape']==shape) for shape in [0,4,16]]
            rows.append([mass,kappa,*[f'{x:.6f}' for x in values]])
            record=d['support']['models'][mi*18+zi*6]
            pure=record['banks']['pure_bbh']['uncovered_target_probability'];mixed=record['banks']['official_mixture']['uncovered_target_probability']
            support.append([mass,kappa,f'{100*pure:.6f}%',f'{100*mixed:.6f}%',f'{100*(pure-mixed):.6f}'])
    contrast_table=table(['Mass model','kappa','shape0 contrast (pp)','shape4 contrast (pp)','shape16 contrast (pp)'],rows)
    support_table=table(['Mass model','kappa','Pure BBH uncovered','Mixture uncovered','Added coverage (pp)'],support)
    shifts=table(['O3 policy','Largest bank median shift (pp)','Largest posterior TV'],[[policy,f"{record['maximum_median_shift_pp']:.6f}",f"{record['maximum_tv']:.8f}"] for policy,record in d['bank_comparisons'].items()])
    reproduction='''```text
python research/population-mass-redshift/run_grid.py
python research/population-mass-redshift/validate_grid.py
python research/population-mass-redshift/plot_grid.py
python research/population-mass-redshift/summarize_grid.py
```
'''
    links='''[Before-fit design](design.json) · [All results](results.json) · [Independent validation](validation.json) · [Vector PDF](mass-redshift-comparison.pdf)

![54-model mass/redshift comparison](mass-redshift-comparison.svg)
'''
    zh='''# 全69来源：质量与红移假设下的54个有效自旋控制模型

2026-10-27。本地研究阶段；没有覆盖已公开的v0.3论文、原始模型或工作台数据。本轮继[官方O3混合注入对照](../population-mixture/README.md)，在计算新权重前保存了固定设计。

**结果：在指定的可分离模型族中，质量/红移假设可让g的后验中位数移动约3.95个百分点，但两种自旋幅度族之间的同形状中位数差仍为正（10.36–17.14个百分点）。这一排序仅限本次模型与保存银行，不能升级为形成渠道识别或模型无关的g测量。**

## 实际计算了什么

保留全部69事件、11,180,578条原始先验表示后验行，以g=P(χ_eff<0)为双星总体量，先验Uniform[0,1]，条件于观测事件数N。模型组合为：

- 质量：α=2+10%高斯峰、α=3+10%高斯峰、α=2无峰；峰为35±5 M☉，范围3≤m₁≤100、m₂≥3 M☉，质量比条件密度仍正比于q。
- 红移：κ=0/2/4，p(z)正比于dVc/dz×(1+z)^(κ−1)，z≤2.3；固定同一Planck15网格及源坐标。
- 自旋幅度：Beta(1.5,8.5)与Beta(2,5)；对称偶形状λ=0/4/16，与既有g控制相同。
- 注入：纯BBH与官方四总体O3混合银行；分别用匹配滤波并集与含cWB并集，保留同一早期O1/O2项和时间/生成数权重。

共有54个模型、4个银行/口径组合、216组条件后验。没有拟合α或κ、改变事件集合、依ESS筛掉事件，或给模型分配后验概率。g仍被假定不随质量/红移变化；2PN守恒及角度测度保持仍是近似条件。

## 同形状幅度模型差与绝对模型依赖

以下差均为同一质量/红移/偶形状下Beta(2,5)减Beta(1.5,8.5)的后验中位数，单位为百分点。108个成对差全部为正；它们不是来自同一联合参数模型的“差值后验”，没有统计显著性的解释。

'''+summary+'''
最后一列是在固定幅度族/偶形状下，九个质量—红移选择之间的最大中位数跨度；它是声明的模型范围，不是可信区间或模型平均误差。

官方混合、匹配滤波口径的全部27个差：

'''+contrast_table+'''
在这一主口径下，单独质量选择的最大中位数跨度为2.551080 pp，单独κ选择为3.115004 pp，同时变化为3.952021 pp。替换注入银行的同模型变化较小，但不能据此认证真实选择误差更小。

'''+shifts+'''
TV是同一g量上归一化后验的总变差距离，范围0–1。不同银行共享真实注入行；这些差不具有独立重复实验的误差条解释。

## 模型也改变了未覆盖支持

下面是源总体在O3提议支持外的概率质量，不是探测器检出比例。偶形状不改变该支持质量；两个幅度族在表中精度下相同，完整54项数值保存在results.json。

'''+support_table+'''
κ从0到4使官方混合仍未覆盖的目标概率从约20.24%升至约46.31%。仅用p_det≤1得到的遗漏选择贡献上界可达到已表示选择分量的约30–1702倍，因而不能从这个平凡上界认证遗漏很小。实际遗漏效率未被测量。

固定z网格的梯形归一化与更精细分段积分的相对差最大8.29e-8；支持表使用连续插值积分重新归一化的目标概率。后验保留旧数值协议，整体目标归一化因子在条件于N的似然中抵消。

## 数值验证与仍未解决的部分

独立脚本重新读取所有69事件的原始缓存，用另一套平直ΛCDM积分恢复PE先验，再从Cartesian自旋及54个目标密度计算积分；选择项直接读取原HDF，后验使用69因子的乘积，没有调用生产脚本的权重或后验函数。最大事件积分相对差2.17e-9，选择相对差2.34e-15，相对峰值密度差1.13e-10，分位数绝对差1.40e-11，支持概率差6.67e-16；均通过计算前记录的门槛。目标宇宙学网格有意复用，不能把实现核验叫作独立数据验证。

最小事件贡献ESS约107.9，最小汇总选择贡献集中度ESS约1508.3；没有据此声称完整MC精度或活动覆盖率。GW191204_171526在所有54模型中仍有零负分量的保存银行积分，原先的缺失分量问题没有因重加权而消失。本轮没有给全部54模型重做原先512次成对抽稀，不沿用六模型抽稀范围充当新模型误差条。

## 与文献的关系和下一步

参数化自旋分布的形状偏差及灵活模型导致更大不确定度已有[Golomb与Talbot](https://arxiv.org/abs/2210.12287v2)研究；质量比—有效自旋相关已有[Adamcewicz等](https://arxiv.org/abs/2307.15278v2)研究（本轮核查其摘要及版本信息）。本项目不把“改变质量模型会影响自旋结果”当新发现。这里的增量是来源完整的69事件、明确先验、两银行/两口径、54个固定模型以及独立数值核验的一致证据链。

尤其需要继续测试g对质量/红移的依赖或更一般参数相关性；本次可分离网格没有测试这些假设。共同自旋参考、支持外事件似然与检出效率、真实注入执行协议、波形系统误差和物理渠道预测仍未完成。论文后续应按此边界修订，不能把所有正中位数差写成渠道证据。

## 复现

先恢复前序全69事件缓存及官方O3混合数据，步骤见[总体复现说明](../REPRODUCING_POPULATION.md)。依赖固定版本在仓库population-requirements.txt，图另需Matplotlib；运行版本在results.json。原始后验、HDF和NPZ检查点保持本地/忽略；event-integrals.json仅包含模型积分及来源哈希。修改密度协议时须对缓存版本化，不能无声复用。

'''+reproduction+links
    en='''# All69-source controls:54 mass/redshift models for negative effective spin

27 October2026. Local research stage following the[official O3 mixture comparison](../population-mixture/README.en.md). The new design was recorded before new weights or fits. Published v0.3 papers, frozen models and workbench data remain unchanged.

**Finding: within the specified separable families, mass/redshift choices move the g posterior median by up to3.95 percentage points, while the same-shape contrast between the two magnitude families remains positive(10.36–17.14 pp). This is conditional ordering on the saved banks, not channel identification or a model-independent measurement of g.**

## Fixed inputs and tested models

All69 events and11,180,578 original-prior posterior rows are retained. The binary-level functional is g=P(χ_eff<0), with Uniform[0,1] prior and inference conditional on N. The grid comprises:

- Three primary-mass choices:α=2+10% Gaussian peak,α=3+10% peak, andα=2 with no peak. The peak is35±5 solar masses;3≤m₁≤100 and m₂≥3, with conditional mass-ratio density proportional to q.
- Redshift exponentsκ=0/2/4 in p(z) proportional to dVc/dz×(1+z)^(κ−1), z≤2.3, using the same Planck15 grid and source coordinates.
- Magnitudes Beta(1.5,8.5) and Beta(2,5), with even-spin shapesλ=0/4/16 as in the previous control.
- Two related O3 banks(pure BBH and official four-population mixture), each under matched-filter union or union with cWB. The early O1/O2 prescription and time/trial weights are fixed.

This yields54 models and216 conditional fits. Noα/κ hyperparameter fit, event/ESS-based rejection, or model-probability weighting is performed. g is still mass/redshift independent.2PN invariance and preserved angular measure remain approximate assumptions.

## Contrasts and absolute model dependence

Each contrast subtracts the Beta(1.5,8.5) posterior median from Beta(2,5) at identical mass,redshift and even shape. All108 contrasts are positive. They are differences between separate conditional model summaries, not a difference-parameter posterior or significance statement.

'''+summary+'''
The last column is the largest span across nine mass/redshift medians at fixed magnitude/shape. It is a declared model range, not a credible interval or model-average error.

All27 contrasts for the official-mixture/matched-filter case:

'''+contrast_table+'''
For that main case, changing mass prescription alone produces a maximum median span of2.551080 pp; changingκ alone gives3.115004 pp; varying both gives3.952021 pp. Bank-replacement changes at identical settings are smaller, but this does not bound actual omitted sensitivity.

'''+shifts+'''
TV is normalized posterior total variation on the same g functional, on0–1. Shared injection rows prevent interpreting bank differences as independent-experiment error bars.

## Model-dependent missing support

The table describes target-source probability outside the O3 proposal, not a detection fraction. Even shapes do not change this support probability. Magnitude families agree at the printed precision; all54 records are retained in results.json.

'''+support_table+'''
Acrossκ=0 to4, the official mixture still misses about20.24%–46.31% of target probability. The trivial p_det≤1 omitted-selection upper bound is about30–1702 times a represented component, so it does not certify negligible omissions. Actual omitted efficiency is unmeasured.

The maximum difference between the fixed-grid trapezoid normalization and a refined piecewise integral is8.29e-8. Support probabilities use the continuously interpolated target with its integral renormalized; posteriors preserve the old numerical protocol. The global target normalization cancels in the likelihood conditional on N.

## Validation and limits

The independent validator rereads all69 raw posterior caches, reconstructs the PE prior with a separate flat-LambdaCDM integral, forms Cartesian-spin densities for all54 targets, sums original HDF injections and multiplies all69 event factors. It does not import the producer weights/posterior. Maximum relative event-integral difference is2.17e-9; relative selection2.34e-15; peak-scaled density1.13e-10; absolute quantile1.40e-11; support probability6.67e-16. All pass the gates recorded before calculation. The fixed target cosmology grid is intentionally shared; this is implementation verification, not independent observations.

Minimum event weight ESS is about107.9 and minimum pooled selection concentration ESS about1508.3. Neither certifies full MC error or actual-campaign coverage. GW191204_171526 still has zero saved negative-component integrals under all54 models; reweighting does not recover unseen components. This study has not repeated the prior512 paired thinnings across all54 models, and does not repurpose those six-model ranges as new-model error bars.

## Prior work and next evidence

Parametric-spin shape bias and enlarged uncertainty with flexible models are already addressed by [Golomb & Talbot](https://arxiv.org/abs/2210.12287v2). Mass-ratio/effective-spin covariance is studied by [Adamcewicz et al.](https://arxiv.org/abs/2307.15278v2) (abstract/version checked here). This does not claim that mass-model sensitivity is new. The increment is the source-complete 69-event, explicit-prior, two-bank/two-policy, 54-model comparison with independent numerical checks.

Mass/redshift-dependent g or more general parameter correlations remain to be tested: this factorized grid cannot settle them. Common spin reference, missing event likelihood, omitted detection probability, actual historical injection protocol, waveform systematics and physical-channel predictions remain unfinished. Positive median contrasts must not be promoted into channel evidence.

## Reproduction

First restore the full69 event caches and official O3 mixture data via the[reproduction guide](../REPRODUCING_POPULATION.md). Scientific versions are recorded in results.json and pinned in population-requirements.txt; plotting also needs Matplotlib. Raw posteriors,HDF and NPZ checkpoints remain local/ignored. event-integrals.json holds only derived integrals/source hashes. Version the cache when its density protocol changes.

'''+reproduction+links
    for old,new in [('All69','All 69'),('all69','all 69'),('full69','full 69'),('October2026','October 2026')]:en=en.replace(old,new)
    (BASE/'README.md').write_text(zh,encoding='utf-8',newline='\n')
    (BASE/'README.en.md').write_text(en,encoding='utf-8',newline='\n')
    print('Bilingual notes generated with identical data tables.')


if __name__=='__main__':main()
