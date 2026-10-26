"""Generate bilingual stage notes from the same independently checked values."""
from pathlib import Path
import hashlib
import json
import numpy as np

BASE=Path(__file__).resolve().parent


def table(headers,rows):
    escape=lambda value:str(value).replace('|',r'\|')
    return '| '+' | '.join(map(escape,headers))+' |\n|'+ '|'.join(['---']*len(headers))+'|\n'+''.join('| '+' | '.join(map(escape,row))+' |\n' for row in rows)


def main():
    d=json.loads((BASE/'results.json').read_text());v=json.loads((BASE/'validation.json').read_text())
    p=json.loads((BASE/'provenance.json').read_text())
    assert v['status']=='pass' and v['result_sha256']==hashlib.sha256((BASE/'results.json').read_bytes()).hexdigest()
    summaries=[];contrasts=[]
    for metric in ['g','r']:
        for policy in ['matched_filter','with_cwb']:
            cases=d['metrics'][metric]['scenarios'][policy];new=cases['mixture_full_union'];old=cases['pure_bbh_baseline']
            summaries.append([metric,policy,f"{np.max(np.abs(new['median_shift_from_pure_bbh_pp'])):.6f}",
                              f"{np.max(np.abs(new['quantile_shift_from_pure_bbh_pp'])):.6f}",f"{max(new['posterior_total_variation_from_pure_bbh']):.8f}"])
            levels=[0,4,16] if metric=='g' else ['Beta(1.5,8.5)','Beta(2,5)']
            for i,level in enumerate(levels):contrasts.append([metric,policy,level,f"{old['model_median_contrasts_pp'][i]:+.6f}",f"{new['model_median_contrasts_pp'][i]:+.6f}"])
    common_table=table(['Metric','O3 policy','Max |median shift| (pp)','Max |5/50/95 shift| (pp)','Max posterior TV'],summaries)
    contrast_table=table(['Metric','O3 policy','Shape / magnitude','Pure BBH (pp)','Mixture (pp)'],contrasts)
    reproducible='''```text
python research/population-mixture/fetch_mixture.py
python research/population-mixture/inspect_mixture.py
python research/population-mixture/study_mixture.py
python research/population-mixture/validate_mixture.py
python research/population-mixture/plot_mixture.py
python research/population-mixture/summarize_mixture.py
```
'''
    links='''[Results](results.json) · [Independent calculation](validation.json) · [Source audit](provenance.json) · [Before-comparison design](design.md) · [Vector PDF](mixture-comparison.pdf)

![Official O3 mixture comparison](mixture-comparison.svg)
'''
    zh='''# 官方 O3 混合注入：全69来源的条件对照

2026-10-27。目标恢复运行后新增的本地阶段研究；已公开的 v0.3 论文与结果保持不变，本结果尚未并入该论文或工作台。

我们从[固定官方发布7890437](https://zenodo.org/records/7890437)读取推荐的四类总体混合注入，替换原来纯BBH的O3选择项。全69事件的原始先验积分、两个六模型族、早期O1/O2项及观测时间比例保持固定。数据中共有1,141,491条接受注入，总生成数86,318,296；284,672条原BBH行全部找到唯一精确连接，40个共用物理/搜索字段逐项相同。

**阶段结论：对这次官方注入表示替换，g 后验变化小于所测试的幅度模型差；原生r的小排序仍依赖搜索口径。混合文件并未解决目标支持不完整。** 这是相关固定银行上的比较，不是独立敏感度测量或真实MC误差校准。

## 后验与模型对比

pp 为百分点。TV 是归一化后验之间的总变差距离（0–1），没有把不同科学指标相减。下表范围只覆盖保存的5/50/95分位数，不表示所有可能分位数的上界。

'''+common_table+'''
g 对比为相同 λ 下 Beta(2,5) 减 Beta(1.5,8.5) 的后验中位数；r 对比为同一幅度族下 σ=0.3 减 σ=1.2。g=P(χ_eff<0)，r=P(cos θ₁<0)，二者均不是已识别的形成渠道比例。

'''+contrast_table+'''
两种搜索口径下，混合文件全部新支持区域的已保存检出贡献均为零，所以完整混合与限制在原BBH提议支持内的选择结果相同。这个数值恒等不等于新区域的真实检出概率为零。原BBH行与混合行并非独立，不能把两估计之差当统计显著性。

## 新覆盖很小，未覆盖仍然存在

在固定目标 m₁≤100、m₂≥3 M☉ 下，BNS/NSBH提议不能增加支持；IMBH提议只在 m₁≥90、m₂≥10 M☉ 的高质量区域补上1.9<z≤2.3。这个质量区域占给定目标质量分布0.306136%，因此总共只增加0.102273个百分点的总体概率质量。混合提议仍未覆盖33.305257%的目标概率质量；这里是源总体质量，不是检出比例或选择函数误差。

新增目标支持中保留了178条注入，匹配滤波及含cWB并集均没有保存的检出行。它们来自有权重的提议和hopeless筛选，不能直接套用“178次独立目标抽样零检出”的二项上界。

后续来源诊断还发现：在只有BBH提议有支持的260,145条重叠行上，官方混合联合密度不严格等于独立BBH联合密度乘生成数比例；相对残差最大6.4486%，中位比值约0.9999995。可能的数值构建原因尚未核实，本研究保留各文件声明的联合密度，没有硬设两者相等。这是事后来源诊断，不是恢复了原始生成配置。

## 验证、文献位置与复现

独立脚本从HDF直接重建Cartesian目标/联合提议密度比、早期与O3积分、69事件完整乘积后验；没有调用生产脚本的权重或后验函数。最大选择相对差8.89e-15、相对峰值密度差1.94e-12、分位数差4.17e-15。支持质量另以分段多项式/Gauss–Legendre积分核对，绝对差4.14e-10。这验证数值实现，复用的事件积分仍是原来的输入，不是独立观测。

使用混合注入和混合权重遵循官方说明及[Essick 2021](https://doi.org/10.3847/2515-5172/ac2ba7)，不作为新方法。Hopeless cut、真实噪声和半解析近似的限制参照[Essick 2023](https://arxiv.org/html/2307.02765v1)。本轮增量是对本项目全部69来源/12个模型的具体对照、精确行连接、后验距离及剩余支持量化。

需要固定原始数据与已有全69事件缓存，依赖仓库 population-requirements.txt；图另需Matplotlib。原HDF、行映射和数值检查点不随Git分发，下载清单保留官方MD5与本地SHA256。仅下载本轮文件不能代替恢复前序研究输入。

'''+reproducible+links
    en='''# Official O3 mixture: conditional all-69-source comparison

27 October 2026. A new local research stage after the goal resumed. Published v0.3 artifacts remain unchanged; this follow-up is not yet incorporated into that paper or the workbench.

The recommended four-population mixture in [fixed official release7890437](https://zenodo.org/records/7890437) replaces only the pure-BBH O3 selection term. All69 original-prior event integrals, both six-model families, the early O1/O2 term and time fraction remain fixed. The mixture contains1,141,491 accepted rows from86,318,296 generated trials. Every one of the284,672 pure-BBH rows has a unique exact match, with40 common physical/search fields identical.

**Stage finding: the change of official injection representation moves g posteriors much less than the tested magnitude-model contrast; the tiny native-r ordering remains search-policy dependent. The mixture does not close the target-support gap.** This compares related fixed banks, not independent sensitivity measurements or calibrated campaign errors.

## Posterior and model comparisons

pp denotes percentage points. TV is normalized posterior total variation on0–1. Different scientific functionals are never subtracted. The quantile-shift column covers only the stored5/50/95 quantiles, not a uniform bound on all quantiles.

'''+common_table+'''
The g contrast subtracts the Beta(1.5,8.5) median from Beta(2,5) at fixed λ. The r contrast is σ=0.3 minus σ=1.2 at fixed magnitude family. Here g=P(χ_eff<0), while r=P(cos θ₁<0) retains mixed native references. Neither is an identified formation-channel fraction.

'''+contrast_table+'''
The saved detected contribution in all newly covered target support is zero for both policies. Consequently, full-union and old-BBH-domain-restricted mixture results agree. This numerical equality does not establish zero physical detection probability. Shared rows make these related estimators; their difference has no assigned significance.

## Added support and remaining omissions

Within the fixed target m₁≤100,m₂≥3 solar masses, BNS/NSBH proposals add no support. IMBH adds only m₁≥90,m₂≥10 at1.9<z≤2.3. This mass region comprises0.306136% of the target mass distribution, adding just0.102273 percentage points of total population probability. The proposal union still misses33.305257% of the target probability mass. These are source-population masses, not detection fractions or selection-error bars.

There are178 saved rows in the newly covered target region and no saved detections under either search policy. Weighted proposal sampling and hopeless selection prevent using this as178 independent direct-target Bernoulli trials.

A subsequent source diagnostic also finds that, in260,145 overlap rows supported only by the BBH proposal, the released mixture joint PDF is not exactly the separate BBH PDF times the generated-trial fraction. Its maximum relative residual is6.4486%, with median ratio0.9999995. The numerical construction history remains unverified. Each official joint PDF is retained as provided; the diagnostic does not recover the actual generation configuration.

## Verification, prior work and reproduction

The independent validator reconstructs Cartesian target/joint-proposal ratios, early/O3 integrals and all69 multiplicative event factors without producer weighting/posterior functions. Maximum relative selection difference is8.89e-15; relative-to-peak density difference1.94e-12; quantile difference4.17e-15. Piecewise-polynomial/Gauss–Legendre support quadratures differ by4.14e-10. This is numerical implementation verification with unchanged event-integral inputs, not independent observations.

Mixture injection use follows the official release and[Essick2021](https://doi.org/10.3847/2515-5172/ac2ba7), not a new method. See[Essick2023](https://arxiv.org/html/2307.02765v1) for hopeless cuts and the distinction between real searches and semianalytic approximations. The incremental result is this specific all69/12-model comparison, row identity, posterior distances and quantified remaining support.

Reproduction needs the fixed raw data and prior all69 caches, plus population-requirements.txt and Matplotlib. Raw HDF, row maps and numerical checkpoints are excluded from Git; official MD5 and local SHA256 are recorded. Fetching this one file does not restore all preceding research inputs.

'''+reproducible+links
    (BASE/'README.md').write_text(zh,encoding='utf-8',newline='\n')
    (BASE/'README.en.md').write_text(en,encoding='utf-8',newline='\n')
    print('Generated bilingual notes and identical data-derived tables.')


if __name__=='__main__':main()
