"""Data-derived bilingual milestone reports; no new scientific calculation."""
from pathlib import Path
import json

BASE=Path(__file__).resolve().parent


def main():
    s=json.loads((BASE/'branching-synthesis.json').read_text());fits=json.loads((BASE/'branching-results.json').read_text())
    assert s['cases']==24
    maximum=s['contrast_maxima']
    for lang in ['zh','en']:
        zh=lang=='zh'
        title='# 五渠道条件推断：24组完整对照' if zh else '# Five-channel conditional inference: 24 complete comparisons'
        intro=('2026-10-27。全部69事件、11,180,578后验行，完整事件与选择积分、两套后验算法及共同正则项拆分均已完成。'
               '结论限定于声明的限制域、固定出生自旋0.1/CE效率2及额外自旋补全；不是宇宙形成比例测量。') if zh else (
               '2026-10-27. All69 events and11,180,578 posterior rows are retained. Event/selection integrals, two posterior algorithms and the common-term decomposition are complete. '
               'Results are conditional on the restricted domain, fixed birth spin0.1/CE efficiency2 and declared spin completions; they are not cosmic formation-fraction measurements.')
        rows=[]
        names=['正则化概率','隐藏自旋补全','O3银行替换','搜索口径'] if zh else ['Regularizer probability','Hidden-spin completion','O3 bank replacement','Search policy']
        for name,key in zip(names,['regularization','hidden_spin','bank','search_policy']):
            r=maximum[key];rows.append(f'| {name} | {r["value_pp"]:.3f} | {r["channel"]} |')
        table=('| 比较轴 | 最大中位数移动 / 百分点 | 渠道 |\n|---|---:|---|\n' if zh else '| Comparison | Largest median movement / pp | Channel |\n|---|---:|---|\n')+'\n'.join(rows)
        description=(
            '最大值分别在各自网格点取得，不能相加或当作可信区间。两套算法的分位数最大差为'
            f'{s["maximum_independent_mcmc_quantile_difference_pp"]:.3f}个百分点；四次重要性重复最大跨度为'
            f'{s["maximum_qmc_replicate_quantile_span_pp"]:.3f}个百分点。归一化抽样半部替换的最大移动为'
            f'{s["maximum_normalizer_half_quantile_shift_pp"]:.3f}个百分点。') if zh else (
            'Maxima occur at separate grid points and are neither additive errors nor credible intervals. The largest independent-algorithm quantile difference is'
            f' {s["maximum_independent_mcmc_quantile_difference_pp"]:.3f}pp; the maximum four-replicate importance span is'
            f' {s["maximum_qmc_replicate_quantile_span_pp"]:.3f}pp. Replacing normalizers with independent draw halves shifts quantiles by at most'
            f' {s["maximum_normalizer_half_quantile_shift_pp"]:.3f}pp.')
        interpretation=(
            '官方混合银行、匹配滤波、均匀基准条件下，CE内禀权重中位数从83.85%（无正则化）到76.17%（1%共同项）。'
            '但共同项的计算效率约8.06%，远高于未正则化CE约0.139%；它在检出模型混合中的责任份额中位数达到12.14%。'
            '这说明“小密度正则项”不必是“小选择修正”，不表示测得12.14%的未知物理渠道。'
            '\n\n`branching-results.json`的`detected_quantiles`是**正则化渠道标签**。共同均匀项按beta分配给标签；'
            '`detected-component-results.json`另给出模型项、共同项和去共同项后的条件检出比例。边缘中位数不必加和为1。') if zh else (
            'For the official-mixture matched-filter case with the uniform base completion, the CE intrinsic median is83.85% without regularization and76.17% with a1% common term. '
            'The computed efficiency of the common term is about8.06%, versus0.139% for unregularized CE, so its median responsibility in the detected model mixture reaches12.14%. '
            'A small density regularizer therefore need not be a small selection correction. This does not measure a12.14% unknown physical channel.'
            '\n\nThe original `detected_quantiles` are **regularized channel labels**; shared U is allocated to labels according to beta. '
            '`detected-component-results.json` separately reports model terms, U, and the conditional detected model terms. Marginal medians need not sum to1.')
        limits=(
            f'最小事件贡献集中量只有{s["minimum_event_concentration"]:.2f}，最小选择集中量{s["minimum_selection_concentration"]:.2f}。'
            '这些值不是PE链或活动ESS，也未用于删除事件。后验求值稳定不代表银行尾部已经校准。历史注入参考/筛选、'
            'PE未见支持、波形误差和真实渠道分量自旋仍未解决。两个条件补全只是压力模型，并非已验证的形成物理。'
            '\n\n下一项实质工作是对本批渠道积分做成对有限银行传播，不能复用旧g/r模型的抽稀范围。随后再检验渠道特定自旋补全和超参数变化。') if zh else (
            f'The minimum event contribution concentration is{s["minimum_event_concentration"]:.2f}, and selection concentration{s["minimum_selection_concentration"]:.2f}. '
            'These are not PE-chain/campaign ESS values and do not trigger deletion. Stable posterior computation does not calibrate bank tails. '
            'Historical injection reference/screening, missing PE support, waveform error and true channel-specific component spins remain unresolved. '
            'Both completions are explicit stress models, not verified formation physics.'
            '\n\nThe next substantive step is paired finite-bank propagation for these channel integrals, without reusing old g/r thinning intervals, then channel-specific spin completions and hyperparameter variations.')
        links=('## 图表、论文与证据' if zh else '## Figures, papers and evidence')+'\n\n'
        links+='- [Posterior figure](channel-posteriors.svg)\n- [Sensitivity and regularizer](channel-sensitivity.svg)\n'
        links+='- [中文 LaTeX PDF](../../output/pdf/blackhole-channel-completion-latex-zh-2026-10-27.pdf) / [English LaTeX PDF](../../output/pdf/blackhole-channel-completion-latex-en-2026-10-27.pdf)\n'
        links+='- [LaTeX projects and shared facts](../paper/channels-2026-10-27/)\n- [Method and literature comparison](measure-method.md)\n- [Complete synthesis](branching-synthesis.json)\n- [Independent posterior validation](observed-mcmc-validation.json)\n'
        commands='''
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
'''
        repro=('## 复现\n\n先恢复来源模型及此前已验证的69事件和注入缓存；环境与来源见本目录README。'
               '在仓库根目录运行以下命令。缓存协议与源码哈希绑定，重跑校准/改动设计后不要混用旧输出；先保留当前版本。'
               '不运行根目录model.py。原直接乘积失败与短链失败均保留。\n') if zh else (
               '## Reproduction\n\nRestore model products and the previously verified69-event/injection caches first; see the directory README. '
               'Run from the repository root. Cache protocols include source hashes; preserve this snapshot before recalibration or design changes, and never mix protocols. '
               'Do not execute root model.py. Original direct-product and short-chain failures remain preserved.\n')
        text='\n\n'.join([title,intro,table,description,interpretation,'## '+('限制与下一步' if zh else 'Limits and next step'),limits,links,repro])+commands
        (BASE/('RESULTS.md' if zh else 'RESULTS.en.md')).write_text(text,encoding='utf-8')


if __name__=='__main__':main()
