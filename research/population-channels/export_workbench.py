"""Export a generic finite-grid mixture study with explicit quantity semantics."""
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np

BASE=Path(__file__).resolve().parent;ROOT=BASE.parents[1];PUBLIC=ROOT/'simulator/public'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    s=json.loads((BASE/'branching-synthesis.json').read_text());fits=json.loads((BASE/'branching-results.json').read_text())
    components=json.loads((BASE/'detected-component-results.json').read_text());mc=json.loads((BASE/'observed-mcmc-validation.json').read_text());assert mc['status']=='pass'
    design=json.loads((BASE/'likelihood-design.json').read_text())
    axes=[dict(id='family',label='隐藏自旋补全',options=[dict(id='uniform_magnitude',label='条件基准 U'),dict(id='beta2_5_magnitude',label='条件基准 B')]),
          dict(id='bank',label='O3 注入银行',options=[dict(id='official_mixture',label='官方四总体混合'),dict(id='pure_bbh',label='纯 BBH')]),
          dict(id='policy',label='搜索口径',options=[dict(id='matched_filter',label='匹配滤波'),dict(id='with_cwb',label='并入 cWB')]),
          dict(id='epsilon',label='共同正则项',options=[dict(id='0',label='0%'),dict(id='0.001',label='0.1%'),dict(id='0.01',label='1%')])]
    quantities=[dict(id='intrinsic',label='域内内禀分支权重 β',definition='固定模型与限制域内的分支权重；不是全部宇宙并合的形成比例。'),
                dict(id='detected_label',label='检出标签权重 η(label)',definition='包含共同正则项的渠道标签；共同项按 β 分配给标签，不是独立事件分类。'),
                dict(id='detected_model',label='检出模型项 η(model)',definition='去掉共同正则项的五个检出模型项；与共同项一起加和为 1。')]
    data=dict(schema='blackhole-mixture-study/1',id='channel-completion-2026-10-27',title='形成渠道：条件补全与选择敏感性',
              scientific_status='conditional_fixed_bank_not_cosmic_rates',value_unit='fraction',quantile_probabilities=[.05,.5,.95],axes=axes,quantities=quantities,
              components=[dict(id=c,label=label) for c,label in zip(s['channels'],['CE · 共同包层','CHE · 化学均匀','GC · 球状星团','NSC · 核星团','SMT · 稳定转移'])],
              observations=dict(events=69,posterior_rows=11180578),context='出生自旋 0.1；CE 效率 2；Dirichlet(1,1,1,1,1) 域内先验。',
              domain='源系 m1 为 5–100、m2≥3 太阳质量；q≤1、0<z≤1.9、|有效自旋|<0.99；分量幅度<0.99。',
              interpretation='这是已计算网格的读取与比较。U/B 是构造隐藏自旋条件的基准律，提升后的幅度边缘一般不再是均匀/Beta。区间为模型内 90% 等尾区间；边缘中位数不必加和为 1。',
              limits=['尚未恢复真实渠道分量自旋、完整注入执行参考、PE未见支持及真实活动覆盖率。',
                      f'最小事件贡献集中量 {s["minimum_event_concentration"]:.2f}；最小选择集中量 {s["minimum_selection_concentration"]:.2f}。这些不是链/活动ESS，也未触发事件删除。',
                      f'两套后验算法分位数最大差 {s["maximum_independent_mcmc_quantile_difference_pp"]:.3f} 个百分点；这不校准保存银行本身。'],
              defaults=dict(family='uniform_magnitude',bank='official_mixture',policy='matched_filter',epsilon='0'),cases=[],papers=[],
              sources=[dict(title='公开渠道模型及许可',url='https://zenodo.org/records/14967688'),dict(title='Colloms 等模型与方法',url='https://arxiv.org/abs/2503.03819v1')])
    for key,r in fits['cases'].items():
        derived=components['cases'][key]
        data['cases'].append(dict(id=key,settings=dict(family=r['family'],bank=r['bank'],policy=r['policy'],epsilon=f'{r["epsilon"]:g}'),
             quantiles=dict(intrinsic=r['intrinsic_quantiles'],detected_label=r['detected_quantiles'],detected_model=derived['model_term_quantiles']),
             intrinsic_mean=np.mean([rep['intrinsic_mean'] for rep in r['replicates']],axis=0).tolist(),
             common_term=derived['uniform_term_quantiles'],common_term_definition='共同正则项在检出混合中的数学责任份额；不是测得的异常源或新增物理渠道。'))
    for lang,label in [('zh','中文研究稿 PDF'),('en','English research note PDF')]:
        name=f'blackhole-channel-completion-latex-{lang}-2026-10-27.pdf';source=ROOT/'output/pdf'/name
        destination=PUBLIC/'research-docs'/name;shutil.copyfile(source,destination)
        data['papers'].append(dict(language=lang,label=label,url='research-docs/'+name,sha256=sha(destination)))
    path=PUBLIC/'datasets/channel-study.json';path.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    sources=[BASE/'branching-synthesis.json',BASE/'branching-results.json',BASE/'detected-component-results.json',BASE/'observed-mcmc-validation.json',BASE/'likelihood-design.json']
    manifest=dict(schema='blackhole-mixture-workbench-export/1',sha256=sha(path),exporter_sha256=sha(Path(__file__)),source_reports=[dict(repo_path=p.relative_to(ROOT).as_posix(),sha256=sha(p)) for p in sources],papers=data['papers'])
    (BASE/'workbench-export.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(cases=len(data['cases']),data_bytes=path.stat().st_size,sha256=manifest['sha256'])))


if __name__=='__main__':main()
