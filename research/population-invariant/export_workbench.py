"""Export validated population curves, semantics and diagnostics for the UI.

This is a display adapter. It does not execute model.py, refit event sources,
alter scientific outputs or infer an actual missing-contribution budget.
"""
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
from scipy.integrate import cumulative_trapezoid

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
CAT=ROOT/'research/population-catalog'
PUBLIC=ROOT/'simulator/public'


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def read(path):return json.loads(path.read_text(encoding='utf-8'))


def main():
    files={
        'native':CAT/'native-catalog-results.json','native_validation':CAT/'native-catalog-validation.json',
        'effective':BASE/'control-results.json','effective_validation':BASE/'control-validation.json',
        'o3':BASE/'paired-o3-results.json','o3_validation':BASE/'paired-o3-validation.json',
        'paired':BASE/'paired-pe-results.json','paired_validation':BASE/'paired-pe-validation.json',
        'support':BASE/'support-results.json','support_validation':BASE/'support-validation.json',
        'tail_plan':BASE/'tail-validation-plan.json'}
    d={k:read(p) for k,p in files.items()}
    for name in ['native','effective','o3','paired','support']:
        assert d[name+'_validation']['result_sha256']==sha(files[name])
    assert d['native']['completed_full_source_events']==d['effective']['event_count']==69
    assert d['paired_validation']['original_step_derivative_gate_passed'] is False
    with np.load(BASE/'control-checkpoints.npz',allow_pickle=False) as cp:
        effective_grid=cp['g'];effective_density=cp['densities']
    with np.load(BASE/'paired-o3-checkpoints.npz',allow_pickle=False) as cp:
        grid=cp['native_grid'];fs=cp['native_fs'];num=cp['native_numerator'];full=cp['full'];full_q=cp['full_quantiles']
    policies=[dict(id='matched',label='O3 匹配滤波并集',scope='不含 cWB；IFAR > 1 年。O1/O2 保持原固定处理。'),
        dict(id='union',label='O3 并集（含 cWB）',scope='加入 cWB；IFAR > 1 年。O1/O2 保持原固定处理。')]
    magnitudes=[dict(id='low',label='Beta(1.5,8.5)',parameters=[1.5,8.5]),dict(id='baseline',label='Beta(2,5)',parameters=[2,5])]
    native_curves={};native_quantiles={}
    for policy,p in [('matched',1),('union',0)]:
        eff=fs*full[p,0,:,0,None]+(1-fs)*full[p,0,:,1,None]
        logp=num-69*np.log(eff);density=np.exp(logp-logp.max(axis=1)[:,None]);density/=np.trapezoid(density,grid,axis=1)[:,None]
        cdf=cumulative_trapezoid(density,grid,axis=1,initial=0);cdf/=cdf[:,-1,None]
        quantiles=np.array([np.interp([.05,.5,.95],c,grid) for c in cdf])
        assert np.max(abs(quantiles-full_q[p,0]))<2e-11
        native_curves[policy]=density;native_quantiles[policy]=quantiles
    metrics=[]
    for metric_id,symbol,label,definition,unit,shapes,curve_grid in [
        ('effective_spin','g','负有效自旋 g','P(χ_eff < 0)','双星总体',[0,4,16],effective_grid),
        ('native_tilt','r','原生倾角 r · 混合参考诊断','P(cos θ₁ < 0)','主黑洞总体',[.3,.6,1.2],grid)]:
        models=[]
        for k in range(6):
            magnitude=magnitudes[k//3];shape=shapes[k%3]
            shape_label=f'λ = {shape:g}' if symbol=='g' else f'σ = {shape:g}'
            model_id=f'{metric_id}-{k}'
            results={}
            for policy,p in [('matched',0),('union',1)]:
                if symbol=='g':q=d['effective']['scenarios'][['matched_filter','with_cwb'][p]]['quantiles'][k];density=effective_density[p,k]
                else:q=native_quantiles[policy][k];density=native_curves[policy][k]
                c=cumulative_trapezoid(density,curve_grid,initial=0);c/=c[-1]
                assert np.max(abs(np.interp([.05,.5,.95],c,curve_grid)-q))<2e-11
                results[policy]=dict(quantiles=np.asarray(q).tolist(),density=density.tolist(),interval_kind='model_conditional_equal_tailed_90_percent')
            compare=k+3 if symbol=='g' and k<3 else k-3 if symbol=='g' else (k//3)*3+(2 if k%3==0 else 0)
            contrast_suffix=f'effective_amplitude:lambda{shape:g}' if symbol=='g' else 'native_width:'+('low_spin' if k<3 else 'baseline_spin')
            models.append(dict(id=model_id,label=magnitude['label']+' · '+shape_label,magnitude_id=magnitude['id'],shape_id=f'{shape:g}',
                comparison_model_id=f'{metric_id}-{compare}',paired_contrast_suffix=contrast_suffix,results=results))
        metrics.append(dict(id=metric_id,symbol=symbol,label=label,definition=definition,population_unit=unit,value_unit='fraction',
            prior=dict(family='Uniform',range=[float(curve_grid[0]),float(curve_grid[-1])]),grid=curve_grid.tolist(),
            shapes=[dict(id=f'{x:g}',label=('λ' if symbol=='g' else 'σ')+f' = {x:g}') for x in shapes],models=models,
            reference=('2PN χ_eff 守恒与各向同性角度测度保持为明确假设；不是精确全 GR，也未恢复实际注入参考。' if symbol=='g' else
                       '53 个事件为 20 Hz、15 个为 10 Hz、1 个为 16 Hz；这是混合原生参考的条件诊断。'),
            limitations=('固定的符号/偶形状族；可信区间重叠。不能把 g 当成形成渠道比例。' if symbol=='g' else
                         '共同 r 先验不修复参考频率差异；不能称为共同物理参考的倾角测量。')))
    contrasts=[]
    for row in d['paired']['contrasts']:
        policy,suffix=row['name'].split(':',1) if row['name'].startswith(('matched:','union:')) else (None,None)
        if policy is None:continue
        if suffix.startswith('effective_amplitude:'):
            metric_id='effective_spin';shape=int(suffix.split('lambda')[1]);k=[0,4,16].index(shape)
            positive=f'{metric_id}-{k+3}';negative=f'{metric_id}-{k}'
            label=f'g 幅度模型差 · λ = {shape}'
            definition='Beta(2,5) 中位数 − Beta(1.5,8.5) 中位数；同一 λ'
        elif suffix.startswith('native_width:'):
            metric_id='native_tilt';offset=0 if suffix.endswith('low_spin') else 3
            positive=f'{metric_id}-{offset}';negative=f'{metric_id}-{offset+2}'
            label='r 宽度模型差 · '+magnitudes[offset//3]['label'];definition='σ = 0.3 的中位数 − σ = 1.2 的中位数；同一幅度族'
        else:continue
        contrasts.append(dict(**row,id=row['name'],metric_id=metric_id,policy_id=policy,positive_model_id=positive,negative_model_id=negative,
            label=label,definition=definition,value_unit='percentage_point',range_kind='specified_artificial_thinning_percentiles',
            scope='固定已保存样本集合上的人工抽稀；不是天体物理可信区间、p 值或真实活动误差。'))
    support=d['support']
    completions=[]
    for row in support['event_completions']:
        results={}
        for policy,p in [('matched',0),('union',1)]:
            results[policy]={f'effective_spin-{k}':dict(quantiles=row['quantiles'][p][k],median_shift_pp=row['median_shift_pp'][p][k]) for k in range(6)}
        completions.append(dict(eta=row['eta'],eta_label=row['eta_label'],results=results))
    before={str(path):sha(path) for path in files.values()}
    papers=[];build_path=ROOT/'research/paper/v0.3/build-manifest.json';build=read(build_path)
    assert build['status']=='compiled'
    docdir=PUBLIC/'research-docs';docdir.mkdir(exist_ok=True)
    for item in build['editions']:
        src=ROOT/item['output'];assert sha(src)==item['pdf_sha256']
        target=docdir/src.name;shutil.copyfile(src,target)
        papers.append(dict(language=item['language'],url=target.relative_to(PUBLIC).as_posix(),sha256=sha(target),pages=item['pages'],version='0.3'))
    result=dict(schema='blackhole-population-study/1',id='gwtc3-full69-v03',title='全69来源：模型、有限样本与未见贡献',date='2026-10-25',
        scientific_status='conditional_not_channel_inference',observations=dict(events=69,posterior_rows=d['effective']['full_source_rows'],reference_counts_hz={'20':53,'10':15,'16':1}),
        policies=policies,magnitudes=magnitudes,metrics=metrics,
        defaults=dict(metric_id='effective_spin',model_id='effective_spin-3',policy_id='matched'),
        fixed_assumptions=['主质量：0.9×m₁⁻² + 0.1×截断高斯（35±5 M☉），范围 3–100 M☉。',
            '质量比密度 ∝ q；m₂ ≥ 3 M☉。红移密度 ∝ (dV_c/dz)(1+z)，z ≤ 2.3。',
            '未拟合质量–自旋或红移–自旋相关性；两个指标使用不同且明确的先验区间。',
            '每事件选一套 Phenom 族分析，不假设全目录配置完全相同。'],
        paired_banks=dict(replicates=512,rows_per_event=6913,o3_retained_rows=128454,o3_accepted_rows=284672,
            design='每事件及 O3 各自无放回等概率抽稀；同一索引在模型、指标和搜索规则间共享。固定 69 个观测与 O1/O2。',
            scope='这是人为指定的有限集合设计，不复现历史 cosmo 拒绝抽样，不校准 PE 链相关性或实际注入活动。',contrasts=contrasts),
        support=dict(kind='assumed_log_correction_envelope',estimated_from_data=False,
            formula='Q(expit(logit(q)−B)) ≤ Q_corrected(q) ≤ Q(expit(logit(q)+B))',
            scope='B 是未知 log 修正的振幅上界。此处只探索假设预算，不估计真实缺失贡献。保守排序保证失效不等于观测到物理翻转。',
            known_o3_domain=support['known_o3_omitted_domain'],completion_metric_id='effective_spin',completion_event=support['completion_event'],
            completions=completions,ordering_thresholds=support['ordering_certificate_thresholds'],
            tail_plan=dict(direct_draws=d['tail_plan']['total_direct_draws_if_all_zero'],executed=False,
                description='假设 0.5% 选择遗漏预算、理想独立直接抽样且零检出时的资源规划；历史检出协议须先核验。')),
        validation=dict(algorithmically_verified=True,physical_missing_budget_calibrated=False,real_campaign_coverage_calibrated=False,
            original_derivative_step_gate_passed=False,
            description='独立源坐标、协方差、索引指纹及完整后验乘积核查；原导数步长失败保留，后续诊断分开记录。'),
        provenance=dict(reports=[dict(id=k,repo_path=p.relative_to(ROOT).as_posix(),sha256=sha(p)) for k,p in files.items()],
            checkpoint_sha256={name:sha(BASE/name) for name in ['control-checkpoints.npz','paired-o3-checkpoints.npz']},
            public_sources=[dict(label='GWTC-2.1 PE v2',url='https://zenodo.org/records/6513631'),dict(label='GWTC-3 PE',url='https://zenodo.org/records/5546663'),dict(label='固定 O3 灵敏度发布',url='https://zenodo.org/records/7890437')],
            exporter_sha256=sha(Path(__file__)),paper_build_sha256=sha(build_path)),papers=papers)
    output=PUBLIC/'datasets/population-study.json'
    output.write_text(json.dumps(result,ensure_ascii=False,separators=(',',':'),allow_nan=False),encoding='utf-8')
    assert before=={str(path):sha(path) for path in files.values()}
    record=dict(schema='blackhole-population-display-export/1',output=output.relative_to(ROOT).as_posix(),sha256=sha(output),bytes=output.stat().st_size,
        models=sum(len(m['models']) for m in metrics),contrasts=len(contrasts),papers=papers,source_reports=result['provenance']['reports'],
        exporter_sha256=sha(Path(__file__)),scope='Display adapter only; source results preserved, no new physical likelihood/selection fit.')
    (BASE/'workbench-export.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:record[k] for k in ['output','sha256','bytes','models','contrasts']},indent=2))


if __name__=='__main__':main()
