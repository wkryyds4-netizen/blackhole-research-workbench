// Probability-functional display contract. No PE fit, resampling or calibrated
// missing-support error is inferred by these presentation calculations.
const finite=x=>typeof x==='number'&&Number.isFinite(x);
const object=x=>x!==null&&typeof x==='object'&&!Array.isArray(x);
const text=x=>typeof x==='string'&&x.length>0;
const hash=x=>typeof x==='string'&&/^[a-f0-9]{64}$/i.test(x);
const assert=(condition,message)=>{if(!condition)throw Error(message)};
const close=(a,b,tolerance=2e-6)=>Math.abs(a-b)<=tolerance;

export function curveCDF(grid,density){
 assert(Array.isArray(grid)&&Array.isArray(density)&&grid.length===density.length&&grid.length>=3,'曲线数组长度无效');
 let total=0;const cdf=[0];
 for(let i=0;i<grid.length;i++){
  assert(finite(grid[i])&&finite(density[i])&&density[i]>=0&&(i===0||grid[i]>grid[i-1]),'曲线必须有限、非负且网格严格递增');
  if(i){total+=(density[i]+density[i-1])*(grid[i]-grid[i-1])/2;cdf.push(total)}
 }
 assert(total>0&&finite(total),'曲线没有正的有限积分');
 return {cdf:cdf.map(v=>v/total),normalization:total};
}

export function quantileFromCDF(grid,cdf,q){
 assert(finite(q)&&q>=0&&q<=1,'分位概率必须在 [0,1] 内');
 if(q===0)return grid[0];if(q===1)return grid.at(-1);
 let lo=0,hi=cdf.length-1;
 while(lo<hi){const mid=(lo+hi)>>1;if(cdf[mid]<q)lo=mid+1;else hi=mid}
 if(lo===0)return grid[0];
 const a=lo-1,b=lo,weight=(q-cdf[a])/(cdf[b]-cdf[a]);
 return grid[a]+weight*(grid[b]-grid[a]);
}

export function medianEnvelope(grid,density,B){
 assert(finite(B)&&B>=0&&B<=20,'假设 B 需在 0–20 内');
 const {cdf}=curveCDF(grid,density),p=1/(1+Math.exp(B));
 return [quantileFromCDF(grid,cdf,p),quantileFromCDF(grid,cdf,1-p)];
}

export function validatePopulationStudy(d){
 assert(object(d)&&d.schema==='blackhole-population-study/1','不支持的群体研究格式');
 assert(text(d.id)&&text(d.title)&&d.scientific_status==='conditional_not_channel_inference','缺少研究标识或条件解释');
 assert(object(d.observations)&&Number.isInteger(d.observations.events)&&d.observations.events>0&&Number.isInteger(d.observations.posterior_rows)&&d.observations.posterior_rows>0,'观测与后验行数无效');
 const unique=(items,label)=>{assert(Array.isArray(items)&&items.length>0&&items.every(x=>object(x)&&text(x.id)&&text(x.label)),label+'缺失');assert(new Set(items.map(x=>x.id)).size===items.length,label+' ID 重复')};
 unique(d.policies,'选择规则');unique(d.magnitudes,'幅度族');unique(d.metrics,'科学指标');
 assert(d.metrics.length<=10&&d.policies.length<=8&&d.magnitudes.length<=20,'研究维度过大');
 const policies=new Set(d.policies.map(x=>x.id)),magnitudeIds=new Set(d.magnitudes.map(x=>x.id));
 const modelIds=new Set(),modelMap=new Map();
 for(const metric of d.metrics){
  assert(text(metric.symbol)&&text(metric.definition)&&text(metric.population_unit)&&text(metric.reference)&&text(metric.limitations),'指标定义或物理参考缺失');
  assert(metric.value_unit==='fraction','此版本仅显示概率积分量，不能将其他单位默认为百分比');
  assert(Array.isArray(metric.grid)&&metric.grid.length>=3&&metric.grid.length<=20001,'指标网格无效');
  assert(metric.grid.every((x,i)=>finite(x)&&x>=0&&x<=1&&(!i||x>metric.grid[i-1])),'指标网格需严格递增且位于 [0,1]');
  assert(metric.prior?.family==='Uniform'&&Array.isArray(metric.prior.range)&&metric.prior.range.length===2&&close(metric.prior.range[0],metric.grid[0],1e-12)&&close(metric.prior.range[1],metric.grid.at(-1),1e-12),'先验与显示支持不一致');
  unique(metric.shapes,'形状参数');unique(metric.models,'模型');assert(metric.models.length<=200,'模型数过大');
  const shapes=new Set(metric.shapes.map(x=>x.id));
  assert(new Set(metric.models.map(x=>x.magnitude_id+'|'+x.shape_id)).size===metric.models.length,'同一幅度与形状组合重复');
  for(const model of metric.models){
   assert(!modelIds.has(model.id),'跨指标模型 ID 重复');modelIds.add(model.id);modelMap.set(model.id,{model,metric});
   assert(magnitudeIds.has(model.magnitude_id)&&shapes.has(model.shape_id)&&text(model.paired_contrast_suffix),'模型参数或成对问题缺失');
   assert(object(model.results),'模型结果缺失');
   for(const policy of policies){
    const r=model.results[policy];assert(object(r)&&r.interval_kind==='model_conditional_equal_tailed_90_percent','区间类型必须明确为模型内可信区间');
    assert(Array.isArray(r.quantiles)&&r.quantiles.length===3&&r.quantiles.every((x,i)=>finite(x)&&x>=metric.grid[0]&&x<=metric.grid.at(-1)&&(!i||x>=r.quantiles[i-1])),'后验分位数无效');
    assert(r.density?.length===metric.grid.length,'密度与指标网格长度不符');
    const {cdf,normalization}=curveCDF(metric.grid,r.density);assert(close(normalization,1,2e-6),'后验密度未归一化');
    [.05,.5,.95].forEach((q,i)=>assert(close(quantileFromCDF(metric.grid,cdf,q),r.quantiles[i],2e-6),'曲线与已保存分位数不一致'));
   }
  }
 }
 for(const {model,metric} of modelMap.values())assert(metric.models.some(x=>x.id===model.comparison_model_id),'比较曲线必须属于同一指标');
 const banks=d.paired_banks;assert(object(banks)&&Number.isInteger(banks.replicates)&&banks.replicates>0&&text(banks.design)&&text(banks.scope),'抽稀设计缺失');
 assert(Array.isArray(banks.contrasts)&&banks.contrasts.length<=1000,'成对比较无效');
 const contrastIds=new Set();
 for(const c of banks.contrasts){
  assert(text(c.id)&&!contrastIds.has(c.id)&&text(c.label)&&text(c.definition)&&text(c.scope),'成对比较元数据缺失或重复');contrastIds.add(c.id);
  const a=modelMap.get(c.positive_model_id),b=modelMap.get(c.negative_model_id);
  assert(a&&b&&a.metric.id===b.metric.id&&a.metric.id===c.metric_id&&policies.has(c.policy_id),'禁止跨指标或未知模型/规则的相减');
  assert(c.value_unit==='percentage_point'&&c.range_kind==='specified_artificial_thinning_percentiles','抽稀范围不能标成可信区间');
  assert(finite(c.full_bank_value_pp)&&close(c.full_bank_value_pp,100*(a.model.results[c.policy_id].quantiles[1]-b.model.results[c.policy_id].quantiles[1]),3e-8),'成对差异与同一指标模型不匹配');
  for(const field of ['event_only_5_50_95_pp','selection_only_5_50_95_pp','joint_5_50_95_pp'])assert(Array.isArray(c[field])&&c[field].length===3&&c[field].every((x,i)=>finite(x)&&(!i||x>=c[field][i-1])),'人工抽稀分位范围无效');
  for(const field of ['event_nonpositive_fraction','joint_nonpositive_fraction'])assert(finite(c[field])&&c[field]>=0&&c[field]<=1,'抽稀条件频率无效');
 }
 for(const {model,metric} of modelMap.values())for(const policy of policies)assert(banks.contrasts.some(c=>c.id===policy+':'+model.paired_contrast_suffix&&c.metric_id===metric.id),'缺少模型关联的成对比较');
 assert(object(d.support)&&d.support.kind==='assumed_log_correction_envelope'&&d.support.estimated_from_data===false&&text(d.support.scope),'此界面只探索声明的假设预算，不接受伪装为已测误差的输入');
 if(d.support.completions!==undefined){
  const metric=d.metrics.find(x=>x.id===d.support.completion_metric_id);assert(metric&&Array.isArray(d.support.completions)&&d.support.completions.length>0&&d.support.completions.length<=100,'积分补入指标或情景缺失');
  assert(object(d.support.completion_event)&&text(d.support.completion_event.event)&&Number.isSafeInteger(d.support.completion_event.rows)&&d.support.completion_event.rows>0,'补入事件信息缺失');
  assert(new Set(d.support.completions.map(x=>x.eta_label)).size===d.support.completions.length,'积分补入标签重复');
  for(const completion of d.support.completions){
   assert(text(completion.eta_label)&&(completion.eta===null||(finite(completion.eta)&&completion.eta>=0))&&object(completion.results),'补入参数无效');
   for(const policy of policies)for(const model of metric.models){
    const r=completion.results[policy]?.[model.id];assert(object(r)&&Array.isArray(r.quantiles)&&r.quantiles.length===3&&r.quantiles.every((x,i)=>finite(x)&&x>=metric.grid[0]&&x<=metric.grid.at(-1)&&(!i||x>=r.quantiles[i-1])),'补入情景分位数无效');
    assert(finite(r.median_shift_pp)&&close(r.median_shift_pp,100*(r.quantiles[1]-model.results[policy].quantiles[1]),3e-8),'补入位移与基线不一致');
   }
  }
 }
 if(d.support.known_o3_domain!==undefined)assert(Array.isArray(d.support.known_o3_domain)&&d.support.known_o3_domain.every(x=>object(x)&&finite(x.population_redshift_tail_mass)&&x.population_redshift_tail_mass>=0&&x.population_redshift_tail_mass<=1&&finite(x.o3_time_weighted_selection_contribution_upper)&&x.o3_time_weighted_selection_contribution_upper>=0),'已知遗漏域摘要无效');
 if(d.support.tail_plan!==undefined)assert(object(d.support.tail_plan)&&Number.isSafeInteger(d.support.tail_plan.direct_draws)&&d.support.tail_plan.direct_draws>0&&d.support.tail_plan.executed===false&&text(d.support.tail_plan.description),'此版本的尾部规划必须明确为未执行');
 assert(object(d.provenance)&&Array.isArray(d.provenance.reports)&&d.provenance.reports.every(x=>text(x.repo_path)&&hash(x.sha256)),'研究来源与哈希缺失');
 assert(Array.isArray(d.fixed_assumptions)&&d.fixed_assumptions.every(text),'固定假设缺失');
 if(d.papers!==undefined)assert(Array.isArray(d.papers)&&d.papers.every(x=>text(x.language)&&/^research-docs\/[a-zA-Z0-9._-]+\.pdf$/.test(x.url)&&hash(x.sha256)&&Number.isInteger(x.pages)),'论文链接或哈希无效');
 validateStudySettings(d,d.defaults);
 return d;
}

export function validateStudySettings(data,settings){
 assert(object(settings),'研究视图配置缺失');
 const metric=data.metrics.find(x=>x.id===settings.metric_id);
 assert(metric&&metric.models.some(x=>x.id===settings.model_id)&&data.policies.some(x=>x.id===settings.policy_id),'未知指标、模型或选择规则');
 if(settings.assumed_B!==null&&settings.assumed_B!==undefined)assert(finite(settings.assumed_B)&&settings.assumed_B>=0&&settings.assumed_B<=1.5,'此视图的假设 B 范围为 0–1.5');
 if(settings.completion_eta!==null&&settings.completion_eta!==undefined)assert(metric.id===data.support.completion_metric_id&&data.support.completions?.some(x=>x.eta_label===settings.completion_eta),'未提供该指标的积分补入情景');
 return {...settings,assumed_B:settings.assumed_B??null,completion_eta:settings.completion_eta??null};
}

export function selectedStudyResult(data,rawSettings){
 const settings=validateStudySettings(data,rawSettings),metric=data.metrics.find(x=>x.id===settings.metric_id),model=metric.models.find(x=>x.id===settings.model_id);
 const comparison=metric.models.find(x=>x.id===model.comparison_model_id),posterior=model.results[settings.policy_id];
 const contrast=data.paired_banks.contrasts.find(x=>x.id===settings.policy_id+':'+model.paired_contrast_suffix);
 let assumed=null;
 if(settings.assumed_B!==null){
  const a=metric.models.find(x=>x.id===contrast.positive_model_id),b=metric.models.find(x=>x.id===contrast.negative_model_id);
  const positive=medianEnvelope(metric.grid,a.results[settings.policy_id].density,settings.assumed_B),negative=medianEnvelope(metric.grid,b.results[settings.policy_id].density,settings.assumed_B);
  assumed={B:settings.assumed_B,estimated_from_data:false,median_bounds:medianEnvelope(metric.grid,posterior.density,settings.assumed_B),
   contrast_outer_bounds_pp:[100*(positive[0]-negative[1]),100*(positive[1]-negative[0])],scope:data.support.scope};
 }
 return {settings,metric,model,comparison,posterior,contrast,assumed};
}

export function populationStudySession(data,settings,sourceHash){
 assert(hash(sourceHash),'缺少导出数据的 SHA-256');
 const selected=selectedStudyResult(data,settings);
 return {schema:'blackhole-population-study-session/1',dataset_sha256:sourceHash,dataset_id:data.id,
  settings:selected.settings,result:{metric_id:selected.metric.id,definition:selected.metric.definition,population_unit:selected.metric.population_unit,
   model_id:selected.model.id,policy_id:selected.settings.policy_id,conditional_quantiles:selected.posterior.quantiles,
   interval_kind:selected.posterior.interval_kind,paired_contrast:selected.contrast,assumed_correction:selected.assumed},
  interpretation:'条件显示与声明的假设探索；未推断形成渠道、实际遗漏预算或活动覆盖率。',study:data};
}
