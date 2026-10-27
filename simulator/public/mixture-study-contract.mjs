// Generic finite-grid mixture display. It never fits data or interprets labels
// as physical membership, and never adds marginal medians to check closure.
const assert=(ok,message)=>{if(!ok)throw Error(message)};
const finite=x=>typeof x==='number'&&Number.isFinite(x);
const text=x=>typeof x==='string'&&x.length>0;
const probability=x=>finite(x)&&x>=0&&x<=1;
const triple=q=>Array.isArray(q)&&q.length===3&&q.every((x,i)=>probability(x)&&(!i||x>=q[i-1]));
export function validateMixtureStudy(d){
 assert(d?.schema==='blackhole-mixture-study/1','不支持的混合研究格式');
 assert(text(d.id)&&text(d.title)&&d.value_unit==='fraction'&&d.scientific_status==='conditional_fixed_bank_not_cosmic_rates','研究定义、单位或条件范围缺失');
 assert(Array.isArray(d.quantile_probabilities)&&d.quantile_probabilities.length===3&&d.quantile_probabilities.every(probability)&&d.quantile_probabilities[0]>0&&d.quantile_probabilities[0]<.5&&d.quantile_probabilities[1]===.5&&Math.abs(d.quantile_probabilities[0]+d.quantile_probabilities[2]-1)<1e-12,'必须声明中位数及等尾区间的概率');
 const unique=(a,max,label)=>{assert(Array.isArray(a)&&a.length>0&&a.length<=max&&a.every(x=>text(x.id)&&text(x.label)),label+'无效');assert(new Set(a.map(x=>x.id)).size===a.length,label+' ID 重复')};
 unique(d.components,20,'混合分量');unique(d.axes,6,'假设轴');unique(d.quantities,10,'科学量');
 assert(d.components.length>=2,'至少需要两个分量');
 assert(d.observations?.events>0&&Number.isInteger(d.observations.events)&&d.observations.posterior_rows>0&&Number.isInteger(d.observations.posterior_rows),'来源样本数无效');
 assert(text(d.domain)&&text(d.context)&&text(d.interpretation)&&Array.isArray(d.limits)&&d.limits.length>0&&d.limits.every(text),'定义和限制不完整');
 let combinations=1;for(const axis of d.axes){unique(axis.options,50,'选项');combinations*=axis.options.length;assert(axis.options.some(x=>x.id===d.defaults?.[axis.id]),'默认设定无效')}
 assert(combinations<=2000&&Array.isArray(d.cases)&&d.cases.length===combinations,'网格不完整或过大');
 const combinationsSeen=new Set();
 for(const r of d.cases){
  assert(text(r.id)&&r.settings,'设定标识缺失');
  for(const axis of d.axes)assert(axis.options.some(x=>x.id===r.settings[axis.id]),'未知选项');
  const k=d.axes.map(a=>r.settings[a.id]).join('|');assert(!combinationsSeen.has(k),'重复网格点');combinationsSeen.add(k);
  for(const q of d.quantities)assert(text(q.definition)&&Array.isArray(r.quantiles?.[q.id])&&r.quantiles[q.id].length===d.components.length&&r.quantiles[q.id].every(triple),'区间、分量或科学量不匹配');
  assert(Array.isArray(r.intrinsic_mean)&&r.intrinsic_mean.length===d.components.length&&r.intrinsic_mean.every(probability)&&Math.abs(r.intrinsic_mean.reduce((a,b)=>a+b,0)-1)<1e-8,'后验均值的组成不闭合');
  assert(triple(r.common_term)&&text(r.common_term_definition),'共同项语义或区间缺失');
 }
 assert(Array.isArray(d.papers)&&d.papers.every(p=>text(p.label)&&/^research-docs\/[a-zA-Z0-9._-]+\.pdf$/.test(p.url)),'论文路径无效');
 return d;
}
export function selectMixtureCase(d,settings){
 const r=d.cases.find(r=>d.axes.every(a=>r.settings[a.id]===settings[a.id]));assert(r,'没有匹配的已计算设定');return r;
}
