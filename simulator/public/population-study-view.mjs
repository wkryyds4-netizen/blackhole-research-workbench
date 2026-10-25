import {validatePopulationStudy,validateStudySettings,selectedStudyResult,populationStudySession,curveCDF,quantileFromCDF} from './population-study-contract.mjs';

const COLORS={selected:'#62d9db',comparison:'#eeb977',prior:'#b7a2ef',event:'#78bb91',o3:'#62b2da',joint:'#e69a69'};
const pct=x=>(100*x).toFixed(2)+'%';
const pp=x=>(Math.abs(x)<.000005?0:x).toFixed(3);
const digest=async bytes=>Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',bytes)),x=>x.toString(16).padStart(2,'0')).join('');

export function createPopulationStudyView({url='datasets/population-study.json'}={}){
 const $=id=>document.getElementById(id),root=$('population-latest');
 if(!root||!$('population-study-status'))throw Error('缺少全来源群体研究容器');
 // This is a fixed application template. Imported labels/results are always
 // assigned with textContent, never interpolated into HTML.
 root.innerHTML=`
 <div class="population-controls study-controls">
  <label>科学指标<select id="study-metric" disabled></select></label>
  <label>自旋幅度族<select id="study-magnitude" disabled></select></label>
  <label>角度形状<select id="study-shape" disabled></select></label>
  <label>O3 搜索口径<select id="study-policy" disabled></select></label>
 </div>
 <div class="study-actions"><button id="study-reset" disabled>重置研究视图</button><button id="study-export" disabled>导出当前研究快照</button><button id="study-import">导入研究快照</button><input id="study-file" type="file" accept=".json,application/json" hidden><span id="study-origin" class="tiny muted"></span></div>
 <p id="study-definition" class="study-definition"></p><p id="study-reference" class="notice"></p>
 <div class="study-interval-types"><span>模型内：90% 可信区间</span><span>抽稀：人工重复的 5%–95% 范围</span><span>未知修正：假设性外包络</span></div>
 <div class="population-study-grid">
  <article class="research-panel"><h2>完整来源的条件后验</h2><div id="study-legend" class="legend"></div><canvas id="study-posterior" aria-label="同一科学指标的模型内条件后验曲线"></canvas><p id="study-interval" class="notice"></p><p id="study-comparison" class="tiny"></p><details class="study-detail"><summary>六个条件模型的数值表</summary><div class="table-scroll"><table id="study-model-table"></table></div></details></article>
  <article class="research-panel"><h2 id="study-contrast-title">成对有限样本检查</h2><p id="study-contrast-definition"></p><canvas id="study-paired" aria-label="人工抽稀的成对中位数差，单位为百分点"></canvas><p id="study-paired-summary" class="notice"></p><div class="table-scroll"><table id="study-paired-table"></table></div><p id="study-paired-scope"></p></article>
  <article class="research-panel"><h2>未知修正的条件边界</h2><p>横轴 B 是未知 log 修正的振幅上界，当前数据没有测出它。阴影展示假设边界，不能当成可信区间。</p><label class="study-check"><input id="study-budget-enable" type="checkbox">探索一个假设预算</label><div class="study-budget"><label for="study-budget">假设 B <output id="study-budget-value">未启用</output></label><input id="study-budget" type="range" min="0" max="1.5" step=".05" value=".5" disabled></div><canvas id="study-envelope" aria-label="成对中位数差随假设 B 的外包络，B未由数据估计"></canvas><p id="study-budget-note" class="notice"></p><p id="study-support-scope"></p></article>
  <article class="research-panel"><h2>已见样本之外还缺什么</h2><p id="study-tail-note"></p><p id="study-plan-note"></p><details class="study-detail" id="study-completion"><summary>单事件未见分量：声明的补入压力族</summary><p id="study-completion-event"></p><label>补入情景<select id="study-completion-eta"></select></label><p id="study-completion-result" class="notice"></p><p>η 是负/正分量积分比，不是缺失后验概率。其余事件与选择保持固定；不声称这些情景来自同一物理似然。此情景不会改写上方基线后验。</p></details><p id="study-limitations"></p></article>
 </div>
 <details class="source-detail study-detail"><summary>固定假设、来源与验证范围</summary><ul id="study-assumptions"></ul><p id="study-validation"></p><div id="study-source-links" class="study-links"></div><div id="study-papers" class="study-links"></div><pre id="study-source"></pre></details>`;
 let data=null,contentHash=null,hashScope='fetched_json_bytes',imported=false,declaredOriginHash=null,loadToken=0;

 function fillSelect(id,items,value){const e=$(id);e.replaceChildren();for(const item of items){const option=document.createElement('option');option.value=item.id;option.textContent=item.label;e.append(option)}e.value=value??items[0]?.id;if(e.selectedIndex<0)e.selectedIndex=0;e.disabled=false}
 function metric(){return data.metrics.find(x=>x.id===$('study-metric').value)}
 function model(){const m=metric();return m.models.find(x=>x.magnitude_id===$('study-magnitude').value&&x.shape_id===$('study-shape').value)||m.models[0]}
 function settings(){return {metric_id:metric().id,model_id:model().id,policy_id:$('study-policy').value,assumed_B:$('study-budget-enable').checked?+$('study-budget').value:null,
  completion_eta:metric().id===data.support.completion_metric_id?($('study-completion-eta').value||null):null}}
 function populateModels(selected){const m=metric(),pick=m.models.find(x=>x.id===selected)||m.models[0];fillSelect('study-magnitude',data.magnitudes.filter(x=>m.models.some(r=>r.magnitude_id===x.id)),pick.magnitude_id);fillSelect('study-shape',m.shapes.filter(x=>m.models.some(r=>r.magnitude_id===pick.magnitude_id&&r.shape_id===x.id)),pick.shape_id)}
 function populateCompletion(selected=null){const enabled=metric().id===data.support.completion_metric_id&&data.support.completions?.length;$('study-completion').hidden=!enabled;
  if(enabled)fillSelect('study-completion-eta',[{id:'',label:'未选择补入情景'},...data.support.completions.map(x=>({id:x.eta_label,label:x.eta===null?'η → ∞（归一化极限）':`η = ${x.eta}`}))],selected??'')}
 function canvas(id){const e=$(id),r=e.getBoundingClientRect();if(r.width<1||r.height<1)return null;const k=Math.min(globalThis.devicePixelRatio||1,2);e.width=Math.round(r.width*k);e.height=Math.round(r.height*k);const c=e.getContext('2d');c.setTransform(k,0,0,k,0,0);c.clearRect(0,0,r.width,r.height);return {c,w:r.width,h:r.height}}
 function line(c,points,color,dashed=false,width=1.8){c.strokeStyle=color;c.lineWidth=width;c.setLineDash(dashed?[5,4]:[]);c.beginPath();points.forEach(([x,y],i)=>i?c.lineTo(x,y):c.moveTo(x,y));c.stroke();c.setLineDash([])}
 function posteriorChart(selected){const v=canvas('study-posterior');if(!v)return;const {c,w,h}=v,m=selected.metric,L=47,R=17,T=26,B=40,pw=w-L-R,ph=h-T-B,lo=m.grid[0],hi=m.grid.at(-1);
  const comparison=selected.comparison.results[selected.settings.policy_id].density,prior=1/(hi-lo),max=Math.max(prior,...selected.posterior.density,...comparison)*1.12;
  const X=x=>L+(x-lo)/(hi-lo)*pw,Y=y=>T+ph*(1-y/max);c.font='10px system-ui';c.fillStyle='#94a7b9';
  for(let j=0;j<=4;j++){const y=T+ph*j/4;c.strokeStyle='#253345';c.beginPath();c.moveTo(L,y);c.lineTo(w-R,y);c.stroke();c.textAlign='right';c.fillText((max*(1-j/4)).toFixed(1),L-6,y+4);c.textAlign='center';c.fillText((100*(lo+(hi-lo)*j/4)).toFixed(hi-lo<.5?1:0)+'%',L+pw*j/4,h-21)}
  c.textAlign='left';c.fillText('概率密度（每单位比例）',L,13);c.textAlign='center';c.fillText(m.symbol+' · '+m.population_unit,L+pw/2,h-4);
  c.save();c.beginPath();c.rect(L,T,pw,ph);c.clip();line(c,m.grid.map(x=>[X(x),Y(prior)]),COLORS.prior,true);line(c,m.grid.map((x,i)=>[X(x),Y(comparison[i])]),COLORS.comparison,true);line(c,m.grid.map((x,i)=>[X(x),Y(selected.posterior.density[i])]),COLORS.selected);c.restore();
 }
 function pairedChart(contrast){const v=canvas('study-paired');if(!v)return;const {c,w,h}=v,L=78,R=18,T=24,B=38,pw=w-L-R,ph=h-T-B;
  const rows=[['完整样本',[contrast.full_bank_value_pp,contrast.full_bank_value_pp,contrast.full_bank_value_pp],'#ffffff'],['仅事件',contrast.event_only_5_50_95_pp,COLORS.event],['仅 O3',contrast.selection_only_5_50_95_pp,COLORS.o3],['联合',contrast.joint_5_50_95_pp,COLORS.joint]];
  let lo=Math.min(0,...rows.map(r=>r[1][0])),hi=Math.max(0,...rows.map(r=>r[1][2])),pad=Math.max((hi-lo)*.12,.06);lo-=pad;hi+=pad;const X=x=>L+(x-lo)/(hi-lo)*pw;c.font='10px system-ui';
  for(let j=0;j<=4;j++){const x=lo+(hi-lo)*j/4;c.strokeStyle='#253345';c.beginPath();c.moveTo(X(x),T);c.lineTo(X(x),h-B);c.stroke();c.fillStyle='#94a7b9';c.textAlign='center';c.fillText(x.toFixed(2),X(x),h-21)}
  c.strokeStyle='#778b9d';c.beginPath();c.moveTo(X(0),T);c.lineTo(X(0),h-B);c.stroke();
  rows.forEach(([label,q,color],i)=>{const y=T+(i+.5)*ph/4;c.fillStyle='#b7c5d3';c.textAlign='right';c.fillText(label,L-10,y+4);line(c,[[X(q[0]),y],[X(q[2]),y]],color);c.strokeStyle=color;for(const x of [q[0],q[2]]){c.beginPath();c.moveTo(X(x),y-4);c.lineTo(X(x),y+4);c.stroke()}c.fillStyle=color;c.beginPath();if(i===0){c.moveTo(X(q[1]),y-5);c.lineTo(X(q[1])+5,y);c.lineTo(X(q[1]),y+5);c.lineTo(X(q[1])-5,y);c.closePath()}else c.arc(X(q[1]),y,3.5,0,2*Math.PI);c.fill()});
  c.fillStyle='#94a7b9';c.textAlign='left';c.fillText('人工抽稀 · 5/50/95 分位数',L,12);c.textAlign='center';c.fillText('中位数差 / pp',L+pw/2,h-4);
 }
 function envelopeChart(selected){const v=canvas('study-envelope');if(!v)return;const {c,w,h}=v,m=selected.metric,contrast=selected.contrast,L=47,R=18,T=27,B=40,pw=w-L-R,ph=h-T-B;
  const positive=m.models.find(x=>x.id===contrast.positive_model_id),negative=m.models.find(x=>x.id===contrast.negative_model_id),ca=curveCDF(m.grid,positive.results[selected.settings.policy_id].density).cdf,cb=curveCDF(m.grid,negative.results[selected.settings.policy_id].density).cdf;
  const rows=Array.from({length:61},(_,i)=>{const b=i*1.5/60,p=1/(1+Math.exp(b));return {b,lo:100*(quantileFromCDF(m.grid,ca,p)-quantileFromCDF(m.grid,cb,1-p)),hi:100*(quantileFromCDF(m.grid,ca,1-p)-quantileFromCDF(m.grid,cb,p))}});
  let lo=Math.min(0,...rows.map(x=>x.lo)),hi=Math.max(0,...rows.map(x=>x.hi)),pad=Math.max((hi-lo)*.08,.05);lo-=pad;hi+=pad;const X=x=>L+x/1.5*pw,Y=y=>T+(hi-y)/(hi-lo)*ph;c.font='10px system-ui';
  for(let j=0;j<=4;j++){const value=lo+(hi-lo)*j/4,y=Y(value);c.strokeStyle='#253345';c.beginPath();c.moveTo(L,y);c.lineTo(w-R,y);c.stroke();c.fillStyle='#94a7b9';c.textAlign='right';c.fillText(value.toFixed(1),L-6,y+4);c.textAlign='center';c.fillText((1.5*j/4).toFixed(2),L+pw*j/4,h-21)}
  c.fillStyle='#b7a2ef25';c.beginPath();rows.forEach((r,i)=>i?c.lineTo(X(r.b),Y(r.hi)):c.moveTo(X(r.b),Y(r.hi)));for(const r of [...rows].reverse())c.lineTo(X(r.b),Y(r.lo));c.closePath();c.fill();
  line(c,rows.map(r=>[X(r.b),Y(r.lo)]),COLORS.prior);line(c,rows.map(r=>[X(r.b),Y(r.hi)]),COLORS.prior,true);line(c,[[L,Y(0)],[w-R,Y(0)]],'#778b9d',true,1);
  if(selected.assumed)line(c,[[X(selected.assumed.B),T],[X(selected.assumed.B),h-B]],COLORS.joint,true,1);
  c.fillStyle='#94a7b9';c.textAlign='left';c.fillText('成对中位数差外包络 / pp',L,13);c.textAlign='center';c.fillText('假设 B（未由数据估计）',L+pw/2,h-4);
 }
 function table(id,headers,rows){const e=$(id);e.replaceChildren();const tr=e.createTHead().insertRow();for(const title of headers){const th=document.createElement('th');th.scope='col';th.textContent=title;tr.append(th)}const body=e.createTBody();for(const row of rows){const r=body.insertRow();for(const value of row){const td=r.insertCell();if(value instanceof Node)td.append(value);else td.textContent=String(value)}}}
 function legend(selected){const e=$('study-legend');e.replaceChildren();for(const [label,color] of [[selected.model.label,COLORS.selected],[selected.comparison.label,COLORS.comparison],['固定均匀先验',COLORS.prior]]){const span=document.createElement('span'),swatch=document.createElement('i');swatch.className='swatch';swatch.style.background=color;span.append(swatch,document.createTextNode(label));e.append(span)}}
 function render(){if(!data||$('population-view').hidden)return;const s=selectedStudyResult(data,settings()),m=s.metric,q=s.posterior.quantiles,other=s.comparison.results[s.settings.policy_id].quantiles,c=s.contrast;
  $('study-definition').textContent=`${m.symbol} = ${m.definition} · ${m.population_unit}。先验：Uniform[${m.prior.range.map(x=>x.toFixed(5)).join(', ')}]。`;
  $('study-reference').textContent=m.reference;
  $('study-interval').textContent=`${s.model.label}：中位数 ${pct(q[1])}；模型内 90% 等尾可信区间 ${pct(q[0])}–${pct(q[2])}。`;
  $('study-comparison').textContent=`同指标对照 ${s.comparison.label}：中位数 ${pct(other[1])}，90% 区间 ${pct(other[0])}–${pct(other[2])}。这些模型不等权合并。`;
  legend(s);posteriorChart(s);pairedChart(c);envelopeChart(s);
  $('study-contrast-title').textContent=c.label;$('study-contrast-definition').textContent=c.definition;
  $('study-paired-summary').textContent=`完整样本的差为 ${pp(c.full_bank_value_pp)} pp。${data.paired_banks.replicates} 次联合人工抽稀中，差非正的频率为 ${(100*c.joint_nonpositive_fraction).toFixed(2)}%；这不是 p 值。`;
  table('study-paired-table',['抽稀设计','5%','50%','95%'],[['仅事件',...c.event_only_5_50_95_pp.map(pp)],['仅 O3',...c.selection_only_5_50_95_pp.map(pp)],['联合',...c.joint_5_50_95_pp.map(pp)]]);
  $('study-paired-scope').textContent=`每事件保留 ${data.paired_banks.rows_per_event.toLocaleString()} 行；O3 保留 ${data.paired_banks.o3_retained_rows.toLocaleString()}/${data.paired_banks.o3_accepted_rows.toLocaleString()} 行。${c.scope} ${data.paired_banks.scope}`;
  $('study-budget').disabled=!$('study-budget-enable').checked;$('study-budget-value').textContent=s.assumed?s.assumed.B.toFixed(2):'未启用';
  $('study-budget-note').textContent=s.assumed?`假设 B = ${s.assumed.B.toFixed(2)}：所选 ${m.symbol} 中位数的条件边界为 ${s.assumed.median_bounds.map(pct).join('–')}；成对差的外包络为 [${s.assumed.contrast_outer_bounds_pp.map(pp).join(', ')}] pp。${s.assumed.contrast_outer_bounds_pp[0]<=0?'充分排序保证在此假设下不成立；这不是观测到物理翻转。':'仅在该假设上界确实成立时，外包络保持正值。'}`:'未指定实际 B。图展示不同假设的数学边界，没有对真实缺失贡献赋零，也没有从抽稀估计它。';
  $('study-support-scope').textContent=data.support.scope;
  const domain=data.support.known_o3_domain?.[0];$('study-tail-note').textContent=domain?`固定目标在 O3 红移支持外的质量约 ${(100*domain.population_redshift_tail_mass).toFixed(4)}%。只用检出概率 ≤ 1 得到的已知遗漏域上界约 ${domain.o3_time_weighted_selection_contribution_upper.toFixed(5)}，无法据此认证选择遗漏很小。`:'未知贡献需要额外似然或检出信息；有限样本重复不能发现集合中从未出现的贡献。';
  $('study-plan-note').textContent=data.support.tail_plan?`后续资源规划：${data.support.tail_plan.direct_draws.toLocaleString()} 条直接注入。${data.support.tail_plan.description} 该规划尚未执行。`:'';
  $('study-limitations').textContent=m.limitations;
  if(m.id===data.support.completion_metric_id){const event=data.support.completion_event;$('study-completion-event').textContent=`${event.event} 的 ${event.rows.toLocaleString()} 行在这些模型上没有负分量权重；这不证明物理似然为零。`;
   const completion=data.support.completions.find(x=>x.eta_label===s.settings.completion_eta),r=completion?.results[s.settings.policy_id]?.[s.model.id];
   $('study-completion-result').textContent=r?`仅该声明情景：${m.symbol} 中位数 ${pct(r.quantiles[1])}，相对基线移动 ${pp(r.median_shift_pp)} pp；情景内 90% 区间 ${pct(r.quantiles[0])}–${pct(r.quantiles[2])}。η 未被样本估计。`:'尚未选择补入量；真实遗漏贡献仍未知。';
  }
  table('study-model-table',['条件模型','中位数','模型内 90% 区间'],m.models.map(model=>{const button=document.createElement('button');button.type='button';button.className='study-model-choice';button.textContent=model.label;button.setAttribute('aria-pressed',String(model.id===s.model.id));button.onclick=()=>{populateModels(model.id);render()};const q=model.results[s.settings.policy_id].quantiles;return [button,pct(q[1]),`${pct(q[0])}–${pct(q[2])}`]}));
 }
 function metadata(){
  $('population-study-title').textContent=data.title;$('population-study-count').textContent=`${data.observations.events} 事件 · ${data.observations.posterior_rows.toLocaleString()} 后验行`;
  $('study-origin').textContent=imported?'本地导入 · 来源声明须另行核查':'全来源研究快照 · v0.3';
  $('study-assumptions').replaceChildren();for(const text of data.fixed_assumptions){const li=document.createElement('li');li.textContent=text;$('study-assumptions').append(li)}
  $('study-validation').textContent=data.validation?.description||'验证范围以数据包中的来源记录为准。';
  $('study-source-links').replaceChildren();for(const source of data.provenance.public_sources||[]){if(!/^https:\/\//.test(source.url))continue;const a=document.createElement('a');a.href=source.url;a.textContent=source.label;a.target='_blank';a.rel='noopener';$('study-source-links').append(a)}
  $('study-papers').replaceChildren();for(const paper of data.papers||[]){const a=document.createElement('a');a.href=paper.url;a.textContent=`${paper.language==='zh'?'中文':'English'} LaTeX PDF v${paper.version} · ${paper.pages} 页`;$('study-papers').append(a)}
  $('study-source').textContent=JSON.stringify({dataset_id:data.id,data_sha256:contentHash,hash_scope:hashScope,imported,declared_origin_sha256:declaredOriginHash,scientific_status:data.scientific_status,provenance:data.provenance,validation:data.validation},null,2);
 }
 function apply(study,hash,{isImport=false,viewSettings=null,scope='fetched_json_bytes',originHash=null}={}){
  validatePopulationStudy(study);const state=validateStudySettings(study,viewSettings||study.defaults);
  data=study;contentHash=hash;imported=isImport;hashScope=scope;declaredOriginHash=originHash;
  fillSelect('study-metric',data.metrics,state.metric_id);populateModels(state.model_id);fillSelect('study-policy',data.policies,state.policy_id);populateCompletion(state.completion_eta);
  $('study-budget-enable').checked=state.assumed_B!==null;$('study-budget').value=state.assumed_B??.5;
  $('study-reset').disabled=$('study-export').disabled=false;metadata();render();
 }
 for(const id of ['study-policy','study-shape','study-budget','study-budget-enable','study-completion-eta'])$(id).addEventListener('input',render);
 $('study-metric').onchange=()=>{populateModels();populateCompletion();render()};
 $('study-magnitude').onchange=()=>{const m=metric(),shape=$('study-shape').value;fillSelect('study-shape',m.shapes.filter(x=>m.models.some(r=>r.magnitude_id===$('study-magnitude').value&&r.shape_id===x.id)),shape);render()};
 $('study-reset').onclick=()=>{if(data)apply(data,contentHash,{isImport:imported,scope:hashScope,originHash:declaredOriginHash})};
 async function exportState(){if(!data)return null;const session=populationStudySession(data,settings(),contentHash);session.hash_scope=hashScope;session.declared_origin_sha256=declaredOriginHash;session.embedded_study_sha256=await digest(new TextEncoder().encode(JSON.stringify(data)));return session}
 $('study-export').onclick=async()=>{try{const session=await exportState();if(!session)return;const url=URL.createObjectURL(new Blob([JSON.stringify(session,null,2)],{type:'application/json'})),a=document.createElement('a');a.href=url;a.download='population-full-source-study.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000)}catch(e){$('population-study-status').textContent='未导出：'+e.message}};
 $('study-import').onclick=()=>$('study-file').click();
 $('study-file').onchange=async()=>{const file=$('study-file').files[0];if(!file)return;const token=++loadToken;try{
  if(file.size>12*1024*1024)throw Error('研究快照超过 12 MB');const parsed=JSON.parse(await file.text()),session=parsed.schema==='blackhole-population-study-session/1',study=session?parsed.study:parsed;
  const hash=await digest(new TextEncoder().encode(JSON.stringify(study)));if(session&&parsed.embedded_study_sha256&&parsed.embedded_study_sha256!==hash)throw Error('内嵌研究内容与其校验值不一致');
  if(token!==loadToken)return;apply(study,hash,{isImport:true,viewSettings:session?parsed.settings:null,scope:'JSON.stringify_embedded_study',originHash:session?parsed.dataset_sha256:null});
  $('population-study-status').textContent='已导入本地快照：格式与曲线一致性检查通过，科学来源和校准声明仍须核查。';
 }catch(e){$('population-study-status').textContent='未导入，保留当前研究：'+e.message}finally{$('study-file').value=''}};
 async function load(){const token=++loadToken;try{const r=await fetch(url);if(!r.ok)throw Error('全来源研究快照尚未载入');const bytes=await r.arrayBuffer();if(bytes.byteLength>12*1024*1024)throw Error('研究快照过大');const hash=await digest(bytes),study=JSON.parse(new TextDecoder().decode(bytes));if(token!==loadToken)return false;apply(study,hash);$('population-study-status').textContent='已载入全来源条件研究。浏览器展示已保存计算及声明的假设探索；实际遗漏预算与形成渠道仍未识别。';return true}catch(e){if(token===loadToken)$('population-study-status').textContent='全来源研究未载入：'+e.message+'。历史快照可单独展开查看。';return false}}
 const ready=load();return {ready,render,exportState,get data(){return data},get contentHash(){return contentHash},get settings(){return data?settings():null}};
}
