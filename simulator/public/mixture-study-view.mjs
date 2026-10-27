import {validateMixtureStudy,selectMixtureCase} from './mixture-study-contract.mjs';
const pct=x=>(100*x).toFixed(2)+'%';
const element=(tag,text,cls)=>{const e=document.createElement(tag);if(text!==undefined)e.textContent=text;if(cls)e.className=cls;return e};
export function createMixtureStudyView(container,url='datasets/channel-study.json'){
 let data,settings,quantity,sourceHash;const selectors=new Map();
 const status=element('p','正在读取物理模型条件研究…','notice');status.setAttribute('role','status');container.append(status);
 function render(){
  if(!data)return;const r=selectMixtureCase(data,settings),q=data.quantities.find(q=>q.id===quantity),interval=Math.round(100*(data.quantile_probabilities[2]-data.quantile_probabilities[0]));
  definition.textContent=q.definition;context.textContent=data.context+' '+data.domain;
  chart.replaceChildren();const ns='http://www.w3.org/2000/svg';const svg=document.createElementNS(ns,'svg');svg.setAttribute('viewBox',`0 0 760 ${data.components.length*42+50}`);svg.setAttribute('role','img');svg.setAttribute('aria-label',q.label+`：模型内${interval}%可信区间`);
  function draw(tag,attrs,text){const n=document.createElementNS(ns,tag);for(const [k,v] of Object.entries(attrs))n.setAttribute(k,String(v));if(text!==undefined)n.textContent=text;svg.append(n)}
  const X=x=>180+540*x;
  for(let tick=0;tick<=100;tick+=20){const x=X(tick/100);draw('line',{x1:x,x2:x,y1:10,y2:data.components.length*42+10,stroke:'#293c4b'});draw('text',{x,y:data.components.length*42+35,fill:'#9aafbf','text-anchor':'middle','font-size':12},tick+'%')}
  const body=table.tBodies[0];body.replaceChildren();
  data.components.forEach((c,i)=>{const [lo,mid,hi]=r.quantiles[quantity][i],y=28+i*42;draw('text',{x:165,y:y+4,fill:'#e0eaf0','text-anchor':'end','font-size':13},c.label);draw('line',{x1:X(lo),x2:X(hi),y1:y,y2:y,stroke:'#65d5dc','stroke-width':3});draw('circle',{cx:X(mid),cy:y,r:5,fill:'#65d5dc'});const row=body.insertRow();for(const value of [c.label,pct(mid),pct(lo)+' – '+pct(hi)])row.insertCell().textContent=value});
  chart.append(svg);common.textContent=`共同正则项的检出混合份额：中位数 ${pct(r.common_term[1])}，${interval}% 区间 ${pct(r.common_term[0])}–${pct(r.common_term[2])}。${r.common_term_definition}`;
  status.textContent=`已计算的 ${data.cases.length} 组条件对照 · ${data.observations.events} 个事件 · ${data.observations.posterior_rows.toLocaleString()} 行。`;
 }
 let definition,context,chart,table,common;
 async function load(){
  try{
   const response=await fetch(url);if(!response.ok)throw Error('HTTP '+response.status);const raw=await response.text();if(raw.length>2000000)throw Error('数据过大');data=validateMixtureStudy(JSON.parse(raw));settings={...data.defaults};quantity=data.quantities[0].id;
   sourceHash=Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',new TextEncoder().encode(raw))),v=>v.toString(16).padStart(2,'0')).join('');
   const heading=element('h2',data.title);container.prepend(heading);container.classList.add('research-panel','mixture-study');
   container.append(element('p',data.interpretation,'intro'));context=element('p','','notice');container.append(context);
   const controls=element('div',undefined,'population-controls');
   for(const axis of data.axes){const label=element('label',axis.label),select=element('select');select.dataset.axis=axis.id;select.setAttribute('aria-label',axis.label);for(const option of axis.options){const o=element('option',option.label);o.value=option.id;select.append(o)}select.value=settings[axis.id];select.addEventListener('change',()=>{settings[axis.id]=select.value;render()});selectors.set(axis.id,select);label.append(select);controls.append(label)}
   const label=element('label','显示科学量'),select=element('select');select.setAttribute('aria-label','显示科学量');for(const q of data.quantities){const o=element('option',q.label);o.value=q.id;select.append(o)}select.addEventListener('change',()=>{quantity=select.value;render()});label.append(select);controls.append(label);container.append(controls);
   definition=element('p','','notice');chart=element('div',undefined,'mixture-chart');table=element('table');const head=table.createTHead().insertRow();for(const name of ['分量','边缘中位数',`模型内 ${Math.round(100*(data.quantile_probabilities[2]-data.quantile_probabilities[0]))}% 区间`]){const th=element('th',name);th.scope='col';head.append(th)}table.createTBody();const scroll=element('div',undefined,'table-scroll');scroll.append(table);common=element('p','','notice');container.append(definition,chart,scroll,common);
   const limits=element('ul');for(const text of data.limits)limits.append(element('li',text));container.append(limits);
   const links=element('p');for(const p of data.papers){const a=element('a',p.label);a.href=p.url;a.target='_blank';a.rel='noopener';links.append(a,document.createTextNode('  '))}container.append(links);
   const button=element('button','导出当前比较设置');button.type='button';button.addEventListener('click',()=>{const r=selectMixtureCase(data,settings),blob=new Blob([JSON.stringify({schema:'blackhole-mixture-view-session/1',study_id:data.id,data_sha256:sourceHash,settings,quantity,quantity_definition:data.quantities.find(q=>q.id===quantity).definition,scientific_status:data.scientific_status,value_unit:data.value_unit,quantile_probabilities:data.quantile_probabilities,domain:data.domain,context:data.context,interpretation:data.interpretation,limits:data.limits,result:r},null,2)],{type:'application/json'});const u=URL.createObjectURL(blob),a=element('a');a.href=u;a.download='channel-study-session.json';a.click();setTimeout(()=>URL.revokeObjectURL(u),1000)});container.append(button);render();
  }catch(error){status.textContent='条件研究加载失败：'+error.message;status.setAttribute('role','alert')}
 }
 const ready=load();return {render,ready,get data(){return data},get settings(){return {...settings}}};
}
