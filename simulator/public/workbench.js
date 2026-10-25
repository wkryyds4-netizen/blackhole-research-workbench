import {createWaveformView} from './waveform-view.mjs';
import {validateScene,sampleScene,opticalCapability,opticalFrame,singleScene,orbitExample} from './scene-contract.mjs';
import {posterior} from './population-contract.mjs';
import {createPopulationStudyView} from './population-study-view.mjs';
import {createAlignmentView} from './alignment-view.mjs';
import {createDynamicsView} from './dynamics-view.mjs';
const $=id=>document.getElementById(id),TAU=Math.PI*2;
let scene=validateScene(singleScene()),sourceMode='kerr',time=0,playing=false,raf=null,last=0,view='scene',lensing=null,population=null,projects=[],az=.4,el=.5;
let lastFrame=null,capability=null,opticalFailure='',datasetCatalog=[],loadToken=0,sceneExtent=5,sceneSeries=[];
function download(data,name){const url=URL.createObjectURL(new Blob([JSON.stringify(data,null,2)],{type:'application/json'}));const a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
function canvas(id){const e=$(id),r=e.getBoundingClientRect();if(r.width<1||r.height<1)return null;const d=Math.min(devicePixelRatio||1,2);if(e.width!==Math.round(r.width*d)||e.height!==Math.round(r.height*d)){e.width=Math.round(r.width*d);e.height=Math.round(r.height*d);}const c=e.getContext('2d');c.setTransform(d,0,0,d,0,0);c.clearRect(0,0,r.width,r.height);return {c,w:r.width,h:r.height};}
function coordinates(f){const d=canvas('generic-coordinates');if(!d)return;const {c,w,h}=d;const extent=sceneExtent,scale=Math.min(w,h)*.35/extent;
 const proj=p=>{const x=p[0]*Math.cos(az)-p[2]*Math.sin(az),z=p[0]*Math.sin(az)+p[2]*Math.cos(az);return [w/2+x*scale,h/2+(z*Math.sin(el)-p[1]*Math.cos(el))*scale]};
 c.strokeStyle='#274252';c.lineWidth=1;for(let i=-2;i<=2;i++){const a=proj([extent*i/2,0,-extent]),b=proj([extent*i/2,0,extent]),a2=proj([-extent,0,extent*i/2]),b2=proj([extent,0,extent*i/2]);c.beginPath();c.moveTo(...a);c.lineTo(...b);c.moveTo(...a2);c.lineTo(...b2);c.stroke();}
 const cols=['#62d9db','#eeb977','#b7a2ef'];
 for(const [axis,v] of [['x',[extent,0,0]],['y',[0,extent,0]],['z',[0,0,extent]]]){const o=proj([0,0,0]),p=proj(v);c.strokeStyle='#496375';c.beginPath();c.moveTo(...o);c.lineTo(...p);c.stroke();c.fillStyle='#94a7b9';c.fillText(axis,...p);}
 f.bodies.forEach((body,i)=>{c.strokeStyle=cols[i%3]+'70';c.lineWidth=1;c.beginPath();let begun=false;const window=(scene.frames.at(-1).t-scene.frames[0].t)*.04;for(const fr of scene.frames){if(fr.t>f.t)break;if(fr.t<f.t-window)continue;const b=fr.bodies.find(b=>b.id===body.id);if(!b)continue;const p=proj(b.position);if(begun)c.lineTo(...p);else{c.moveTo(...p);begun=true}}c.stroke()});
 f.bodies.forEach((b,i)=>{const col=cols[i%3],p=proj(b.position),r=Math.max(5,Math.sqrt(b.mass)*16);c.beginPath();c.fillStyle='#010408';c.strokeStyle=col;c.lineWidth=2;c.arc(...p,r,0,TAU);c.fill();c.stroke();c.fillStyle=col;c.font='12px system-ui';c.textAlign='center';c.fillText(b.id,p[0],p[1]+r+17);if(Math.hypot(...b.spin)>1e-8){const s=proj(b.position.map((v,k)=>v+b.spin[k]*extent*.7));c.beginPath();c.moveTo(...p);c.lineTo(...s);c.stroke();c.beginPath();c.arc(...s,3,0,TAU);c.fill();}});
}
function renderScene(){
 if(view!=='scene')return;
 const f=sampleScene(scene,time),cap=opticalCapability(f);lastFrame=f;capability=cap;
 const use=cap.supported&&lensing&&$('generic-mode').value==='auto';
 $('generic-optics').hidden=!use;$('generic-coordinates').hidden=!!use;lensing?.setActive(!!use);
 if(use)lensing.setFrame(opticalFrame(f),{playing,clockOrigin:scene.frames[0].t});else coordinates(f);
 $('capability').textContent=use?cap.mode:'坐标与自旋';
 $('capability-reason').textContent=use?(f.bodies.length===1?(scene.provenance.kind==='numerical'?'将该帧质量、自旋代入定常 Kerr 参考图像；共同视界刚形成时尚非定常，图像不是该模拟的真实时变光学。':'采用单黑洞定常度规进行数值光线积分；盘的发光模型仍是示意。'):'两中心光学近似，仅回放所给轨迹；未求解时变双黑洞度规。'):(opticalFailure||($('generic-mode').value==='coordinates'?'已选择坐标视图：半径与自旋箭头经显示缩放。':cap.reason||'透镜加载中，暂用坐标视图。'));
 const seconds=scene.units.time==='s'?time:time*scene.units.reference_mass_msun*4.9254909476412675e-6;
 $('generic-clock').textContent=`t = ${time.toFixed(3)} ${scene.units.time} · ${seconds.toExponential(3)} s`;
 $('generic-caption').textContent=use?'拖动旋转 · 滚轮缩放 · 环境盘为示意':'拖动旋转 · 线性显示插值 · 坐标轨迹不等于光学观测';
 $('body-count').textContent=f.bodies.length;$('time-unit').textContent=scene.units.time;$('reference-mass').textContent=scene.units.reference_mass_msun+' M☉';
 $('generic-source').textContent=`${scene.title}。${scene.provenance.kind}。来源：${scene.provenance.source}。版本：${scene.provenance.version}。坐标：${scene.provenance.frame}。质量参考系：${scene.units.mass_frame}。`;
 $('generic-source').textContent += [scene.provenance.clock,scene.provenance.physical_scale,scene.provenance.license,scene.provenance.archive_note].filter(Boolean).join(' · ');
 $('generic-merger').hidden=scene.provenance.source_common_horizon_time===undefined;
 drawSceneSeries(f);
 const duration=scene.frames.at(-1).t-scene.frames[0].t;
 $('generic-scrub').disabled=duration===0;$('generic-play').disabled=duration===0;$('generic-scrub').value=duration?(time-scene.frames[0].t)/duration*1000:0;
}
function stop(){playing=false;cancelAnimationFrame(raf);raf=null;$('generic-play').textContent='▶ 播放';}
function drawSceneSeries(f){
 const quantity=$('scene-quantity').value,d=canvas('scene-series');if(!d)return;const {c,w,h}=d,L=50,R=16,T=28,B=35,pw=w-L-R,ph=h-T-B;
 const ids=[...new Set(sceneSeries.flatMap(x=>x.bodies.map(b=>b.id)))],colors=['#62d9db','#eeb977','#b7a2ef'];
 const value=(r,id)=>quantity==='separation'?r.separation:(()=>{const b=r.bodies.find(b=>b.id===id);return b?(quantity==='mass'?b.mass:Math.hypot(...b.spin)):null})();
 const lines=quantity==='separation'?['A–B']:ids;let ymax=0;for(const row of sceneSeries)for(const id of lines){const v=value(row,id);if(v!==null)ymax=Math.max(ymax,v)}ymax=ymax>0?ymax*1.1:1;
 const start=scene.frames[0].t,end=scene.frames.at(-1).t,X=t=>L+(t-start)/(end-start||1)*pw,Y=v=>T+ph*(1-v/ymax);
 c.font='11px system-ui';for(let j=0;j<=4;j++){const y=T+ph*j/4;c.strokeStyle='#253345';c.beginPath();c.moveTo(L,y);c.lineTo(w-R,y);c.stroke();c.fillStyle='#94a7b9';c.textAlign='right';c.fillText((ymax*(1-j/4)).toPrecision(3),L-5,y+4);c.textAlign='center';c.fillText((start+(end-start)*j/4).toFixed(0),L+pw*j/4,h-18)}
 lines.forEach((id,i)=>{c.strokeStyle=colors[i%3];c.fillStyle=c.strokeStyle;c.textAlign='left';c.fillText(id,L+i*65,15);c.lineWidth=1.5;c.beginPath();let active=false;for(const row of sceneSeries){const v=value(row,id);if(v===null){active=false;continue}if(active)c.lineTo(X(row.t),Y(v));else c.moveTo(X(row.t),Y(v));active=true}c.stroke()});
 if(scene.provenance.source_common_horizon_time!==undefined){c.strokeStyle='#526879';c.setLineDash([3,3]);c.beginPath();c.moveTo(X(0),T);c.lineTo(X(0),T+ph);c.stroke();c.setLineDash([])}
 c.strokeStyle='#ffffff';c.beginPath();c.moveTo(X(f.t),T);c.lineTo(X(f.t),T+ph);c.stroke();c.fillStyle='#94a7b9';c.textAlign='center';c.fillText('t / '+scene.units.time,L+pw/2,h-2);
 const sep=f.bodies.length===2?Math.hypot(...f.bodies[0].position.map((v,k)=>v-f.bodies[1].position[k])):null;
 $('scene-values').textContent=f.bodies.map(b=>`${b.id}: m/Mref=${b.mass.toFixed(6)}, |χ|=${Math.hypot(...b.spin).toFixed(6)}`).join(' · ')+(sep===null?' · 当前没有双体分离度':` · d/Mref=${sep.toFixed(4)}`);
 const p=scene.provenance,gap=p.gaps?.find(g=>f.t>g.start&&f.t<g.end);
 $('scene-data-note').textContent=(p.deprecated?`旧版算例，官方已由 ${p.superseded_by} 替代；用于接口验证。`:'')+(gap?'当前位于原始采样缺口：保持最后 A/B 帧，未推算缺失轨迹。':p.source_common_horizon_time!==undefined?'t=0 为共同视界首个样本。曲线在 A/B→C 处断开；视界质量变化不是直接的辐射能量曲线。':p.kind==='numerical'?'数值来自当前场景；坐标、时间原点与物理解释以来源记录为准。':'数值来自当前场景；示例的质量和轨道并非相对论演化解。');
}
$('scene-quantity').onchange=()=>drawSceneSeries(lastFrame);
$('generic-merger').onclick=()=>{stop();time=0;renderScene()};
$('scene-series').onclick=e=>{const r=e.target.getBoundingClientRect(),u=Math.max(0,Math.min(1,(e.clientX-r.left-50)/(r.width-66)));stop();time=scene.frames[0].t+u*(scene.frames.at(-1).t-scene.frames[0].t);renderScene()};
async function loadDatasets(){try{let r=await fetch('datasets/catalog.json');if(r.status===404)r=await fetch('datasets/catalog.public.json');if(!r.ok)throw Error('数据目录读取失败');datasetCatalog=await r.json();const group=document.createElement('optgroup');group.label='SXS 数值相对论 · 可用快照';for(const entry of datasetCatalog){const o=document.createElement('option');o.value=entry.id;o.textContent=entry.title;group.append(o)}$('generic-preset').append(group);const requested=new URLSearchParams(location.search).get('dataset');if(datasetCatalog.some(d=>d.id===requested)){$('generic-preset').value=requested;$('generic-preset').dispatchEvent(new Event('change'))}}catch(e){$('import-status').textContent=e.message}}
function cacheScene(){sceneExtent=5;sceneSeries=scene.frames.map(f=>{for(const b of f.bodies)sceneExtent=Math.max(sceneExtent,Math.hypot(...b.position)+2*b.mass);return {t:f.t,bodies:f.bodies,separation:f.bodies.length===2?Math.hypot(...f.bodies[0].position.map((v,k)=>v-f.bodies[1].position[k])):null}})}
function loadScene(data,mode){const valid=validateScene(data);loadToken++;stop();scene=valid;sourceMode=mode;time=scene.frames[0].t;cacheScene();$('generic-spin').disabled=!['kerr','schwarzschild'].includes(mode);$('generic-spin-value').textContent=$('generic-spin').disabled?'由数据提供':(+$('generic-spin').value).toFixed(2);renderScene();}
function tick(now){if(!playing||view!=='scene')return;const dt=last?Math.min(.1,(now-last)/1000):0;last=now;const duration=scene.frames.at(-1).t-scene.frames[0].t;time+=dt*duration/(+$('generic-duration').value);if(time>=scene.frames.at(-1).t){time=scene.frames.at(-1).t;stop();}renderScene();if(playing)raf=requestAnimationFrame(tick);}
$('generic-play').onclick=()=>{if(playing){stop();renderScene();return;}if(time>=scene.frames.at(-1).t)time=scene.frames[0].t;playing=true;last=0;$('generic-play').textContent='Ⅱ 暂停';raf=requestAnimationFrame(tick)};
$('generic-scrub').oninput=()=>{stop();time=scene.frames[0].t+$('generic-scrub').value/1000*(scene.frames.at(-1).t-scene.frames[0].t);renderScene()};
$('generic-preset').onchange=async()=>{const mode=$('generic-preset').value;const entry=datasetCatalog.find(e=>e.id===mode);if(entry){const token=++loadToken;stop();$('import-status').textContent='正在载入已校验的本地数据…';try{const r=await fetch(entry.url);if(!r.ok)throw Error('数据读取失败');const data=await r.json();if(token!==loadToken)return;loadScene(data,mode);$('import-status').textContent='已载入数值相对论轨迹。参考质量 65 M☉ 仅用于显示换算。'}catch(e){if(token===loadToken){$('generic-preset').value=sourceMode;$('import-status').textContent=e.message}}return;}if(mode==='kerr'||mode==='schwarzschild'){$('generic-spin').value=mode==='kerr'?.7:0;$('generic-spin-value').textContent=(+$('generic-spin').value).toFixed(2);loadScene(singleScene(+$('generic-spin').value),mode)}else loadScene(orbitExample(mode==='spin'),mode);$('import-status').textContent='已载入内置示例。'};
$('generic-spin').oninput=()=>{$('generic-spin-value').textContent=(+$('generic-spin').value).toFixed(2);loadScene(singleScene(+$('generic-spin').value),sourceMode)};
$('generic-mode').onchange=renderScene;$('generic-quality').onchange=()=>lensing?.setQuality($('generic-quality').value);$('generic-disk').onchange=()=>lensing?.setDisk($('generic-disk').checked);
$('generic-oblique').onclick=()=>{az=.4;el=.5;lensing?.cameraView(false);renderScene()};$('generic-top').onclick=()=>{el=Math.PI/2;lensing?.cameraView(true);renderScene()};
let drag=null;$('generic-coordinates').onpointerdown=e=>{drag={x:e.clientX,y:e.clientY,az,el};e.target.setPointerCapture(e.pointerId)};$('generic-coordinates').onpointermove=e=>{if(drag){az=drag.az+(e.clientX-drag.x)*.007;el=Math.max(.1,Math.min(1.57,drag.el+(e.clientY-drag.y)*.007));coordinates(lastFrame)}};$('generic-coordinates').onpointerup=$('generic-coordinates').onpointercancel=()=>drag=null;
$('generic-export').onclick=()=>download(scene,'blackhole-scene.json');$('generic-import').onclick=()=>$('generic-file').click();
$('generic-file').onchange=async()=>{const f=$('generic-file').files[0];if(!f)return;try{if(f.size>5_000_000)throw Error('文件上限为 5 MB');loadScene(JSON.parse(await f.text()),'imported');$('generic-preset').value='';$('import-status').textContent='已载入：'+scene.title+'。单位与来源通过格式检查；数据物理有效性仍由提供者负责。';}catch(e){$('import-status').textContent='未导入，保留当前场景：'+e.message}finally{$('generic-file').value=''}};
function switchView(v){stop();dynamics.stop();view=v;lensing?.setActive(v==='scene');for(const name of ['scene','alignment','waveform','dynamics','population','projects'])$(name+'-view').hidden=name!==v;document.querySelectorAll('[data-view]').forEach(b=>{b.classList.toggle('active',b.dataset.view===v);b.setAttribute('aria-pressed',b.dataset.view===v?'true':'false')});if(v==='scene')renderScene();if(v==='population'){drawPopulation();populationStudy.render()}if(v==='alignment')alignment.render();if(v==='dynamics')dynamics.render();if(v==='waveform')waveform.render();}
document.querySelectorAll('[data-view]').forEach(b=>b.onclick=()=>switchView(b.dataset.view));
function chart(id,series,{ymax,ylabel,threshold}={}){const d=canvas(id);if(!d)return;const {c,w,h}=d,L=46,R=16,T=26,B=36,pw=w-L-R,ph=h-T-B;let max=ymax||Math.max(...series.flatMap(s=>s.y.filter(Number.isFinite)))*1.12;max=Math.max(max,1e-6);const X=x=>L+x*pw,Y=y=>T+ph-y/max*ph;c.font='11px system-ui';for(let i=0;i<=4;i++){const x=L+pw*i/4,y=T+ph*i/4;c.strokeStyle='#253345';c.lineWidth=.7;c.beginPath();c.moveTo(L,y);c.lineTo(w-R,y);c.stroke();c.fillStyle='#94a7b9';c.textAlign='right';c.fillText((max*(1-i/4)).toFixed(max>10?0:1),L-6,y+4);c.textAlign='center';c.fillText((i/4).toFixed(2),x,h-18);}c.textAlign='left';c.fillText(ylabel||'概率密度',L,13);c.textAlign='center';c.fillText('偏对齐代理分量比例 f',L+pw/2,h-2);c.save();c.beginPath();c.rect(L,T,pw,ph);c.clip();for(const s of series){c.strokeStyle=s.color;c.lineWidth=1.7;c.setLineDash(s.dashed?[4,3]:[]);c.beginPath();s.y.forEach((v,i)=>{const x=X(population.f_grid[i]),y=Y(v);i?c.lineTo(x,y):c.moveTo(x,y)});c.stroke();}if(threshold!==undefined){c.strokeStyle='#eeb977';c.setLineDash([5,4]);c.beginPath();c.moveTo(L,Y(threshold));c.lineTo(w-R,Y(threshold));c.stroke();}c.restore();}
function popSettings(){return {key:['tilt','spin','mass','redshift'].map(x=>$('pop-'+x).value).join('-'),prior:$('pop-prior').value,calibration:+$('pop-calibration').value};}
let scenarioPoints=[];
function drawPopulation(){if(!population||view!=='population')return;const settings=popSettings(),r=posterior(population,settings.key,settings),uncorrected=posterior(population,settings.key,{...settings,selection:false});
 $('pop-calibration-value').textContent=settings.calibration.toFixed(2);
 let prior=population.f_grid.map(f=>settings.prior==='central'?f*(1-f):settings.prior==='edge'?1/Math.sqrt(f*(1-f)):1);const norm=prior.reduce((a,b)=>a+b,0)/prior.length;prior=prior.map(x=>x/norm);
 chart('population-posterior',[{y:r.density,color:'#62d9db'},{y:uncorrected.density,color:'#eeb977',dashed:true},{y:prior,color:'#b7a2ef',dashed:true}]);
 const fmt=x=>(x*100).toFixed(1)+'%';$('pop-interval').textContent=`当前条件中位数 ${fmt(r.interval[1])}；90% 等尾区间 ${fmt(r.interval[0])}–${fmt(r.interval[2])}。忽略选择的对照中位数 ${fmt(uncorrected.interval[1])}。`;
 chart('population-mc',[{y:r.mc.map(x=>x.sd),color:'#62d9db'}],{ylabel:'σMC(log L)',threshold:1});
 $('pop-coverage').textContent=`当前后验质量中 ${fmt(r.passing)} 通过三项 MC 筛查。这个比例衡量积分精度覆盖，不能用于模型证据或渠道判别。`;
 const eps=.001,up=posterior(population,settings.key,{...settings,calibration:settings.calibration*Math.exp(eps)}),down=posterior(population,settings.key,{...settings,calibration:settings.calibration*Math.exp(-eps)}),slope=(up.mean-down.mean)/(2*eps);
 $('pop-sensitivity').textContent=`局部选择敏感性：ξA/ξI 增加约 1%，后验均值约改变 ${(slope).toFixed(2)} 个百分点（其他假设固定，一阶近似）。k 是人为压力测试，不是校准误差的测量值。`;
 const d=canvas('population-scenarios');if(!d)return;const {c,w,h}=d,L=42,R=15,T=24,B=35;c.font='10px system-ui';c.fillStyle='#94a7b9';c.textAlign='left';c.fillText('中位数 f',L,13);for(let j=0;j<=4;j++){const y=T+(h-T-B)*(1-j/4);c.strokeStyle='#253345';c.beginPath();c.moveTo(L,y);c.lineTo(w-R,y);c.stroke();c.textAlign='right';c.fillText((j/4).toFixed(2),L-7,y+4);}scenarioPoints=[];
 Object.keys(population.scenarios).forEach((key,i)=>{const out=posterior(population,key,settings),x=L+i/80*(w-L-R),y=T+(1-out.interval[1])*(h-T-B),selected=key===settings.key;c.strokeStyle=selected?'#eeb977':'#62d9db';c.fillStyle=c.strokeStyle;c.lineWidth=selected?2:1;c.beginPath();c.arc(x,y,selected?5:2.7,0,TAU);if(out.passing>=.9)c.fill();c.stroke();scenarioPoints.push({key,x,y});});c.fillStyle='#94a7b9';c.textAlign='center';c.fillText('固定假设组合 1 → 81（非后验样本）',L+(w-L-R)/2,h-8);
}
for(const id of ['pop-tilt','pop-spin','pop-mass','pop-redshift','pop-prior','pop-calibration'])$(id).addEventListener('input',drawPopulation);
$('population-scenarios').onclick=e=>{if(!scenarioPoints.length)return;const rect=e.target.getBoundingClientRect(),x=e.clientX-rect.left,y=e.clientY-rect.top;const nearest=scenarioPoints.reduce((a,b)=>Math.hypot(x-a.x,y-a.y)<Math.hypot(x-b.x,y-b.y)?a:b);nearest.key.split('-').forEach((v,i)=>$('pop-'+['tilt','spin','mass','redshift'][i]).value=v);drawPopulation()};
$('pop-reset').onclick=()=>{['1','1','0','1'].forEach((v,i)=>$('pop-'+['tilt','spin','mass','redshift'][i]).value=v);$('pop-prior').value='uniform';$('pop-calibration').value='1';drawPopulation()};
$('pop-export').onclick=()=>{if(!population)return;const s=popSettings();download({schema:'population-session/1',input_sha256:population.input_sha256,settings:s,result:posterior(population,s.key,s),interpretation:population.interpretation,selection_audit:population.selection_audit},'population-research-state.json')};
async function loadResearch(){try{const r=await fetch('research-data.json');if(!r.ok)throw Error('研究数据加载失败');population=await r.json();$('population-status').textContent='已载入固定数据快照；推断与诊断在本地进行。';drawPopulation()}catch(e){$('population-status').textContent=e.message}}
async function loadProjects(){try{const r=await fetch('projects.json');projects=await r.json();for(const p of projects){const card=document.createElement('article');card.className='project-card';const h=document.createElement('h2'),a=document.createElement('a');a.href=p.url;a.target='_blank';a.rel='noopener';a.textContent=p.name;h.append(a);card.append(h);for(const [text,cls] of [[p.role,''],[p.action,'project-state'],[p.limits,''],[`许可：${p.license} · ${p.status}`,'']]){const el=document.createElement('p');el.textContent=text;el.className=cls;card.append(el)}$('project-cards').append(card)}}catch(e){$('project-cards').textContent='项目记录加载失败：'+e.message}}
async function initOptics(){try{const {createLensing}=await import('./lensing/renderer.js');lensing=createLensing($('generic-optics'),{onError:e=>{lensing=null;opticalFailure=e.message;renderScene()}});lensing.setQuality($('generic-quality').value);lensing.setDisk($('generic-disk').checked);lensing.setActive(view==='scene');renderScene()}catch(e){opticalFailure='光学不可用，使用坐标视图：'+e.message;renderScene()}}
const alignment=createAlignmentView();
const dynamics=createDynamicsView();
const waveform=createWaveformView();
const populationStudy=createPopulationStudyView();
$('population-legacy').addEventListener('toggle',drawPopulation);
const resize=new ResizeObserver(()=>{if(view==='scene')renderScene();else if(view==='population'){drawPopulation();populationStudy.render()}else if(view==='alignment')alignment.render();else if(view==='dynamics')dynamics.render();else if(view==='waveform')waveform.render()});document.querySelectorAll('canvas').forEach(c=>resize.observe(c));
document.addEventListener('visibilitychange',()=>{if(document.hidden){stop();lensing?.setActive(false)}else{if(view==='scene')renderScene()}});
window.researchWorkbench={get scene(){return scene},get frame(){return lastFrame},get capability(){return capability},get lensing(){return lensing},get population(){return population},populationStudy,alignment,dynamics,waveform,loadScene,switchView,posterior,popSettings};
cacheScene();renderScene();loadResearch();loadProjects();loadDatasets();initOptics();const initialView=new URLSearchParams(location.search).get('view');if(['alignment','dynamics','waveform','population','projects'].includes(initialView))switchView(initialView);
