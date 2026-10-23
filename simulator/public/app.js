'use strict';
(() => {
const $=id=>document.getElementById(id), TAU=2*Math.PI;
const colors={text:'#b2c2d2',grid:'#253345',cyan:'#62d9db',gold:'#eeb977',purple:'#b7a2ef',gray:'#8195aa'};
const defaults={m1:36,m2:29,s1:0,s2:0,distance:440,inclination:35,model:'IMRPhenomD'};
const presets={default:defaults,equal:{...defaults,m1:30,m2:30,inclination:0},aligned:{...defaults,s1:.6,s2:.6},unequal:{...defaults,m1:60,m2:20,s1:.4,s2:-.2}};
const fields=[['m1','主黑洞质量 m₁',20,100,1,'M☉'],['m2','伴黑洞质量 m₂',10,100,1,'M☉'],['s1','主黑洞自旋 χ₁',-.8,.8,.05,''],['s2','伴黑洞自旋 χ₂',-.8,.8,.05,''],['distance','光度距离',100,2000,10,'Mpc'],['inclination','物理倾角 ι',0,180,1,'°']];
const state={p:{...defaults},wave:null,reference:null,time:0,playing:false,tab:'simulation',camera:{az:.3,el:.6},event:null,fit:null,busy:false};
let requestId=0,fitId=0,debounce,raf,last=0,audioContext,audioSource;
let lensing=null;
try {const s=JSON.parse(localStorage.getItem('bbh-lab-v1'));if(s?.parameters)state.p={...defaults,...s.parameters};}catch{}
for(const [id,label,min,max,step,unit] of fields){
 const block=document.createElement('div');block.innerHTML=`<label class="slider-label" for="${id}"><span>${label}</span><output id="${id}-value"></output></label><input type="range" id="${id}" min="${min}" max="${max}" step="${step}" aria-label="${label}">`;
 $('parameter-fields').append(block);
 $(id).addEventListener('input',()=>{state.p[id]=+$(id).value;if(id==='m1')state.p.m2=Math.max(Math.ceil(state.p.m1/4),Math.min(state.p.m1,state.p.m2));if(id==='m2')state.p.m1=Math.min(100,Math.max(state.p.m2,Math.min(state.p.m1,state.p.m2*4)));requestId++;fitId++;state.fit=null;state.busy=true;$('play').disabled=true;syncInputs();pause();clearTimeout(debounce);setStatus('参数已修改，准备计算',true);debounce=setTimeout(updateWave,350);});
}
function syncInputs(){for(const [id,, , ,step,unit] of fields){$(id).value=state.p[id];$(id+'-value').textContent=(step<1?state.p[id].toFixed(2):Math.round(state.p[id]))+' '+unit;}$('model').value=state.p.model;}
function setStatus(s,busy=false){$('status').textContent=s;$('status').classList.toggle('busy',busy);}
function stopAudio(){if(audioSource){try{audioSource.stop();}catch{}audioSource=null;}$('audio').textContent='♫ 听波形 · 原速';}
function pause(){state.playing=false;$('play').textContent='▶ 播放';stopAudio();if(raf)cancelAnimationFrame(raf);raf=null;lensing?.update(state);}
function save(){try{localStorage.setItem('bbh-lab-v1',JSON.stringify({parameters:state.p}));}catch{}}
async function post(path,data){const r=await fetch(path,{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(data)});const out=await r.json();if(!r.ok)throw Error(out.error||'计算失败');return out;}
async function updateWave(){
 const id=++requestId;fitId++;state.fit=null;pause();state.busy=true;$('play').disabled=true;$('error').hidden=true;setStatus('正在生成物理波形…',true);
 try{
  const w=await post('/api/waveform',state.p);if(id!==requestId)return;
  state.wave=w;state.p={...w.parameters};state.time=w.t[0]+(w.t.at(-1)-w.t[0])*.68;state.busy=false;
  $('play').disabled=false;$('chirp').textContent=w.meta.chirp_mass.toFixed(2)+' M☉';$('chieff').textContent=w.meta.chi_eff.toFixed(3);$('q').textContent=w.meta.q.toFixed(3);$('peak').textContent=w.meta.peak_strain.toExponential(2);
  $('footer-model').textContent=`${state.p.model} · Ripple ${w.meta.version} · ${w.sample_rate} Hz`;$('engine-version').textContent=`当前引擎：Ripple ${w.meta.version}；波形模型：${state.p.model}。`;
  setStatus('计算完成 · 物理时间已同步');syncInputs();save();drawAll();$('fit-summary').textContent='参数已更新，点击对齐当前模板。';
 }catch(e){if(id!==requestId)return;state.busy=false;$('error').textContent=e.message;$('error').hidden=false;setStatus('计算未完成');}
}
function dimensions(canvas){const r=canvas.getBoundingClientRect();if(r.width<1||r.height<1)return null;const dpr=Math.min(devicePixelRatio||1,2);if(canvas.width!==Math.round(r.width*dpr)||canvas.height!==Math.round(r.height*dpr)){canvas.width=Math.round(r.width*dpr);canvas.height=Math.round(r.height*dpr);}const c=canvas.getContext('2d');c.setTransform(dpr,0,0,dpr,0,0);c.clearRect(0,0,r.width,r.height);return {c,w:r.width,h:r.height};}
function sampleIndex(t){const w=state.wave;return Math.max(0,Math.min(w.t.length-1,Math.round((t-w.t[0])*w.sample_rate)));}
function finite(v){return v!==null&&Number.isFinite(v);}
function drawScene(){
 const optical=!!lensing&&$('scene-mode').value==='lensing';
 $('lensing').hidden=!optical;$('orbit').hidden=optical;
 lensing?.setActive(optical&&state.tab==='simulation');
 if(optical)lensing.update(state);else drawOrbit();
 if(!state.wave)return;
 const w=state.wave,i=sampleIndex(state.time),t=w.t[i];
 $('stage').textContent=t<-.012?'旋近':t<.008?'并合附近':'振铃';
 $('scene-time').textContent=`t = ${t.toFixed(4)} s`;
 $('scene-frequency').textContent=finite(w.frequency[i])?`fGW ≈ ${w.frequency[i].toFixed(1)} Hz`:t<-.012?'fGW = —（低幅度区间）':'旋近频率读数已结束';
 $('scene-caption').textContent=optical?'拖动旋转 · 滚轮缩放 · 光学无自旋近似 · '+(t>=-.012?'并合/遗迹为视觉过渡':'轨道由波形相位驱动'):'拖动旋转视角 · 轨道为示意，非数值相对论轨迹';
}
function lensingError(error){
 lensing?.dispose();lensing=null;if(window.bbhLab)window.bbhLab.lensing=null;$('scene-mode').value='schematic';$('scene-mode').querySelector('[value="lensing"]').disabled=true;
 $('optical-status').textContent='透镜不可用，已切换轨道示意：'+error.message;drawScene();
}
async function initLensing(){
 try{const {createLensing}=await import('./lensing/renderer.js');lensing=createLensing($('lensing'),{onError:lensingError,onStatus:s=>$('optical-status').textContent=s});window.bbhLab.lensing=lensing;lensing.setQuality($('optical-quality').value);lensing.setDisk($('optical-disk').checked);drawScene();}
 catch(e){lensingError(e);}
}
$('scene-mode').onchange=drawScene;
$('optical-quality').onchange=()=>lensing?.setQuality($('optical-quality').value);
$('optical-disk').onchange=()=>lensing?.setDisk($('optical-disk').checked);
function drawOrbit(){
 const d=dimensions($('orbit'));if(!d)return;const {c,w,h}=d;const wave=state.wave;if(!wave){c.fillStyle=colors.text;c.textAlign='center';c.fillText('正在加载经过物理模型计算的波形',w/2,h/2);return;}
 const i=sampleIndex(state.time),t=wave.t[i],p=wave.parameters,phase=wave.phase[i]/2;
 let sep=wave.separation[i];if(!finite(sep)){const before=wave.separation.slice(0,i).findLast(finite);sep=before??12;}
 const base=Math.min(w*.34,h*.38),scale=base/12;
 const az=state.camera.az,el=state.camera.el;
 function project(x,y,z=0){const a=x*Math.cos(az)-y*Math.sin(az),b=x*Math.sin(az)+y*Math.cos(az);return [w*.5+a*scale,h*.52+b*Math.sin(el)*scale-z*Math.cos(el)*scale,b*Math.cos(el)+z*Math.sin(el)];}
 // Orbit plane is a coordinate guide, not a literal spacetime surface.
 c.lineWidth=.7;c.strokeStyle='#48697f33';
 for(let r=4;r<=24;r+=4){c.beginPath();for(let k=0;k<=100;k++){let q=project(r*Math.cos(k*TAU/100),r*Math.sin(k*TAU/100));k?c.lineTo(q[0],q[1]):c.moveTo(q[0],q[1]);}c.stroke();}
 for(let k=0;k<12;k++){const a=k*TAU/12,q=project(24*Math.cos(a),24*Math.sin(a));c.beginPath();c.moveTo(w/2,h*.52);c.lineTo(q[0],q[1]);c.stroke();}
 const blend=Math.max(0,Math.min(1,-t/.012)),a=Math.min(26,Math.max(5,sep))*blend;
 if(t<0){
  for(const [frac,col] of [[p.m2/(p.m1+p.m2),colors.cyan],[-p.m1/(p.m1+p.m2),colors.gold]]){
   c.beginPath();c.strokeStyle=col+'70';c.lineWidth=1.4;
   for(let k=0;k<=90;k++){const angle=phase-k/90*TAU*1.4,r=a*(1+k/250)*frac;const q=project(r*Math.cos(angle),r*Math.sin(angle));k?c.lineTo(q[0],q[1]):c.moveTo(q[0],q[1]);}c.stroke();
  }
  const holes=[{mass:p.m1,frac:p.m2/(p.m1+p.m2),spin:p.s1,col:colors.cyan,label:'m₁'},{mass:p.m2,frac:-p.m1/(p.m1+p.m2),spin:p.s2,col:colors.gold,label:'m₂'}].map(b=>({...b,pos:project(a*b.frac*Math.cos(phase),a*b.frac*Math.sin(phase))})).sort((x,y)=>x.pos[2]-y.pos[2]);
  holes.forEach(b=>hole(b.pos[0],b.pos[1],Math.max(10,base*.10)*Math.sqrt(b.mass/32),b.col,b.label,b.spin));
 }else{hole(w/2,h*.52,base*.16,colors.cyan,'并合遗迹',0);}
 function hole(x,y,r,col,label,spin){
  const glow=c.createRadialGradient(x,y,r*.6,x,y,r*2.6);glow.addColorStop(0,col+'60');glow.addColorStop(.5,col+'1a');glow.addColorStop(1,col+'00');c.fillStyle=glow;c.beginPath();c.arc(x,y,r*2.6,0,TAU);c.fill();
  const disk=c.createRadialGradient(x-r*.3,y-r*.4,0,x,y,r);disk.addColorStop(0,'#1a2c3a');disk.addColorStop(.55,'#060c15');disk.addColorStop(.93,'#02070c');disk.addColorStop(1,col);c.fillStyle=disk;c.beginPath();c.arc(x,y,r,0,TAU);c.fill();
  c.strokeStyle=col+'b0';c.lineWidth=1;c.beginPath();c.ellipse(x,y,r*1.25,r*.32,-.2,0,TAU);c.stroke();
  if(Math.abs(spin)>.01){const sy=y-Math.sign(spin)*(r+15+Math.abs(spin)*20);c.strokeStyle=col;c.beginPath();c.moveTo(x,y);c.lineTo(x,sy);c.lineTo(x-4,sy+Math.sign(spin)*7);c.moveTo(x,sy);c.lineTo(x+4,sy+Math.sign(spin)*7);c.stroke();}
  c.font='12px system-ui';c.textAlign='center';c.fillStyle=col;c.fillText(label,x,y+r+22);
 }
 const stage=t<-.012?'旋近':t<.008?'并合附近':'振铃';$('stage').textContent=stage;$('scene-time').textContent=`t = ${t.toFixed(4)} s`;$('scene-frequency').textContent=finite(wave.frequency[i])?`fGW ≈ ${wave.frequency[i].toFixed(1)} Hz`:t<-.012?'fGW = —（低幅度区间）':'旋近频率读数已结束';
}
function drawPolar(){const d=dimensions($('polar'));if(!d||!state.wave)return;const {c,w,h}=d,wave=state.wave,i=sampleIndex(state.time);let peak=0;for(const x of wave.envelope)peak=Math.max(peak,x);const mode=$('polar-mode').value,hp=mode==='cross'?0:wave.hp[i]/peak*.32,hc=mode==='plus'?0:wave.hc[i]/peak*.32,r=Math.min(w,h)*.32,cx=w/2,cy=h/2;
 c.strokeStyle=colors.grid;c.setLineDash([3,4]);c.beginPath();c.arc(cx,cy,r,0,TAU);c.stroke();c.setLineDash([]);
 c.strokeStyle=colors.grid;c.beginPath();c.moveTo(cx-r*1.5,cy);c.lineTo(cx+r*1.5,cy);c.moveTo(cx,cy-r*1.5);c.lineTo(cx,cy+r*1.5);c.stroke();
 for(let k=0;k<28;k++){const x=r*Math.cos(TAU*k/28),y=r*Math.sin(TAU*k/28),xx=x+.5*(hp*x+hc*y),yy=y+.5*(hc*x-hp*y);c.fillStyle=k%7===0?colors.gold:colors.cyan;c.beginPath();c.arc(cx+xx,cy+yy,k%7===0?3.6:2.4,0,TAU);c.fill();}
 c.fillStyle=colors.text;c.font='11px system-ui';c.textAlign='right';c.fillText('x',cx+r*1.5,cy-5);c.textAlign='left';c.fillText('y',cx+5,cy-r*1.4);
}
function chart(canvas,series,{xmin,xmax,ylabel='应变 / 10⁻²¹',cursor=null,shade=null}={}){
 const d=dimensions(canvas);if(!d)return;const {c,w,h}=d;const left=w<400?50:57,right=14,top=22,bottom=37,pw=w-left-right,ph=h-top-bottom;
 if(xmin===undefined){xmin=Math.min(...series.map(s=>s.t[0]));xmax=Math.max(...series.map(s=>s.t.at(-1)));}
 let max=0;for(const s of series)for(let i=0;i<s.y.length;i++)if(s.t[i]>=xmin&&s.t[i]<=xmax&&finite(s.y[i]))max=Math.max(max,Math.abs(s.y[i]));max=Math.max(max*1.15,1e-4);
 const X=x=>left+(x-xmin)/(xmax-xmin)*pw,Y=y=>top+ph/2-y/max*ph/2;
 c.font='11px system-ui';c.lineWidth=.7;
 for(let k=0;k<5;k++){const val=max*(1-k/2),y=Y(val);c.strokeStyle=colors.grid;c.beginPath();c.moveTo(left,y);c.lineTo(w-right,y);c.stroke();c.fillStyle=colors.text;c.textAlign='right';c.fillText(Math.abs(val)<1e-10?'0':val.toFixed(max<1?2:1),left-6,y+4);}
 const ticks=w<420?3:5;for(let k=0;k<ticks;k++){const val=xmin+(xmax-xmin)*k/(ticks-1),x=X(val);c.strokeStyle=colors.grid;c.beginPath();c.moveTo(x,top);c.lineTo(x,h-bottom);c.stroke();c.fillStyle=colors.text;c.textAlign=k===0?'left':k===ticks-1?'right':'center';c.fillText(val.toFixed(2),x,h-bottom+17);}
 c.fillStyle=colors.text;c.textAlign='left';c.fillText(ylabel,left,12);c.textAlign='center';c.fillText('时间 / s',left+pw/2,h-3);
 c.save();c.beginPath();c.rect(left,top,pw,ph);c.clip();
 if(shade){c.fillStyle='#62d9db09';c.fillRect(X(shade[0]),top,X(shade[1])-X(shade[0]),ph);}
 for(const s of series){c.strokeStyle=s.color;c.lineWidth=s.width||1.3;c.setLineDash(s.dashed?[4,3]:[]);c.beginPath();let active=false;for(let i=0;i<s.y.length;i++){if(s.t[i]<xmin-1/100||s.t[i]>xmax+1/100||!finite(s.y[i]))continue;const x=X(s.t[i]),y=Y(s.y[i]);active?c.lineTo(x,y):c.moveTo(x,y);active=true;}c.stroke();}c.setLineDash([]);
 if(cursor!==null&&cursor>=xmin&&cursor<=xmax){const x=X(cursor);c.strokeStyle='#e7f3ff';c.lineWidth=1;c.beginPath();c.moveTo(x,top);c.lineTo(x,h-bottom);c.stroke();}
 c.restore();canvas.plot={xmin,xmax,left,pw};
}
function drawWave(){if(!state.wave)return;const w=state.wave;const series=[];if($('show-plus').checked)series.push({t:w.t,y:w.hp,color:colors.cyan});if($('show-cross').checked)series.push({t:w.t,y:w.hc,color:colors.gold});if(state.reference)series.push({t:state.reference.t,y:state.reference.hp,color:colors.purple,dashed:true});let xmin=w.t[0],xmax=w.t.at(-1);if(state.reference){xmin=Math.min(xmin,state.reference.t[0]);xmax=Math.max(xmax,state.reference.t.at(-1));}if($('window').value!=='full'){const duration=+$('window').value;xmin=-duration*.85;xmax=duration*.15;}chart($('wave'),series,{xmin,xmax,cursor:state.time});}
function drawEvent(){if(!state.event)return;const det=$('detector').value,e=state.event.detectors[det],fit=state.fit;const series=[{t:e.t,y:e.strain,color:colors.gray,width:1}];if(fit)series.push({t:fit.t,y:fit.template,color:colors.cyan,width:1.6});chart($('event-wave'),series,{xmin:-.35,xmax:.15,ylabel:'白化振幅',shade:[-.2,.08]});chart($('residual'),fit?[{t:fit.t,y:fit.residual,color:colors.gold}]:[],{xmin:-.35,xmax:.15,ylabel:'白化残差'});drawSpectrogram(e.spectrogram);}
function drawSpectrogram(s){const d=dimensions($('spectrogram'));if(!d)return;const {c,w,h}=d,left=57,right=14,top=18,bottom=35,pw=w-left-right,ph=h-top-bottom;const flat=s.power.flat().sort((a,b)=>a-b),lo=flat[Math.floor(flat.length*.1)],hi=flat[Math.floor(flat.length*.995)],dt=s.t[1]-s.t[0],df=s.f[1]-s.f[0];
 c.save();c.beginPath();c.rect(left,top,pw,ph);c.clip();for(let j=0;j<s.f.length;j++)for(let k=0;k<s.t.length;k++){const v=Math.max(0,Math.min(1,(s.power[j][k]-lo)/(hi-lo)));c.fillStyle=`rgb(${Math.round(12+v*170)},${Math.round(24+v*190)},${Math.round(38+v*145)})`;const x=left+(s.t[k]-dt/2+.5)/.75*pw,y=top+(512-(s.f[j]+df/2))/480*ph;c.fillRect(x,y,dt/.75*pw+1,df/480*ph+1);}c.restore();
 c.font='11px system-ui';c.fillStyle=colors.text;c.textAlign='left';c.fillText('频率 / Hz',left,12);c.textAlign='right';for(const f of [64,256,480])c.fillText(f,left-7,top+(512-f)/480*ph+4);for(const t of [-.5,-.25,0,.25]){c.textAlign=t===-.5?'left':t===.25?'right':'center';c.fillText(t.toFixed(2),left+(t+.5)/.75*pw,h-18);}c.textAlign='center';c.fillText('时间 / s',left+pw/2,h-2);
}
function drawAll(){lensing?.setActive(state.tab==='simulation'&&$('scene-mode').value==='lensing');if(state.tab==='simulation'){drawScene();drawPolar();drawWave();if(state.wave)$('scrub').value=1000*(state.time-state.wave.t[0])/(state.wave.t.at(-1)-state.wave.t[0]);}else if(state.tab==='observation')drawEvent();}
function tick(now){if(!state.playing||!state.wave)return;const dt=last?Math.min((now-last)/1000,.15):0;last=now;state.time+=dt*+$('speed').value;if(state.time>=state.wave.t.at(-1)){state.time=state.wave.t.at(-1);pause();}drawAll();if(state.playing)raf=requestAnimationFrame(tick);}
function play(){if(!state.wave||state.busy)return;if(state.time>=state.wave.t.at(-1))state.time=state.wave.t[0];state.playing=true;last=0;$('play').textContent='Ⅱ 暂停';raf=requestAnimationFrame(tick);}
$('play').onclick=()=>state.playing?pause():play();$('restart').onclick=()=>{pause();if(state.wave)state.time=state.wave.t[0];drawAll();};
$('scrub').oninput=()=>{pause();if(state.wave)state.time=state.wave.t[0]+($('scrub').value/1000)*(state.wave.t.at(-1)-state.wave.t[0]);drawAll();};
$('speed').onchange=()=>{stopAudio();};$('window').onchange=drawWave;$('polar-mode').onchange=drawPolar;$('show-plus').onchange=drawWave;$('show-cross').onchange=drawWave;
$('wave').onpointerdown=e=>{if(!state.wave)return;const q=$('wave').plot,x=e.clientX-$('wave').getBoundingClientRect().left;pause();state.time=Math.max(state.wave.t[0],Math.min(state.wave.t.at(-1),q.xmin+(x-q.left)/q.pw*(q.xmax-q.xmin)));drawAll();};
$('preset').onchange=()=>{state.p={...presets[$('preset').value],model:state.p.model};syncInputs();clearTimeout(debounce);updateWave();};$('model').onchange=()=>{state.p.model=$('model').value;clearTimeout(debounce);updateWave();};
$('save-reference').onclick=()=>{if(!state.wave||state.busy)return;state.reference=structuredClone(state.wave);$('ref-legend').hidden=false;$('clear-reference').disabled=false;const p=state.reference.parameters;$('reference-note').textContent=`已保存：${p.m1} + ${p.m2} M☉，χ₁=${p.s1}，χ₂=${p.s2}，${p.distance} Mpc，${p.model}。紫色虚线为其 h+。`;drawWave();};
$('clear-reference').onclick=()=>{state.reference=null;$('ref-legend').hidden=true;$('clear-reference').disabled=true;$('reference-note').textContent='质量改变时间尺度，自旋改变相位演化；距离和倾角影响观测振幅。';drawWave();};
$('export').onclick=()=>{if(!state.wave||state.busy)return;const out={...state.wave,reference_parameters:state.reference?.parameters??null,event_fit:state.fit,notes:['Detector-frame masses; hp/hc in 1e-21','Trajectory is schematic','Time aligned at face-on envelope peak; input phase_c=0','20–25 and 900–1000 Hz frequency tapers']};const blob=new Blob([JSON.stringify(out,null,2)],{type:'application/json'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download='bbh-waveform-'+state.wave.parameters.model+'.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);};
$('audio').onclick=async()=>{if(!state.wave||state.busy)return;if(audioSource){pause();drawAll();return;}pause();audioContext??=new AudioContext();await audioContext.resume();const w=state.wave,buf=audioContext.createBuffer(1,w.hp.length,w.sample_rate),dest=buf.getChannelData(0);let peak=0;for(const x of w.hp)peak=Math.max(peak,Math.abs(x));for(let i=0;i<dest.length;i++)dest[i]=w.hp[i]/Math.max(peak,1e-30)*.2*Math.min(1,i/100,(dest.length-1-i)/100);audioSource=audioContext.createBufferSource();audioSource.buffer=buf;audioSource.connect(audioContext.destination);audioSource.start();$('audio').textContent='■ 停止声音';$('speed').value='1';state.time=w.t[0];play();};
let drag=null;$('orbit').onpointerdown=e=>{drag={x:e.clientX,y:e.clientY,az:state.camera.az,el:state.camera.el};$('orbit').setPointerCapture(e.pointerId);};$('orbit').onpointermove=e=>{if(!drag)return;state.camera.az=drag.az+(e.clientX-drag.x)*.007;state.camera.el=Math.max(.12,Math.min(Math.PI/2,drag.el+(e.clientY-drag.y)*.007));drawOrbit();};$('orbit').onpointerup=$('orbit').onpointercancel=()=>drag=null;
$('camera-top').onclick=()=>{state.camera.el=Math.PI/2;lensing?.cameraView(true);drawScene();};$('camera-oblique').onclick=()=>{state.camera={az:.3,el:.6};lensing?.cameraView(false);drawScene();};
document.querySelectorAll('[data-tab]').forEach(b=>b.onclick=()=>{pause();state.tab=b.dataset.tab;document.querySelectorAll('[data-tab]').forEach(x=>{x.classList.toggle('active',x===b);x.setAttribute('aria-pressed',x===b?'true':'false');});$('workspace').hidden=state.tab==='methods';$('methods').hidden=state.tab!=='methods';$('simulation').hidden=state.tab!=='simulation';$('observation').hidden=state.tab!=='observation';if(state.tab==='observation')loadEvent();drawAll();});
async function loadEvent(){if(state.event){drawEvent();return;}try{const r=await fetch('/api/event');const data=await r.json();if(!r.ok)throw Error(data.error);state.event=data;drawEvent();}catch(e){$('fit-summary').textContent='真实数据未加载：'+e.message;}}
$('detector').onchange=()=>{fitId++;state.fit=null;$('fit-summary').textContent='探测器已切换，点击对齐当前模板。';drawEvent();};
$('fit').onclick=async()=>{if(!state.wave||state.busy)return;const id=++fitId;const detector=$('detector').value;$('fit').disabled=true;$('fit-summary').textContent='正在进行独立时间、相位与振幅对齐…';try{const fit=await post('/api/fit',{parameters:state.wave.parameters,detector});if(id!==fitId)return;state.fit=fit;$('fit-summary').textContent=`${detector} · 窗口相关系数 ${fit.correlation.toFixed(3)} · 模板平移 ${fit.time_shift_ms.toFixed(1)} ms · 拟合窗口 −0.20 至 +0.08 s`;drawEvent();}catch(e){if(id===fitId)$('fit-summary').textContent=e.message;}finally{$('fit').disabled=false;}};
const resize=new ResizeObserver(()=>drawAll());document.querySelectorAll('canvas').forEach(c=>resize.observe(c));
document.addEventListener('visibilitychange',()=>{if(document.hidden)pause();});
window.bbhLab={state,updateWave,drawAll};syncInputs();drawAll();updateWave();initLensing();
})();
