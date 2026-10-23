const finite=x=>typeof x==='number'&&Number.isFinite(x);
export function validateDynamics(d){
 if(d?.schema!=='blackhole-dynamics/1'||!Array.isArray(d.t)||d.t.length<2||d.t.length>50000||d.t.some((t,i)=>!finite(t)||(i&&t<=d.t[i-1])))throw Error('代理时间数据无效');
 for(const [name,n] of [['q_copr',4],['chiA',3],['chiB',3]])if(!Array.isArray(d[name])||d[name].length!==d.t.length||d[name].some(v=>!Array.isArray(v)||v.length!==n||v.some(x=>!finite(x))))throw Error('代理向量数据无效：'+name);
 if(d.q_copr.some(q=>Math.abs(Math.hypot(...q)-1)>1e-5))throw Error('旋转四元数未归一化');
 if([...d.chiA,...d.chiB].some(s=>Math.hypot(...s)>1+1e-8))throw Error('无量纲自旋超过有效范围');
 for(const name of ['orbphase','h22_real','h22_imag'])if(!Array.isArray(d[name])||d[name].length!==d.t.length||d[name].some(v=>!finite(v)))throw Error('代理标量数据无效：'+name);
 if(!finite(d.dynamics_valid_until_M)||d.dynamics_valid_until_M<d.t[0]||d.dynamics_valid_until_M>d.t.at(-1))throw Error('缺少有效动力学显示区间');
 return d;
}
// Scalar-first unit quaternion, active rotation from coprecessing to inertial.
export function rotateQuaternion(q,v){const [w,x,y,z]=q,[a,b,c]=v;const uv=[y*c-z*b,z*a-x*c,x*b-y*a],uuv=[y*uv[2]-z*uv[1],z*uv[0]-x*uv[2],x*uv[1]-y*uv[0]];return v.map((r,i)=>r+2*(w*uv[i]+uuv[i]));}
export function nearestDynamicsFrame(d,t){let lo=0,hi=d.t.length-1;while(lo<hi){const m=Math.floor((lo+hi)/2);if(d.t[m]<t)lo=m+1;else hi=m}return lo&&Math.abs(d.t[lo-1]-t)<Math.abs(d.t[lo]-t)?lo-1:lo;}
export function createDynamicsView(){const $=id=>document.getElementById(id);let data=null,index=0,playing=false,raf=null,previous=0,time=0,hash=null;
 function canvas(id){const e=$(id),r=e.getBoundingClientRect();if(r.width<1)return null;const k=Math.min(devicePixelRatio||1,2);e.width=Math.round(r.width*k);e.height=Math.round(r.height*k);const c=e.getContext('2d');c.setTransform(k,0,0,k,0,0);return {c,w:r.width,h:r.height};}
 function orientation(){const d=canvas('dynamics-orientation');if(!d)return;const {c,w,h}=d,scale=Math.min(w,h)*.32,az=.65,el=.45;const p=v=>[w/2+(v[0]*Math.cos(az)-v[1]*Math.sin(az))*scale,h/2+((v[0]*Math.sin(az)+v[1]*Math.cos(az))*Math.sin(el)-v[2]*Math.cos(el))*scale];
 c.strokeStyle='#2d4054';c.beginPath();c.arc(w/2,h/2,scale,0,2*Math.PI);c.stroke();c.font='12px system-ui';c.textAlign='center';
 if(data.t[index]>data.dynamics_valid_until_M+1e-7){c.fillStyle='#94a7b9';c.fillText('已超出动力学解释区间',w/2,h/2);$('dynamics-frame-note').textContent='保留完整波形；不将晚期延伸自旋展示为物理黑洞自旋。';return;}
 const q=data.q_copr[index],phi=data.orbphase[index];
 c.beginPath();for(let j=0;j<=80;j++){const a=j/80*2*Math.PI,pt=p(rotateQuaternion(q,[Math.cos(a),Math.sin(a),0]));j?c.lineTo(...pt):c.moveTo(...pt)}c.stroke();
 for(const [name,v,col,dashed] of [['χA',data.chiA[index],'#62d9db',false],['χB',data.chiB[index],'#eeb977',false],['z copr',rotateQuaternion(q,[0,0,1]),'#b7a2ef',false],['x coorb',rotateQuaternion(q,[Math.cos(phi),Math.sin(phi),0]),'#94a7b9',true]]){const pt=p(v);c.strokeStyle=col;c.fillStyle=col;c.lineWidth=2;c.setLineDash(dashed?[4,4]:[]);c.beginPath();c.moveTo(w/2,h/2);c.lineTo(...pt);c.stroke();c.setLineDash([]);c.beginPath();c.arc(...pt,3,0,2*Math.PI);c.fill();c.fillText(name,pt[0],pt[1]-9)}
 $('dynamics-frame-note').textContent='球半径为 1；χ 箭头按无量纲模长显示，框架方向为单位向量。圆环仅标记框架平面。';
 }
 function wave(){const d=canvas('dynamics-wave');if(!d)return;const {c,w,h}=d,L=52,R=14,T=25,B=40,pw=w-L-R,ph=h-T-B,t=data.t;let max=0;for(const a of [data.h22_real,data.h22_imag])for(const v of a)max=Math.max(max,Math.abs(v));max=max*1.1||1;const X=x=>L+(x-t[0])/(t.at(-1)-t[0])*pw,Y=y=>T+ph/2-y/max*ph/2;
 c.fillStyle='#94a7b9';c.font='10px system-ui';for(let j=0;j<=4;j++){const y=T+ph*j/4;c.strokeStyle='#253345';c.beginPath();c.moveTo(L,y);c.lineTo(w-R,y);c.stroke();c.textAlign='right';c.fillText((max*(1-j/2)).toFixed(2),L-6,y+4);c.textAlign='center';c.fillText((t[0]+(t.at(-1)-t[0])*j/4).toFixed(0),L+pw*j/4,h-20)}
 [[data.h22_real,'#62d9db'],[data.h22_imag,'#eeb977']].forEach(([a,col])=>{c.strokeStyle=col;c.lineWidth=1;c.beginPath();a.forEach((v,i)=>i?c.lineTo(X(t[i]),Y(v)):c.moveTo(X(t[i]),Y(v)));c.stroke()});
 c.strokeStyle='#ffffff';c.beginPath();c.moveTo(X(t[index]),T);c.lineTo(X(t[index]),T+ph);c.stroke();c.setLineDash([3,4]);c.strokeStyle='#b7a2ef';c.beginPath();c.moveTo(X(data.dynamics_valid_until_M),T);c.lineTo(X(data.dynamics_valid_until_M),T+ph);c.stroke();c.setLineDash([]);c.fillStyle='#94a7b9';c.textAlign='center';c.fillText('t / M · 虚线为动力学显示截止',L+pw/2,h-3);
 }
 function render(){if(!data||$('dynamics-view').hidden)return;index=nearestDynamicsFrame(data,time);orientation();wave();$('dynamics-time').value=(time-data.t[0])/(data.t.at(-1)-data.t[0])*1000;$('dynamics-clock').textContent=`最近采样帧 t/M = ${data.t[index].toFixed(3)} · ${index+1} / ${data.t.length}`;$('dynamics-values').textContent=data.t[index]<=data.dynamics_valid_until_M+1e-7?`|χA|=${Math.hypot(...data.chiA[index]).toFixed(6)} · |χB|=${Math.hypot(...data.chiB[index]).toFixed(6)} · 框架轨道相位=${data.orbphase[index].toFixed(4)} rad`:'当前仅解释波形模式。';}
 function stop(){playing=false;cancelAnimationFrame(raf);raf=null;$('dynamics-play').textContent='▶ 播放';}
 function tick(now){if(!playing||$('dynamics-view').hidden){stop();return}const dt=previous?Math.min(.1,(now-previous)/1000):0;previous=now;time+=dt*(data.t.at(-1)-data.t[0])/30;if(time>=data.t.at(-1)){time=data.t.at(-1);stop()}render();if(playing)raf=requestAnimationFrame(tick)}
 $('dynamics-play').onclick=()=>{if(playing){stop();return}if(!data)return;if(time>=data.t.at(-1))time=data.t[0];previous=0;playing=true;$('dynamics-play').textContent='Ⅱ 暂停';raf=requestAnimationFrame(tick)};
 $('dynamics-time').oninput=()=>{if(!data)return;stop();time=data.t[0]+(data.t.at(-1)-data.t[0])*(+$('dynamics-time').value)/1000;render()};$('dynamics-cutoff').onclick=()=>{stop();time=data.dynamics_valid_until_M;render()};
 $('dynamics-export').onclick=()=>{const url=URL.createObjectURL(new Blob([JSON.stringify({data_sha256:hash,output:data},null,2)],{type:'application/json'}));const a=document.createElement('a');a.href=url;a.download='nrsur7dq4-output.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000)};
 async function load(){try{const r=await fetch('datasets/dynamics.json');if(!r.ok)throw Error('代理输出尚未就绪');const raw=await r.arrayBuffer();hash=Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',raw)),x=>x.toString(16).padStart(2,'0')).join('');data=validateDynamics(JSON.parse(new TextDecoder().decode(raw)));time=data.t[0];$('dynamics-source').textContent=JSON.stringify({metadata:data.metadata,dynamics_valid_until_M:data.dynamics_valid_until_M,data_sha256:hash},null,2);$('dynamics-status').textContent=`已载入实际 NRSur7dq4 模型输出。t/M > ${data.dynamics_valid_until_M} 的延伸自旋没有物理黑洞自旋含义，方向球在该范围停止显示；波形继续保留。`;for(const id of ['dynamics-play','dynamics-time','dynamics-cutoff','dynamics-export'])$(id).disabled=false;render()}catch(e){$('dynamics-status').textContent=e.message}}
 document.addEventListener('visibilitychange',()=>{if(document.hidden)stop()});load().then(()=>{if(data){const p=data.metadata.params;$('dynamics-params').textContent=`本次固定输入：q=${p.q}，χA=(${p.chiA0.join(', ')})，χB=(${p.chiB0.join(', ')})。自旋参考为模型起始时刻；未指定总质量、距离或探测器。`}});return {render,stop,get data(){return data}};
}
