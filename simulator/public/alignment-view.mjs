const METRICS={relative_position_difference_M:{label:'相对位置向量差的范数 / Mref',rms:'rms_relative_position_M'},separation_difference_M:{label:'分离度绝对差 / Mref',rms:'rms_separation_difference_M'},planar_phase_residual_rad:{label:'平面相位差 / rad',rms:'rms_planar_phase_residual_rad'}};
const COLORS=['#62d9db','#eeb977','#b7a2ef'];
const finite=x=>typeof x==='number'&&Number.isFinite(x);
export function validateAlignment(data){
 if(data?.schema!=='blackhole-alignment-study/1'||!Array.isArray(data.time_M)||data.time_M.length<2||data.time_M.length>30000)throw Error('比较数据格式无效');
 if(data.time_M.some((t,i)=>!finite(t)||(i&&t<=data.time_M[i-1])))throw Error('比较时间必须严格递增');
 if(!Array.isArray(data.methods)||data.methods.length<1||data.methods.length>10)throw Error('缺少比较方法');
 for(const m of data.methods){if(typeof m.label!=='string'||!finite(m.delta_t_M)||!finite(m.phase_rotation_rad))throw Error('对齐参数无效');for(const [key,meta] of Object.entries(METRICS)){if(!Array.isArray(m.series?.[key])||m.series[key].length!==data.time_M.length||m.series[key].some(v=>!finite(v))||!finite(m.validation?.[meta.rms]))throw Error('残差曲线或指标无效')}}
 for(const window of [data.training_window_M,data.validation_window_M])if(!Array.isArray(window)||window.length!==2||!window.every(finite)||window[1]<=window[0])throw Error('训练或验证窗口无效');
 if(data.training_window_M[1]>=data.validation_window_M[0])throw Error('训练和验证窗口重叠');
 return data;
}
export function createAlignmentView(){
 const $=id=>document.getElementById(id);let data=null,index=0,contentHash=null;
 const format=x=>Math.abs(x)<.001&&x!==0?x.toExponential(3):x.toFixed(6);
 function row(parent,values){const tr=document.createElement('tr');for(const value of values){const td=document.createElement('td');td.textContent=value;tr.append(td)}parent.append(tr)}
 function table(){const key=$('alignment-metric').value,meta=METRICS[key];const body=$('alignment-table').tBodies[0];body.replaceChildren();$('alignment-table').tHead.rows[0].cells[3].textContent=meta.label+' · RMS';for(const m of data.methods)row(body,[m.label,format(m.delta_t_M),format(m.phase_rotation_rad),format(m.validation[meta.rms])]);const sensitivity=$('alignment-sensitivity').tBodies[0];sensitivity.replaceChildren();$('alignment-sensitivity').tHead.rows[0].cells[2].textContent=meta.label+' · 验证 RMS';for(const m of data.training_window_sensitivity)row(sensitivity,[`${Math.round(m.training_fraction*100)}% · ${m.training_window_M.map(x=>x.toFixed(0)).join('–')} M`,format(m.delta_t_M),format(m.validation[meta.rms])]);}
 function render(){
  if(!data||$('alignment-view').hidden)return;const e=$('alignment-chart'),r=e.getBoundingClientRect();if(r.width<1||r.height<1)return;const dpr=Math.min(devicePixelRatio||1,2);e.width=Math.round(r.width*dpr);e.height=Math.round(r.height*dpr);const c=e.getContext('2d');c.setTransform(dpr,0,0,dpr,0,0);const w=r.width,h=r.height,L=64,R=18,T=30,B=40,pw=w-L-R,ph=h-T-B,key=$('alignment-metric').value,meta=METRICS[key];
  let lo=0,hi=0;for(const m of data.methods)for(const value of m.series[key]){lo=Math.min(lo,value);hi=Math.max(hi,value)}const pad=(hi-lo||1)*.06;lo-=lo<0?pad:0;hi+=pad;const t=data.time_M,X=x=>L+(x-t[0])/(t.at(-1)-t[0])*pw,Y=y=>T+(hi-y)/(hi-lo)*ph;
  c.font='11px system-ui';c.fillStyle='#94a7b9';c.textAlign='left';c.fillText(meta.label,L,15);for(let j=0;j<=4;j++){const y=T+ph*j/4,v=hi-(hi-lo)*j/4;c.strokeStyle='#253345';c.beginPath();c.moveTo(L,y);c.lineTo(w-R,y);c.stroke();c.fillStyle='#94a7b9';c.textAlign='right';c.fillText(Math.abs(v)<.001&&v!==0?v.toExponential(1):v.toFixed(3),L-7,y+4);c.textAlign='center';c.fillText((t[0]+(t.at(-1)-t[0])*j/4).toFixed(0),L+pw*j/4,h-20)}
  data.methods.forEach((m,j)=>{c.strokeStyle=COLORS[j%COLORS.length];c.lineWidth=1.6;c.beginPath();m.series[key].forEach((value,i)=>i?c.lineTo(X(t[i]),Y(value)):c.moveTo(X(t[i]),Y(value)));c.stroke()});
  c.strokeStyle='#ffffff';c.setLineDash([3,4]);c.beginPath();c.moveTo(X(t[index]),T);c.lineTo(X(t[index]),T+ph);c.stroke();c.setLineDash([]);c.fillStyle='#94a7b9';c.textAlign='center';c.fillText('Lev3 原始坐标时间 / Mref',L+pw/2,h-3);
  $('alignment-cursor').textContent=`t/Mref = ${t[index].toFixed(2)} · `+data.methods.map(m=>`${m.label}: ${format(m.series[key][index])}`).join(' · ');table();
 }
 $('alignment-metric').onchange=render;
 $('alignment-chart').onpointermove=e=>{if(!data)return;const r=e.target.getBoundingClientRect(),u=Math.max(0,Math.min(1,(e.clientX-r.left-64)/(r.width-82)));index=Math.round(u*(data.time_M.length-1));render()};
 $('alignment-export').onclick=()=>{if(!data)return;const url=URL.createObjectURL(new Blob([JSON.stringify({schema:'alignment-session/1',data_sha256:contentHash,view:{metric:$('alignment-metric').value,cursor_index:index},study:data},null,2)],{type:'application/json'}));const a=document.createElement('a');a.href=url;a.download='alignment-study.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000)};
 async function load(){try{const r=await fetch('datasets/alignment.json');if(!r.ok)throw Error('对齐研究结果读取失败');const raw=await r.arrayBuffer();contentHash=Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',raw)),b=>b.toString(16).padStart(2,'0')).join('');data=validateAlignment(JSON.parse(new TextDecoder().decode(raw)));index=data.time_M.length-1;
  $('alignment-legend').replaceChildren();data.methods.forEach((m,j)=>{const span=document.createElement('span'),i=document.createElement('i');i.className='swatch';i.style.background=COLORS[j%3];span.append(i,document.createTextNode(m.label));$('alignment-legend').append(span)});
  $('alignment-windows').textContent=`早期训练窗口：${data.training_window_M.map(x=>x.toFixed(2)).join('–')} Mref；独立晚期验证窗口：${data.validation_window_M.map(x=>x.toFixed(2)).join('–')} Mref。图中显示验证区间。`;
  $('alignment-transform').textContent='比较 Rz(φ) rLev2(t+Δt) 与 rLev3(t)。只使用一个常量时间偏移和一个绕原始 +z 轴的旋转；共同视界规则使用晚期标志，属于回溯比较。';
  $('alignment-status').textContent='已载入原始视界数据计算结果；合成偏移恢复、无外推和窗口隔离检查已保存。';const {time_M,methods,...source}=data;$('alignment-source').textContent=JSON.stringify({...source,data_sha256:contentHash},null,2);$('alignment-export').disabled=false;render();
 }catch(error){$('alignment-status').textContent=error.message}}
 load();return {render,get data(){return data},get contentHash(){return contentHash}};
}

