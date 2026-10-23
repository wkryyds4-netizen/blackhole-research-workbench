const finite=x=>typeof x==='number'&&Number.isFinite(x);
const object=x=>x!==null&&typeof x==='object';
const text=x=>typeof x==='string'&&x.length>0&&x.length<1000;
const COLORS={surrogate:'#62d9db',extrapolation:'#eeb977'};
const LABELS={surrogate:'NRSur7dq4 对 NR',extrapolation:'NR 外推阶之间'};
const format=x=>x===0?'0':Math.abs(x)<.001?x.toExponential(4):x.toFixed(6);

export function validateWaveformStudy(d){
 if(d?.schema!=='blackhole-waveform-study/1'||!Array.isArray(d.rows)||d.rows.length<1||d.rows.length>5000)throw Error('波形研究格式或数值行无效');
 const ids=new Set();
 for(const r of d.rows){
  if(!text(r.id)||ids.has(r.id)||!text(r.label)||!['surrogate','extrapolation','sampling'].includes(r.comparison)||!text(r.psd))throw Error('波形比较行标识无效或重复');ids.add(r.id);
  for(const key of ['mass_msun','fmin_hz','fmax_hz','taper_fraction','dt_M','mismatch','match','time_shift_M','phase_shift_rad','start_frequency_hz'])if(!finite(r[key]))throw Error('波形比较数值无效：'+key);
  if(r.mass_msun<=0||r.dt_M<=0||r.fmin_hz<=0||r.fmax_hz<=r.fmin_hz||r.start_frequency_hz<0||r.taper_fraction<0||r.taper_fraction>1||r.mismatch<0||r.mismatch>1||r.match<0||r.match>1||Math.abs(r.match+r.mismatch-1)>1e-7)throw Error('波形比较数值超出范围');
  if(typeof r.at_time_boundary!=='boolean'||typeof r.start_band_supported!=='boolean')throw Error('缺少时间边界或低频覆盖标记');
  if(r.peak_window_weight!==undefined&&(!finite(r.peak_window_weight)||r.peak_window_weight<0||r.peak_window_weight>1))throw Error('并合峰窗权重无效');
 }
 if(!text(d.baseline_id)||!ids.has(d.baseline_id))throw Error('波形参考配置缺失');
 if(!d.rows.some(r=>Math.abs(r.dt_M-1)<1e-10))throw Error('缺少 dt/M=1 主扫描');
 const w=d.waveform,t=w?.time_M;
 if(!Array.isArray(t)||t.length<2||t.length>50000||t.some((x,i)=>!finite(x)||(i>0&&x<=t[i-1])))throw Error('波形显示时间必须有限且严格递增');
 for(const key of ['nr_re','surrogate_re','residual_abs'])if(!Array.isArray(w[key])||w[key].length!==t.length||w[key].some(x=>!finite(x)||(key==='residual_abs'&&x<0)))throw Error('波形显示数组无效：'+key);
 for(const key of ['provenance','protocol','limits'])if(!object(d[key]))throw Error('缺少研究来源或边界：'+key);
 return d;
}

export function selectWaveformRows(d,{psd,fmin,taper}){
 return d.rows.filter(r=>r.comparison!=='sampling'&&Math.abs(r.dt_M-1)<1e-10&&r.psd===psd&&Math.abs(r.fmin_hz-Number(fmin))<1e-9&&Math.abs(r.taper_fraction-Number(taper))<1e-10);
}

export function groupWaveformRows(rows){
 const groups=new Map();for(const r of rows){const key=JSON.stringify([r.comparison,r.label,r.fmax_hz]);if(!groups.has(key))groups.set(key,[]);groups.get(key).push(r)}return [...groups.values()].map(group=>group.sort((a,b)=>a.mass_msun-b.mass_msun));
}

export function createWaveformView({url='datasets/waveform-comparison.json'}={}){
 const $=id=>document.getElementById(id);let data=null,contentHash=null;
 const required=['waveform-view','waveform-status','waveform-psd','waveform-fmin','waveform-taper','waveform-export','waveform-summary','waveform-chart','waveform-legend','waveform-wave','waveform-baseline','waveform-table','waveform-source'];
 for(const id of required)if(!$(id))throw Error('波形面板缺少元素：'+id);
 const filters=()=>({psd:$('waveform-psd').value,fmin:+$('waveform-fmin').value,taper:+$('waveform-taper').value});
 function colorOf(r){if(r.comparison==='surrogate')return COLORS.surrogate;const labels=[...new Set(data.rows.filter(x=>x.comparison==='extrapolation').map(x=>x.label))].sort();return ['#eeb977','#b7a2ef','#e08fb8'][labels.indexOf(r.label)%3]||COLORS.extrapolation;}
 function canvas(id){const e=$(id),r=e.getBoundingClientRect();if(r.width<1||r.height<1)return null;const k=Math.min(globalThis.devicePixelRatio||1,2);e.width=Math.round(r.width*k);e.height=Math.round(r.height*k);const c=e.getContext('2d');c.setTransform(k,0,0,k,0,0);return {c,w:r.width,h:r.height};}
 function fillSelect(id,values,selected,label){const e=$(id);e.replaceChildren();for(const v of values){const o=document.createElement('option');o.value=String(v);o.textContent=label(v);e.append(o)}e.value=String(selected);if(e.selectedIndex<0)e.selectedIndex=0;e.disabled=false;}
 function line(c,points,color,dashed=false){c.strokeStyle=color;c.lineWidth=1.6;c.setLineDash(dashed?[6,4]:[]);c.beginPath();points.forEach(([x,y],i)=>i?c.lineTo(x,y):c.moveTo(x,y));c.stroke();c.setLineDash([]);}
 function chart(rows){const v=canvas('waveform-chart');if(!v)return;const {c,w,h}=v,L=68,R=20,T=29,B=50,pw=w-L-R,ph=h-T-B;c.font='11px system-ui';c.fillStyle='#94a7b9';
  if(!rows.length){c.fillText('此筛选组合没有主扫描结果。',L,T);return;}
  let loM=Math.min(...rows.map(r=>r.mass_msun)),hiM=Math.max(...rows.map(r=>r.mass_msun));if(loM===hiM){loM-=5;hiM+=5}
  const plotted=rows.map(r=>Math.max(1e-12,r.mismatch));let lo=Math.floor(Math.log10(Math.min(...plotted))),hi=Math.ceil(Math.log10(Math.max(...plotted)));if(hi===lo)hi=lo+1;
  const X=x=>L+(x-loM)/(hiM-loM)*pw,Y=y=>T+(hi-Math.log10(Math.max(1e-12,y)))/(hi-lo)*ph;
  c.textAlign='left';c.fillText('1 − match · 对数坐标',L,14);
  const stride=Math.max(1,Math.ceil((hi-lo)/6));
  for(let exponent=lo;exponent<=hi;exponent+=stride){const y=Y(10**exponent);c.strokeStyle='#253345';c.beginPath();c.moveTo(L,y);c.lineTo(w-R,y);c.stroke();c.textAlign='right';c.fillText(`1e${exponent}`,L-8,y+4)}
  [...new Set(rows.map(r=>r.mass_msun))].sort((a,b)=>a-b).forEach(m=>{c.textAlign='center';c.fillText(String(m),X(m),h-29)});
  for(const group of groupWaveformRows(rows)){const color=colorOf(group[0]);line(c,group.map(r=>[X(r.mass_msun),Y(r.mismatch)]),color,group[0].comparison==='extrapolation');for(const r of group){const x=X(r.mass_msun),y=Y(r.mismatch);c.strokeStyle=color;c.fillStyle=color;c.beginPath();c.arc(x,y,r.start_band_supported?3:5,0,2*Math.PI);r.start_band_supported?c.fill():c.stroke();if(r.at_time_boundary){c.fillStyle='#f18a8a';c.fillRect(x-3,y-3,6,6)}}}
  c.fillStyle='#94a7b9';c.textAlign='center';c.fillText('探测器系总质量 / M☉',L+pw/2,h-8);
 }
 function wave(){const v=canvas('waveform-wave');if(!v)return;const {c,w,h}=v,a=data.waveform,t=a.time_M,L=62,R=17,T=29,B=36,pw=w-L-R;let lo=Math.max(t[0],-200),hi=Math.min(t.at(-1),100);if(hi<=lo){lo=t[0];hi=t.at(-1)}const indices=t.map((x,i)=>x>=lo&&x<=hi?i:-1).filter(i=>i>=0);if(indices.length<2)return;const upperH=(h-T-B)*.59,lowerT=T+upperH+39,lowerH=h-B-lowerT;
  const X=x=>L+(x-lo)/(hi-lo)*pw;let amplitude=0,residual=0;for(const i of indices){amplitude=Math.max(amplitude,Math.abs(a.nr_re[i]),Math.abs(a.surrogate_re[i]));residual=Math.max(residual,a.residual_abs[i])}amplitude=amplitude*1.08||1;residual=residual*1.08||1;const Y=y=>T+upperH/2-y/amplitude*upperH/2,Z=y=>lowerT+lowerH-y/residual*lowerH;
  c.font='10px system-ui';c.fillStyle='#94a7b9';c.textAlign='left';c.fillText('Re(h₂₂) · r h/M · 固定参考配置',L,14);
  for(const value of [-amplitude,0,amplitude]){c.strokeStyle='#253345';c.beginPath();c.moveTo(L,Y(value));c.lineTo(w-R,Y(value));c.stroke();c.textAlign='right';c.fillText(value.toFixed(3),L-7,Y(value)+3)}
  line(c,indices.map(i=>[X(t[i]),Y(a.nr_re[i])]),'#62d9db');line(c,indices.map(i=>[X(t[i]),Y(a.surrogate_re[i])]),'#eeb977',true);
  c.textAlign='left';c.fillStyle='#94a7b9';c.fillText('|复模式残差|（独立纵轴）',L,lowerT-10);line(c,indices.map(i=>[X(t[i]),Z(a.residual_abs[i])]),'#b7a2ef');c.textAlign='right';c.fillText(residual.toExponential(1),L-7,lowerT+4);c.fillText('0',L-7,lowerT+lowerH+3);
  for(let j=0;j<=4;j++){const x=lo+(hi-lo)*j/4;c.textAlign='center';c.fillText(x.toFixed(0),X(x),h-20)}c.fillText('NR 辐射波形参考时间 / M',L+pw/2,h-4);
 }
 function table(){const e=$('waveform-table');e.replaceChildren();const head=e.createTHead(),tr=head.insertRow();for(const title of ['ID / 比较','M☉','PSD','频段 / Hz','Tukey α','dt/M','失配','match','Δt/M','相位 / rad','起始频率 / Hz','低频覆盖','搜索边界']){const th=document.createElement('th');th.scope='col';th.textContent=title;tr.append(th)}const body=e.createTBody();for(const r of data.rows){const tr=body.insertRow();for(const value of [`${r.id} · ${r.label}`,r.mass_msun,r.psd,`${r.fmin_hz}–${r.fmax_hz}`,r.taper_fraction,r.dt_M,format(r.mismatch),format(r.match),format(r.time_shift_M),format(r.phase_shift_rad),r.start_frequency_hz.toFixed(3),r.start_band_supported?'起始频率 ≤ 下限':'起始频率高于下限',r.at_time_boundary?'触及边界':'否']){const td=tr.insertCell();td.textContent=String(value)}}}
 function render(){if(!data||$('waveform-view').hidden)return;const rows=selectWaveformRows(data,filters());chart(rows);wave();const gaps=rows.filter(r=>!r.start_band_supported).length,bounds=rows.filter(r=>r.at_time_boundary).length,weights=[...new Set(rows.map(r=>r.peak_window_weight).filter(finite).map(x=>x.toFixed(3)))];$('waveform-summary').textContent=`主扫描 dt/M=1：当前 ${rows.length} 行；${gaps} 行记录起始局部频率高于频带下限，${bounds} 行达到时间搜索边界。${weights.length?`当前窗在 t=0 处的振幅权重：${weights.join(' / ')}（1 表示不压低，0 表示完全去除）。`:''}空心圆标记低频未覆盖，红方块标记搜索边界。低频起始条件通过也不证明波形长度充分；需结合窗口、网格和频段敏感性。零失配在图上画于 1e−12，表格保留精确数值。`;}
 for(const id of ['waveform-psd','waveform-fmin','waveform-taper'])$(id).onchange=render;
 $('waveform-export').onclick=()=>{if(!data)return;const objectUrl=URL.createObjectURL(new Blob([JSON.stringify({schema:'blackhole-waveform-session/1',data_sha256:contentHash,view:filters(),study:data},null,2)],{type:'application/json'}));const a=document.createElement('a');a.href=objectUrl;a.download='waveform-study.json';a.click();setTimeout(()=>URL.revokeObjectURL(objectUrl),1000)};
 async function load(){try{const r=await fetch(url);if(!r.ok)throw Error('波形比较结果尚未就绪');const bytes=await r.arrayBuffer();if(bytes.byteLength>20*1024*1024)throw Error('波形研究文件过大');contentHash=Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',bytes)),v=>v.toString(16).padStart(2,'0')).join('');data=validateWaveformStudy(JSON.parse(new TextDecoder().decode(bytes)));const baseline=data.rows.find(r=>r.id===data.baseline_id),main=data.rows.filter(r=>Math.abs(r.dt_M-1)<1e-10),unique=key=>[...new Set(main.map(r=>r[key]))];
  fillSelect('waveform-psd',unique('psd'),baseline.psd,v=>v==='white'?'白 PSD':v==='aligo_historical'?'LAL 历史解析 PSD':v);fillSelect('waveform-fmin',unique('fmin_hz').sort((a,b)=>a-b),baseline.fmin_hz,v=>`${v} Hz`);fillSelect('waveform-taper',unique('taper_fraction').sort((a,b)=>a-b),baseline.taper_fraction,v=>`α=${v}`);
  $('waveform-status').textContent='比较对象是数学复模式 h₂₂，并非探测器应变。所选 NR 是否属于代理训练集尚未核实，不能称作独立测试集验证。外推阶差异是一项处理敏感性，不是误差条。长窗可能压低并合峰，应同时检查 Tukey α 的变化。';
  $('waveform-baseline').textContent=`波形图固定参考：${baseline.label}；${baseline.mass_msun} M☉，${baseline.psd}，${baseline.fmin_hz}–${baseline.fmax_hz} Hz，α=${baseline.taper_fraction}。显示未加窗的 NR（青）与全局时间/相位对齐代理（橙虚线），未作振幅缩放；紫色为复模式残差幅度。筛选仅改变上方失配图。`;
  $('waveform-legend').replaceChildren();const legendRows=[...new Map(main.filter(r=>r.comparison!=='sampling').map(r=>[r.label,r])).values()];for(const r of legendRows){const span=document.createElement('span'),swatch=document.createElement('i');swatch.className='swatch';swatch.style.background=colorOf(r);span.append(swatch,document.createTextNode(r.label+(r.comparison==='extrapolation'?'（虚线）':'（实线）')+'　'));$('waveform-legend').append(span)}$('waveform-source').textContent=JSON.stringify({data_sha256:contentHash,baseline_id:data.baseline_id,provenance:data.provenance,protocol:data.protocol,limits:data.limits},null,2);table();$('waveform-export').disabled=false;render();return true;
 }catch(error){$('waveform-status').textContent=error.message;return false}}
 const ready=load();return {render,ready,get data(){return data},get contentHash(){return contentHash}};
}
