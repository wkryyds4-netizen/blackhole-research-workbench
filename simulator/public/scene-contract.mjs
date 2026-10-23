// Explicit data contract. Imported files remain local; no code or URLs are executed.
const finite=x=>typeof x==='number'&&Number.isFinite(x);
const vec=v=>Array.isArray(v)&&v.length===3&&v.every(finite);
const fail=s=>{throw Error(s)};
export function validateScene(d){
 if(d?.schema!=='blackhole-scene/1')fail('需要 blackhole-scene/1 格式');
 if(typeof d.title!=='string'||!d.title.trim()||d.title.length>160)fail('标题需为 1–160 字符');
 if(!['s','GM/c^3'].includes(d.units?.time)||d.units?.length!=='GM/c^2'||!finite(d.units?.reference_mass_msun)||d.units.reference_mass_msun<=0||!['source','detector','geometric'].includes(d.units?.mass_frame))fail('需要显式时间、长度、参考质量及质量参考系');
 if(!['analytic','numerical','illustration','waveform-derived'].includes(d.provenance?.kind)||typeof d.provenance?.source!=='string'||!d.provenance.source.trim()||typeof d.provenance?.frame!=='string'||!d.provenance.frame.trim()||typeof d.provenance?.version!=='string'||!d.provenance.version.trim())fail('需要来源类别、来源版本与坐标系说明');
 if(!Array.isArray(d.frames)||d.frames.length<1||d.frames.length>30000)fail('帧数必须为 1–30000');
 let prev=-Infinity;
 for(const f of d.frames){
  if(!finite(f.t)||f.t<=prev)fail('时间必须有限且严格递增');prev=f.t;
  if(!Array.isArray(f.bodies)||f.bodies.length<1||f.bodies.length>12)fail('每帧需有 1–12 个天体');
  const ids=new Set();for(const b of f.bodies){if(typeof b.id!=='string'||!b.id||ids.has(b.id))fail('天体 ID 需唯一');ids.add(b.id);if(!finite(b.mass)||b.mass<=0||!vec(b.position)||!vec(b.spin)||Math.hypot(...b.spin)>1+1e-9)fail('天体质量、坐标或无量纲自旋无效');}
 }
 if(d.provenance.source_common_horizon_time!==undefined&&!finite(d.provenance.source_common_horizon_time))fail('共同视界原始时刻必须有限');
 if(d.provenance.gaps!==undefined){if(!Array.isArray(d.provenance.gaps)||d.provenance.gaps.some(g=>!finite(g?.start)||!finite(g?.end)||g.end<=g.start||g.start<d.frames[0].t||g.end>d.frames.at(-1).t))fail('采样缺口必须是场景范围内的有限时间区间');}
 return structuredClone(d);
}
export function sampleScene(d,t){
 const frames=d.frames;let lo=0,hi=frames.length-1;
 while(lo<hi){const mid=Math.ceil((lo+hi)/2);if(frames[mid].t<=t)lo=mid;else hi=mid-1;}
 const a=frames[lo],b=frames[Math.min(lo+1,frames.length-1)];
 // Topology changes (e.g. common horizon) are discrete; never morph IDs.
 const same=a.bodies.length===b.bodies.length&&a.bodies.every((x,i)=>x.id===b.bodies[i].id);
 const u=same&&b.t>a.t?Math.max(0,Math.min(1,(t-a.t)/(b.t-a.t))):0;
 return {t:Math.max(frames[0].t,Math.min(t,frames.at(-1).t)),bodies:a.bodies.map((x,i)=>{const y=same?b.bodies[i]:x;return {id:x.id,mass:x.mass*(1-u)+y.mass*u,position:x.position.map((v,k)=>v*(1-u)+y.position[k]*u),spin:x.spin.map((v,k)=>v*(1-u)+y.spin[k]*u)}})};
}
export function opticalCapability(f){
 if(f.bodies.length>2)return {supported:false,reason:'当前透镜仅支持一或两个天体；使用通用坐标视图'};
 if(f.bodies.some(b=>b.mass<.15||b.mass>1.01))return {supported:false,reason:'半径超出当前透镜积分尺度；使用坐标视图'};
 if(f.bodies.some(b=>Math.hypot(...b.position)>22))return {supported:false,reason:'轨迹超出当前光线积分场景范围；使用坐标视图'};
 if(f.bodies.length===2){
  if(f.bodies.some(b=>Math.hypot(...b.spin)>1e-8))return {supported:false,reason:'两个自旋黑洞不在光学模型能力内；保留自旋数据并显示坐标'};
  const a=f.bodies[0],b=f.bodies[1];const sep=Math.hypot(...a.position.map((v,k)=>v-b.position[k]));
  if(sep>25||sep<2*(a.mass+b.mass))return {supported:false,reason:'分离度超出两中心显示范围；使用坐标视图'};
  return {supported:true,mode:'双中心近似',spin:0};
 }
 const s=f.bodies[0].spin;
 if(Math.abs(s[0])+Math.abs(s[2])>1e-8||s[1]>1e-8||-s[1]>.998)return {supported:false,reason:'当前 Kerr 渲染轴固定为 −y，且 χ≤0.998；使用坐标视图'};
 return {supported:true,mode:-s[1]>0?'单 Kerr':'单 Schwarzschild',spin:-s[1]};
}
export function opticalFrame(f){
 const c=opticalCapability(f);if(!c.supported)throw Error(c.reason);
 const [a,b]=f.bodies;
 return {time:f.t,primary:a.position,secondary:b?.position??a.position,primaryRs:2*a.mass,secondaryRs:b?2*b.mass:0,binary:!!b,opticalSpin:c.spin};
}
export function singleScene(spin=.7){return {schema:'blackhole-scene/1',title:spin?'单 Kerr 黑洞 · 定常光学实验':'单 Schwarzschild 黑洞 · 定常光学实验',units:{time:'GM/c^3',length:'GM/c^2',reference_mass_msun:65,mass_frame:'geometric'},provenance:{kind:'analytic',source:'Kerr metric; algometrix/blackhole-sim optical integrator',version:'0ee5541ec4ffcd9ca708f37c080bf00b7f3ef6d1',frame:'Cartesian Kerr-Schild; spin along -y; illustrative emissive disk'},frames:[{t:0,bodies:[{id:'BH',mass:1,position:[0,0,0],spin:[0,-spin,0]}]}]};}
export function orbitExample(spinning=false){
 const frames=[],q=.6,a=10,e=.25,M1=1/(1+q),M2=q/(1+q);
 for(let j=0;j<=160;j++){const t=j/160*2*Math.PI*Math.sqrt(a**3),mean=j/160*2*Math.PI;let E=mean;for(let k=0;k<8;k++)E-=(E-e*Math.sin(E)-mean)/(1-e*Math.cos(E));const x=a*(Math.cos(E)-e),z=a*Math.sqrt(1-e*e)*Math.sin(E);frames.push({t,bodies:[{id:'A',mass:M1,position:[x*M2,0,z*M2],spin:spinning?[.3,-.2,.1]:[0,0,0]},{id:'B',mass:M2,position:[-x*M1,0,-z*M1],spin:[0,0,0]}]});}
 return {schema:'blackhole-scene/1',title:spinning?'自旋向量接口演示 · 自动使用坐标视图':'牛顿偏心轨道 · 接口演示',units:{time:'GM/c^3',length:'GM/c^2',reference_mass_msun:65,mass_frame:'geometric'},provenance:{kind:'illustration',source:'Newtonian Kepler ellipse, e=0.25; no radiation reaction; not an eccentric GR waveform',version:'kepler-demo/1',frame:'Barycentric Cartesian, orbital plane xz; linear display interpolation'},frames};
}
