// Lengths: GM_total/c². Positions describe a visualization, not NR horizons.
export function sceneAt(wave,time){
 const i=Math.max(0,Math.min(wave.t.length-1,Math.round((time-wave.t[0])*wave.sample_rate)));
 const t=wave.t[i],p=wave.parameters,M=p.m1+p.m2,phase=wave.phase[i]/2;
 let separation=wave.separation[i];
 if(!Number.isFinite(separation)){let j=i;while(j>=0&&!Number.isFinite(wave.separation[j]))j--;separation=j>=0?wave.separation[j]:12;}
 // The last 12 ms use an explicitly illustrative transition. No remnant fit.
 const x=Math.max(0,Math.min(1,(t+.012)/.012)),merge=x*x*(3-2*x);
 const r=separation*(1-merge),dx=r*Math.cos(phase),dz=r*Math.sin(phase);
 const primary=[dx*p.m2/M,0,dz*p.m2/M],secondary=[-dx*p.m1/M,0,-dz*p.m1/M];
 return {index:i,time:t,phase,separation:r,primary,secondary,
  primaryRs:(2*p.m1/M)*(1-merge)+2*merge,
  secondaryRs:(2*p.m2/M)*(1-merge),binary:t<0&&merge<.999,
  illustrativeMerger:t>=-.012,opticalSpin:0};
}
