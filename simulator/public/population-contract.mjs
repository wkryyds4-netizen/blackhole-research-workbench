export function posterior(data,key,{selection=true,prior='uniform',calibration=1}={}){
 const r=data.scenarios[key];if(!r)throw Error('未知假设组合');
 if(!Number.isFinite(calibration)||calibration<=0)throw Error('选择效率倍率需为正数');
 const f=data.f_grid,n=r.n_events,[a,b,aa,bb,ab]=r.selection;
 const logp=f.map((x,i)=>r.event_log_likelihood[i]-(selection?n*Math.log(x*a*calibration+(1-x)*b):0)+(prior==='central'?Math.log(x*(1-x)):prior==='edge'?-.5*Math.log(x*(1-x)):0));
 const max=Math.max(...logp),p=logp.map(x=>Math.exp(x-max)),sum=p.reduce((a,b)=>a+b,0);p.forEach((_,i)=>p[i]/=sum);
 function quantile(q){let c=0;for(let i=0;i<p.length;i++){const old=c;c+=p[i];if(c>=q)return i===0?f[0]:f[i-1]+(f[i]-f[i-1])*(q-old)/p[i];}return f.at(-1);}
 const mc=f.map((x,i)=>{const xi=x*a*calibration+(1-x)*b,second=x*x*aa*calibration**2+(1-x)**2*bb+2*x*(1-x)*ab*calibration,inj=xi*xi/second,sd=Math.sqrt(Math.max(0,r.mc_event_var[i])+(selection?n*n*Math.max(0,second/xi**2-1/r.total_generated):0));return {sd,inj,pass:r.min_event_ess[i]>=50&&(!selection||inj>=4*n)&&sd<=1}});
 const passing=p.reduce((s,w,i)=>s+(mc[i].pass?w:0),0),mean=f.reduce((s,x,i)=>s+x*p[i],0);
 return {p,density:p.map(x=>x*p.length),interval:[.05,.5,.95].map(quantile),passing,mean,mc};
}
