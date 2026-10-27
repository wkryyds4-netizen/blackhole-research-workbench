import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {validateMixtureStudy,selectMixtureCase} from './public/mixture-study-contract.mjs';
const root=new URL('../',import.meta.url),read=p=>readFileSync(new URL(p,root));
const sha=b=>createHash('sha256').update(b).digest('hex');
const bytes=read('simulator/public/datasets/channel-study.json'),d=validateMixtureStudy(JSON.parse(bytes));
const manifest=JSON.parse(read('research/population-channels/workbench-export.json'));
assert.equal(sha(bytes),manifest.sha256);
for(const p of manifest.source_reports)assert.equal(sha(read(p.repo_path)),p.sha256);
for(const p of manifest.papers)assert.equal(sha(read('simulator/public/'+p.url)),p.sha256);
assert.equal(sha(read('research/population-channels/export_workbench.py')),manifest.exporter_sha256);
assert.equal(d.cases.length,24);assert.equal(d.observations.events,69);assert.equal(d.observations.posterior_rows,11180578);
for(const r of d.cases)assert.equal(selectMixtureCase(d,r.settings),r);
const fit=JSON.parse(read('research/population-channels/branching-results.json'));
const terms=JSON.parse(read('research/population-channels/detected-component-results.json'));
for(const r of d.cases){
 assert.deepEqual(r.quantiles.intrinsic,fit.cases[r.id].intrinsic_quantiles);
 assert.deepEqual(r.quantiles.detected_label,fit.cases[r.id].detected_quantiles);
 assert.deepEqual(r.quantiles.detected_model,terms.cases[r.id].model_term_quantiles);
 assert.deepEqual(r.common_term,terms.cases[r.id].uniform_term_quantiles);
}
// Marginal medians are not a composition; enforcing their sum would be wrong.
const r=selectMixtureCase(d,d.defaults);assert.ok(Math.abs(r.quantiles.intrinsic.reduce((s,q)=>s+q[1],0)-1)>.001);
assert.ok(Math.abs(r.intrinsic_mean.reduce((s,x)=>s+x,0)-1)<1e-8);
const otherInterval=structuredClone(d);otherInterval.quantile_probabilities=[.025,.5,.975];validateMixtureStudy(otherInterval);
for(const mutate of [x=>x.value_unit='percentage_point',x=>x.quantile_probabilities=[.05,.6,.95],x=>x.cases.pop(),x=>x.cases[0].intrinsic_mean[0]=0,x=>x.cases[0].quantiles.intrinsic[0]=[.9,.8,.7],x=>x.cases[0].common_term_definition='',x=>x.papers[0].url='javascript:alert(1)']){
 const copy=structuredClone(d);mutate(copy);assert.throws(()=>validateMixtureStudy(copy));
}
console.log('Mixture study:24 grid points, quantity semantics, means vs medians, provenance and safe paper paths passed.');
