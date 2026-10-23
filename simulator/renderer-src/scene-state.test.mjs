import test from 'node:test';
import assert from 'node:assert/strict';
import {sceneAt} from './scene-state.mjs';
const wave={t:[-.03,-.02,-.01,0,.01],sample_rate:100,parameters:{m1:60,m2:20},phase:[0,Math.PI,2*Math.PI,3*Math.PI,4*Math.PI],separation:[12,10,8,null,null]};
test('Inspiral respects mass-weighted barycenter and waveform phase',()=>{
 for(const t of [-.03,-.02]){const s=sceneAt(wave,t);for(let k=0;k<3;k++)assert.ok(Math.abs(60*s.primary[k]+20*s.secondary[k])<1e-12);assert.ok(Math.abs(Math.hypot(...s.primary.map((x,k)=>x-s.secondary[k]))-wave.separation[s.index])<1e-12);assert.equal(s.phase,wave.phase[s.index]/2);}
});
test('Backward scrubbing reproduces the same state without accumulated dynamics',()=>{const a=sceneAt(wave,-.02);sceneAt(wave,.01);assert.deepEqual(sceneAt(wave,-.02),a);});
test('Missing inspiral estimate and merger produce finite coordinates',()=>{for(const t of [-1,-.01,0,1]){const s=sceneAt(wave,t);assert.ok([...s.primary,...s.secondary,s.primaryRs,s.secondaryRs].every(Number.isFinite));}assert.equal(sceneAt(wave,0).binary,false);assert.equal(sceneAt(wave,0).opticalSpin,0);assert.equal(sceneAt(wave,-.01).illustrativeMerger,true);});
