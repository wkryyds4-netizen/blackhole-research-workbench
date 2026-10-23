import * as THREE from 'three';
import {OrbitControls} from 'three/addons/controls/OrbitControls.js';
import {BlackHolePass} from './vendor/blackhole-sim/src/render/blackHolePass';
import {RenderPipeline} from './vendor/blackhole-sim/src/render/pipeline';
import {Starfield} from './vendor/blackhole-sim/src/render/starfield';
import {defaultSettings} from './vendor/blackhole-sim/src/settings';
import {sceneAt} from './scene-state.mjs';
import {innermostStableCircularOrbit} from './vendor/blackhole-sim/src/physics/kerr';

export function createLensing(canvas,{onError=()=>{},onStatus=()=>{}}={}){
 const renderer=new THREE.WebGLRenderer({canvas,antialias:false,alpha:false,powerPreference:'high-performance'});
 if(!renderer.extensions.has('EXT_color_buffer_float')){renderer.dispose();throw Error('设备不支持 HDR 浮点渲染');}
 renderer.setPixelRatio(Math.min(devicePixelRatio||1,1.5));
 const camera=new THREE.PerspectiveCamera(48,1,.05,200),renderCamera=camera.clone();
 const controls=new OrbitControls(camera,canvas);controls.enableDamping=true;controls.dampingFactor=.12;controls.enablePan=false;controls.minDistance=14;controls.maxDistance=45;controls.maxPolarAngle=Math.PI-.08;controls.minPolarAngle=.08;
 const sky=new Starfield(renderer,{...defaultSettings().sky,nebulaIntensity:.6},512);
 let quality='medium',disk=true,active=true,disposed=false,frame=null,wave=null,time=NaN,playing=false,scene=null,clockOrigin=0,accum=0,frames=0,width=0,height=0;
 const bh=new BlackHolePass(sky.texture,quality),pipeline=new RenderPipeline(renderer,renderCamera,bh,quality);
 const primary=new THREE.Vector3(),secondary=new THREE.Vector3();
 const diagnostics={frames:0,quality,accum:0,scene:null,active:true};
 function schedule(){if(!disposed&&active&&!document.hidden&&frame===null&&scene)frame=requestAnimationFrame(render);}
 function invalidate(){accum=0;pipeline.resetAccumulation();schedule();}
 function cameraView(top=false){const r=25,el=top?Math.PI/2-.03:.24,az=.45;camera.position.set(r*Math.sin(az)*Math.cos(el),r*Math.sin(el),r*Math.cos(az)*Math.cos(el));controls.target.set(0,0,0);camera.lookAt(controls.target);controls.update();invalidate();}
 function resize(){const r=canvas.getBoundingClientRect();if(r.width<1||r.height<1)return false;if(width!==r.width||height!==r.height){width=r.width;height=r.height;renderer.setSize(width,height,false);camera.aspect=width/height;camera.fov=THREE.MathUtils.radToDeg(2*Math.atan(Math.tan(THREE.MathUtils.degToRad(24))*Math.max(1,1.8/camera.aspect)));camera.updateProjectionMatrix();const size=renderer.getDrawingBufferSize(new THREE.Vector2());pipeline.setSize(size.x,size.y);invalidate();}return true;}
 function render(){
  frame=null;if(disposed||!active||document.hidden||!scene)return;
  try{
   if(!resize())return;
   const moving=controls.update();
   primary.fromArray(scene.primary);secondary.fromArray(scene.secondary).sub(primary);
   // Upstream's primary is at zero. Shift the camera and companion together;
   // the user camera remains centered on the binary barycenter.
   renderCamera.copy(camera);renderCamera.position.sub(primary);renderCamera.updateMatrixWorld(true);
   const spin=scene.binary?0:(scene.opticalSpin??0);
   bh.setPrimary({rs:scene.primaryRs,spin});
   bh.setSecondary(scene.binary?{pos:secondary,rs:scene.secondaryRs}:null);
   bh.setDisc({inner:innermostStableCircularOrbit(spin,'prograde')*scene.primaryRs,outer:8*scene.primaryRs,brightness:disk?.4:0});
   pipeline.render(scene.time-clockOrigin,{accumMode:playing||moving?'off':'progressive',bloomStrength:.65,fade:1});
   accum=playing||moving?0:accum+1;frames++;
   Object.assign(diagnostics,{frames,accum,scene,active,quality});
   if(playing||moving||accum<128)schedule();
  }catch(e){fail(e);}
 }
 function fail(e){dispose();onError(e);}
 function lost(e){e.preventDefault();fail(Error('WebGL 上下文丢失，已切换轨道示意'));}
 function dispose(){if(disposed)return;disposed=true;active=false;cancelAnimationFrame(frame);frame=null;canvas.removeEventListener('webglcontextlost',lost);document.removeEventListener('visibilitychange',visibility);observer.disconnect();controls.dispose();sky.dispose();if(!renderer.getContext().isContextLost())renderer.forceContextLoss();renderer.dispose();}
 function visibility(){if(document.hidden){cancelAnimationFrame(frame);frame=null;}else invalidate();}
 controls.addEventListener('change',invalidate);
 canvas.addEventListener('webglcontextlost',lost);
 document.addEventListener('visibilitychange',visibility);
 const observer=new ResizeObserver(invalidate);observer.observe(canvas);
 cameraView();onStatus('透镜就绪 · 两中心近似');
 return {
  diagnostics,
  update(state){const changed=wave!==state.wave||time!==state.time||playing!==state.playing;wave=state.wave;time=state.time;playing=state.playing;if(wave&&changed){clockOrigin=wave.t[0];scene=sceneAt(wave,time);invalidate();}},
  setFrame(value,options={}){if(value.binary&&value.opticalSpin)throw Error('Kerr binary optics are unsupported');if(!Number.isFinite(value.opticalSpin??0)||(value.opticalSpin??0)<0||(value.opticalSpin??0)>.998)throw Error('Optical spin out of range');scene=structuredClone(value);wave=null;clockOrigin=options.clockOrigin??0;playing=options.playing??false;invalidate();},
  setActive(value){if(active===value)return;active=value;diagnostics.active=value;if(!value){cancelAnimationFrame(frame);frame=null;}else invalidate();},
  setQuality(value){if(!['low','medium','high'].includes(value))return;quality=value;pipeline.setQuality(value);invalidate();},
  setDisk(value){disk=value;invalidate();},
  cameraView,dispose
 };
}
