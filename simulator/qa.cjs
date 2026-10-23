const {chromium}=require('playwright');
const fs=require('fs'),path=require('path');
(async()=>{
 const out=path.join(__dirname,'qa-output');fs.mkdirSync(out,{recursive:true});
 const browser=await chromium.launch({headless:true,channel:process.env.PLAYWRIGHT_CHANNEL==='chromium'?undefined:(process.env.PLAYWRIGHT_CHANNEL||'msedge')});
 const page=await browser.newPage({viewport:{width:1440,height:1100},reducedMotion:'reduce'});
 const errors=[];page.on('pageerror',e=>errors.push(e.message));
 await page.goto('http://127.0.0.1:8765');
 await page.waitForFunction(()=>window.bbhLab?.state.wave&&!window.bbhLab.state.busy);
 const report={checks:[],errors};
 function ok(name,value){if(!value)throw Error(name);report.checks.push(name);}
 const baseline=await page.evaluate(()=>({peak:bbhLab.state.wave.meta.peak_strain,p:bbhLab.state.wave.parameters,time:bbhLab.state.time}));
 await page.screenshot({path:path.join(out,'desktop.png'),fullPage:true});
 await page.locator('#save-reference').click();ok('Reference saved',await page.locator('#clear-reference').isEnabled());
 await page.locator('#distance').fill('880');await page.locator('#distance').dispatchEvent('input');
 await page.waitForFunction(()=>!bbhLab.state.busy&&bbhLab.state.wave.parameters.distance===880);
 const peak=await page.evaluate(()=>bbhLab.state.wave.meta.peak_strain);ok('Distance changes physical amplitude',Math.abs(peak*2/baseline.peak-1)<1e-6);
 await page.locator('#scrub').fill('850');await page.locator('#scrub').dispatchEvent('input');
 const time=await page.locator('#scene-time').textContent();ok('Scrub updates physical clock',time.includes('t ='));
 await page.locator('#restart').click();await page.locator('#play').click();await page.waitForTimeout(350);await page.locator('#play').click();
 ok('Playback advances time',await page.evaluate(()=>bbhLab.state.time>bbhLab.state.wave.t[0]));
 await page.locator('#model').selectOption('IMRPhenomXAS');await page.waitForFunction(()=>!bbhLab.state.busy&&bbhLab.state.wave.parameters.model==='IMRPhenomXAS');ok('Second model computes',true);
 await page.locator('#preset').selectOption('equal');await page.waitForFunction(()=>!bbhLab.state.busy&&bbhLab.state.wave.parameters.m1===30);ok('Equal-mass boundary is finite',await page.evaluate(()=>bbhLab.state.wave.hp.every(Number.isFinite)));
 await page.locator('#model').selectOption('IMRPhenomD');await page.waitForFunction(()=>!bbhLab.state.busy&&bbhLab.state.wave.parameters.model==='IMRPhenomD');
 await page.locator('#preset').selectOption('default');await page.waitForFunction(()=>!bbhLab.state.busy&&bbhLab.state.wave.parameters.m1===36);
 // Exercise a stale response arriving inside the slider debounce interval.
 await page.route('**/api/waveform',async route=>{const body=route.request().postDataJSON();if(body.distance===600)await new Promise(r=>setTimeout(r,180));await route.continue();});
 await page.evaluate(()=>{bbhLab.state.p.distance=600;bbhLab.updateWave();});
 await page.waitForTimeout(30);await page.locator('#distance').fill('700');await page.locator('#distance').dispatchEvent('input');
 await page.waitForFunction(()=>!bbhLab.state.busy&&bbhLab.state.wave.parameters.distance===700);
 ok('Stale responses cannot overwrite new slider input',await page.evaluate(()=>bbhLab.state.p.distance===700));await page.unroute('**/api/waveform');
 await page.locator('#preset').selectOption('default');await page.waitForFunction(()=>!bbhLab.state.busy&&bbhLab.state.wave.parameters.distance===440);
 await page.locator('#audio').click();await page.waitForTimeout(100);ok('Audio playback starts from user action',(await page.locator('#audio').textContent()).includes('停止'));await page.locator('#audio').click();
 const downloadPromise=page.waitForEvent('download');await page.locator('#export').click();const download=await downloadPromise;await download.saveAs(path.join(out,'export.json'));const data=JSON.parse(fs.readFileSync(path.join(out,'export.json')));ok('Export contains waveform and provenance',data.hp.length>1000&&data.meta.version==='0.2.1');
 await page.locator('[data-tab="observation"]').click();await page.waitForFunction(()=>bbhLab.state.event);
 await page.locator('#fit').click();await page.waitForFunction(()=>bbhLab.state.fit);ok('H1 fit computes',await page.evaluate(()=>bbhLab.state.fit.correlation>0.5));
 await page.screenshot({path:path.join(out,'observation.png'),fullPage:true});
 await page.locator('#detector').selectOption('L1');await page.locator('#fit').click();await page.waitForFunction(()=>bbhLab.state.fit?.detector==='L1');ok('L1 fit computes',await page.evaluate(()=>bbhLab.state.fit.correlation>0.5));
 await page.locator('[data-tab="methods"]').click();ok('Methods visible',await page.locator('#methods').isVisible());
 await page.locator('[data-tab="simulation"]').click();await page.locator('#clear-reference').click();
 report.layouts=[];
 for(const width of [1440,1024,736,390,320]){
  await page.setViewportSize({width,height:1100});await page.waitForTimeout(200);
  const r=await page.evaluate(()=>({scroll:document.documentElement.scrollWidth,width:innerWidth,canvases:[...document.querySelectorAll('canvas')].filter(c=>c.getBoundingClientRect().width>0).map(c=>({id:c.id,width:c.getBoundingClientRect().width,height:c.getBoundingClientRect().height}))}));
  report.layouts.push(r);ok('No horizontal overflow at '+width,r.scroll<=width+1);ok('Polarization remains available at '+width,r.canvases.some(c=>c.id==='polar'));
  if(width===390)await page.screenshot({path:path.join(out,'mobile.png'),fullPage:true});
 }
 await page.setViewportSize({width:1440,height:1100});await page.locator('#scrub').fill('900');await page.locator('#scrub').dispatchEvent('input');await page.screenshot({path:path.join(out,'desktop-final.png'),fullPage:true});
 ok('No browser script errors',errors.length===0);fs.writeFileSync(path.join(out,'validation.json'),JSON.stringify(report,null,2));await browser.close();console.log(JSON.stringify(report,null,2));
})().catch(e=>{console.error(e);process.exit(1)});
