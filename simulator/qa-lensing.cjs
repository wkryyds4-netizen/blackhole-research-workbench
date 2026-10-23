const {chromium}=require('playwright');
const fs=require('fs'),path=require('path');
(async()=>{
 const b=await chromium.launch({channel:process.env.PLAYWRIGHT_CHANNEL==='chromium'?undefined:(process.env.PLAYWRIGHT_CHANNEL||'msedge'),headless:true});const page=await b.newPage({viewport:{width:1440,height:1150}});
 const report={checks:[],errors:[],warnings:[],httpFailures:[]};page.on('response',r=>{if(r.status()>=400)report.httpFailures.push({url:r.url(),status:r.status()})});page.on('pageerror',e=>report.errors.push(e.message));page.on('console',m=>{if(m.type()==='error')report.errors.push(m.text()+' '+m.location().url);if(m.type()==='warning')report.warnings.push(m.text())});
 const dir=path.join(__dirname,'qa-output');fs.mkdirSync(dir,{recursive:true});
 const ok=(name,v)=>{if(!v)throw Error(name);report.checks.push(name)};
 await page.goto('http://127.0.0.1:8765/');await page.waitForFunction(()=>bbhLab.state.wave&&!bbhLab.state.busy&&bbhLab.lensing?.diagnostics.frames>3,{},{timeout:60000});
 await page.waitForTimeout(800);await page.screenshot({path:path.join(dir,'lensing-desktop.png'),fullPage:true});
 const original=await page.evaluate(()=>JSON.stringify(bbhLab.state.wave));
 for(const q of ['low','high','medium']){await page.locator('#optical-quality').selectOption(q);await page.waitForTimeout(180);ok('Quality '+q,await page.evaluate(q=>bbhLab.lensing.diagnostics.quality===q,q));}
 await page.locator('#optical-disk').uncheck();await page.locator('#camera-top').click();await page.waitForTimeout(400);await page.screenshot({path:path.join(dir,'lensing-no-disk.png')});
 ok('Optical controls preserve physical arrays',original===await page.evaluate(()=>JSON.stringify(bbhLab.state.wave)));
 ok('Camera leaves physical inclination intact',await page.locator('#inclination').inputValue()==='35');
 await page.locator('#optical-disk').check();await page.locator('#camera-oblique').click();
 await page.locator('#scrub').fill('450');await page.locator('#scrub').dispatchEvent('input');await page.waitForTimeout(150);
 const first=await page.evaluate(()=>JSON.stringify(bbhLab.lensing.diagnostics.scene));
 ok('Single clock drives optics',await page.evaluate(()=>Math.abs(bbhLab.lensing.diagnostics.scene.time-bbhLab.state.time)<=1/bbhLab.state.wave.sample_rate));
 await page.locator('#scrub').fill('1000');await page.locator('#scrub').dispatchEvent('input');await page.waitForTimeout(200);ok('Merger produces one remnant',await page.evaluate(()=>bbhLab.lensing.diagnostics.scene.binary===false));await page.screenshot({path:path.join(dir,'lensing-remnant.png')});
 await page.locator('#scrub').fill('450');await page.locator('#scrub').dispatchEvent('input');await page.waitForTimeout(200);ok('Backward seek is deterministic',first===await page.evaluate(()=>JSON.stringify(bbhLab.lensing.diagnostics.scene)));
 await page.locator('#scene-mode').selectOption('schematic');await page.waitForTimeout(100);const n=await page.evaluate(()=>bbhLab.lensing.diagnostics.frames);await page.waitForTimeout(250);ok('Inactive renderer stops',n===await page.evaluate(()=>bbhLab.lensing.diagnostics.frames));ok('Canvas fallback visible',await page.locator('#orbit').isVisible());
 await page.locator('#scene-mode').selectOption('lensing');await page.waitForFunction(()=>bbhLab.lensing.diagnostics.accum>=128,{},{timeout:30000});const settled=await page.evaluate(()=>bbhLab.lensing.diagnostics.frames);await page.waitForTimeout(250);ok('Paused accumulation converges and stops',settled===await page.evaluate(()=>bbhLab.lensing.diagnostics.frames));
 await page.locator('#play').click();await page.waitForTimeout(250);await page.locator('#play').click();ok('Playback remains synchronized',await page.evaluate(()=>Math.abs(bbhLab.lensing.diagnostics.scene.time-bbhLab.state.time)<.01));
 for(const width of [1024,736,390,320]){await page.setViewportSize({width,height:1100});await page.waitForTimeout(250);ok('Layout '+width,await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1));if(width===390)await page.screenshot({path:path.join(dir,'lensing-mobile.png'),fullPage:true});}
 await page.setViewportSize({width:1440,height:1150});await page.locator('#scrub').fill('680');await page.locator('#scrub').dispatchEvent('input');await page.waitForTimeout(350);
 await page.screenshot({path:path.join(dir,'lensing-desktop.png'),fullPage:true});
 // Simulate the actual WebGL context loss event; application must recover.
 await page.evaluate(()=>{const gl=document.querySelector('#lensing').getContext('webgl2');gl.getExtension('WEBGL_lose_context').loseContext()});await page.waitForTimeout(300);
 ok('Context loss recovers to schematic',await page.locator('#orbit').isVisible()&&(await page.locator('#scene-mode').inputValue())==='schematic');
 fs.writeFileSync(path.join(dir,'lensing-validation.json'),JSON.stringify(report,null,2));
 console.log(JSON.stringify({errors:report.errors,httpFailures:report.httpFailures}));
 ok('No script, HTTP or WebGL errors',report.errors.length===0&&report.httpFailures.length===0);
 fs.writeFileSync(path.join(dir,'lensing-validation.json'),JSON.stringify(report,null,2));console.log(JSON.stringify(report,null,2));await b.close();
})().catch(e=>{console.error(e);process.exit(1)});
