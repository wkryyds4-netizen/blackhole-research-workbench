const {chromium}=require('playwright');
const fs=require('fs'),path=require('path');
(async()=>{
 const browser=await chromium.launch({channel:process.env.PLAYWRIGHT_CHANNEL==='chromium'?undefined:(process.env.PLAYWRIGHT_CHANNEL||'msedge'),headless:true});const p=await browser.newPage({viewport:{width:1440,height:1100}});
 const report={checks:[],errors:[]},out=path.join(__dirname,'../research/adapter-results');fs.mkdirSync(out,{recursive:true});
 const check=(name,condition)=>{if(!condition)throw Error(name);report.checks.push(name)};p.on('pageerror',e=>report.errors.push(e.message));
 try{
  await p.goto('http://127.0.0.1:8765/workbench.html?dataset=sxs-0003-lev5');
  await p.waitForFunction(()=>window.researchWorkbench?.scene.title.includes('0003'));
  check('All seven local SXS datasets registered',await p.locator('#generic-preset optgroup option').count()===7);
  check('Precessing spin preserved and optical fallback enforced',await p.evaluate(()=>Math.hypot(...researchWorkbench.frame.bodies[0].spin)>.49&&!researchWorkbench.capability.supported));
  check('Deprecated source visibly marked',(await p.locator('#scene-data-note').innerText()).includes('2327'));
  await p.locator('#generic-scrub').fill('680');await p.locator('#generic-scrub').dispatchEvent('input');
  await p.locator('#scene-quantity').selectOption('spin');
  check('Time changes actual spin vector',await p.evaluate(()=>JSON.stringify(researchWorkbench.frame.bodies[0].spin)!==JSON.stringify(researchWorkbench.scene.frames[0].bodies[0].spin)));
  await p.screenshot({path:path.join(out,'sxs-precession.png'),fullPage:true});
  await p.locator('#generic-merger').click();
  check('Common horizon is an exact discrete C transition',await p.evaluate(()=>researchWorkbench.frame.t===0&&researchWorkbench.frame.bodies.length===1&&researchWorkbench.frame.bodies[0].id==='C'));
  await p.locator('#scene-quantity').selectOption('mass');await p.screenshot({path:path.join(out,'sxs-common-horizon.png'),fullPage:true});
  const dl=p.waitForEvent('download');await p.locator('#generic-export').click();const file=await dl;await file.saveAs(path.join(out,'scene-roundtrip.json'));
  const saved=JSON.parse(fs.readFileSync(path.join(out,'scene-roundtrip.json')));check('Export retains source hashes, license and coordinate convention',saved.provenance.source_files.every(f=>f.sha256.length===64)&&saved.provenance.license==='CC-BY-4.0'&&saved.provenance.rotation_matrix.length===3);
  await p.locator('#generic-preset').selectOption('sxs-0389-lev3');await p.waitForFunction(()=>researchWorkbench.scene.title.includes('0389'));
  check('Resolution preset switches dataset',await p.evaluate(()=>researchWorkbench.scene.provenance.version.includes('Lev3')));
  await p.locator('#generic-preset').selectOption('sxs-0001-lev5');await p.waitForFunction(()=>researchWorkbench.scene.title.includes('0001'));
  await p.evaluate(()=>{const s=researchWorkbench.scene,g=s.provenance.gaps[0],t=(g.start+g.end)/2,e=document.getElementById('generic-scrub');e.value=(t-s.frames[0].t)/(s.frames.at(-1).t-s.frames[0].t)*1000;e.dispatchEvent(new Event('input'))});
  check('Native gap reported without inventing merger interpolation',(await p.locator('#scene-data-note').innerText()).includes('采样缺口')&&await p.evaluate(()=>researchWorkbench.frame.bodies.length===2));
  await p.locator('#generic-merger').click();check('Gap ends at exact common horizon',await p.evaluate(()=>researchWorkbench.frame.bodies[0].id==='C'));
  await p.route('**/datasets/sxs-0389-lev2.json',route=>route.fulfill({status:503,body:'unavailable'}));
  await p.locator('#generic-preset').selectOption('sxs-0389-lev2');await p.waitForFunction(()=>document.getElementById('import-status').textContent.includes('失败'));
  check('Failed fetch preserves prior valid scene',await p.evaluate(()=>researchWorkbench.scene.title.includes('0001')&&document.getElementById('generic-preset').value==='sxs-0001-lev5'));
  for(const width of [390,320]){await p.setViewportSize({width,height:1000});check('Numerical dataset mobile layout '+width,await p.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1))}
  await p.screenshot({path:path.join(out,'sxs-mobile.png'),fullPage:true});check('No page errors',report.errors.length===0);
 }finally{fs.writeFileSync(path.join(out,'browser-validation.json'),JSON.stringify(report,null,2));console.log(JSON.stringify(report,null,2));await browser.close()}
})().catch(e=>{console.error(e);process.exitCode=1});

