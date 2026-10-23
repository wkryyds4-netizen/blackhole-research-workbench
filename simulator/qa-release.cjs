// Public-distribution smoke test; never requires the excluded local datasets.
const {chromium}=require('playwright');
const assert=require('node:assert/strict');
(async()=>{
 const browser=await chromium.launch({channel:process.env.PLAYWRIGHT_CHANNEL==='chromium'?undefined:(process.env.PLAYWRIGHT_CHANNEL||'msedge'),headless:true});
 const page=await browser.newPage({viewport:{width:1280,height:900}}),errors=[];
 page.on('pageerror',e=>errors.push(e.message));
 const base=process.env.WORKBENCH_URL||'http://127.0.0.1:8765';
 try {
  await page.route('**/datasets/catalog.json',r=>r.fulfill({status:404,body:'Local catalog excluded from public release'}));
  await page.goto(base+'/workbench.html?view=waveform');
  await page.waitForFunction(()=>window.researchWorkbench?.waveform?.data&&document.querySelectorAll('#generic-preset optgroup option').length===4);
  assert.equal(await page.locator('#waveform-table tbody tr').count(),145);
  assert.equal(await page.locator('#generic-preset option[value^="sxs-2327"]').count(),0);
  await page.locator('[data-view="scene"]').click();
  const entries=await page.locator('#generic-preset optgroup option').evaluateAll(o=>o.map(e=>({id:e.value,title:e.textContent})));
  for(const entry of entries){await page.locator('#generic-preset').selectOption(entry.id);await page.waitForFunction(title=>window.researchWorkbench.scene.title===title,entry.title);assert.ok((await page.locator('#import-status').innerText()).includes('已载入'));}
  await page.locator('[data-view="population"]').click();
  await page.waitForFunction(()=>window.researchWorkbench.population!==null);
  for(const view of ['alignment','dynamics','waveform']){await page.locator(`[data-view="${view}"]`).click();assert.ok(await page.locator(`#${view}-view`).isVisible());}
  await page.setViewportSize({width:390,height:900});assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1));
  const download=page.waitForEvent('download');await page.locator('#waveform-export').click();assert.ok((await download).suggestedFilename().endsWith('.json'));
  assert.deepEqual(errors,[]);console.log('Public release passed: four redistributable scenes,145 waveform rows,population,tabs,export,mobile layout; no page errors.');
 } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1});
