const {chromium}=require('playwright');
const fs=require('fs'),path=require('path');
(async()=>{const b=await chromium.launch({channel:process.env.PLAYWRIGHT_CHANNEL==='chromium'?undefined:(process.env.PLAYWRIGHT_CHANNEL||'msedge'),headless:true});const p=await b.newPage({viewport:{width:1440,height:1100}}),out=path.join(__dirname,'../research/stage3-results');const report={checks:[],errors:[]};const check=(n,v)=>{if(!v)throw Error(n);report.checks.push(n)};p.on('pageerror',e=>report.errors.push(e.message));
try{
 await p.goto('http://127.0.0.1:8765/workbench.html?view=waveform');await p.waitForFunction(()=>window.researchWorkbench?.waveform?.data);
 check('Waveform deep link opens correct view',await p.locator('#waveform-view').isVisible());
 check('All145 configurations loaded',await p.locator('#waveform-table tbody tr').count()===145);
 check('Scientific boundaries visible',(await p.locator('#waveform-status').innerText()).includes('训练集'));
 const reference=await p.locator('#waveform-baseline').innerText();const initial=await p.locator('#waveform-chart').evaluate(c=>c.toDataURL());
 await p.locator('#waveform-psd').selectOption('white');await p.locator('#waveform-fmin').selectOption('30');await p.locator('#waveform-taper').selectOption('0.1');
 check('Filters update plotted curves',initial!==await p.locator('#waveform-chart').evaluate(c=>c.toDataURL()));
 check('Waveform display remains fixed reference',reference===await p.locator('#waveform-baseline').innerText());
 const dl=p.waitForEvent('download');await p.locator('#waveform-export').click();await (await dl).saveAs(path.join(out,'waveform-session.json'));const saved=JSON.parse(fs.readFileSync(path.join(out,'waveform-session.json')));
 check('Export includes hash, filters, provenance and all records',saved.data_sha256.length===64&&saved.study.rows.length===145&&saved.view.psd==='white'&&saved.view.taper===.1&&saved.study.provenance.independent_validation===false);
 await p.locator('#waveform-psd').selectOption('aligo_historical');await p.locator('#waveform-fmin').selectOption('20');await p.locator('#waveform-taper').selectOption('0.02');
 await p.screenshot({path:path.join(out,'waveform-desktop.png'),fullPage:true});
 for(const width of [1024,736,390,320]){await p.setViewportSize({width,height:1000});check(`No overflow at width ${width}`,await p.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1));if(width===390)await p.screenshot({path:path.join(out,'waveform-mobile.png'),fullPage:true})}
 await p.setViewportSize({width:1440,height:1100});for(const view of ['dynamics','alignment','population','scene','projects','waveform']){await p.locator(`[data-view="${view}"]`).click();check(`Tab ${view} remains usable`,await p.locator(`#${view}-view`).isVisible())}
 check('No browser runtime errors',report.errors.length===0);
}finally{fs.writeFileSync(path.join(out,'browser-validation.json'),JSON.stringify(report,null,2));console.log(JSON.stringify(report,null,2));await b.close()}})().catch(e=>{console.error(e);process.exitCode=1});
