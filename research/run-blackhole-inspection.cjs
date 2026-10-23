const {chromium}=require('playwright');
const fs=require('fs'),path=require('path');
(async()=>{
 const out=path.join(__dirname,'blackhole-sim-analysis');fs.mkdirSync(out,{recursive:true});
 const browser=await chromium.launch({headless:true,channel:process.env.PLAYWRIGHT_CHANNEL==='chromium'?undefined:(process.env.PLAYWRIGHT_CHANNEL||'msedge')});
 const page=await browser.newPage({viewport:{width:1100,height:760},deviceScaleFactor:1});
 const messages=[];page.on('pageerror',e=>messages.push({type:'pageerror',text:e.message}));page.on('console',m=>{if(m.type()==='error'||m.type()==='warning')messages.push({type:m.type(),text:m.text().slice(0,2000)});});
 await page.goto('http://127.0.0.1:5174/',{waitUntil:'networkidle',timeout:60000});
 await page.locator('#app canvas').waitFor({timeout:60000});
 await page.waitForTimeout(4000);
 const meta=await page.evaluate(()=>{const canvas=document.querySelector('#app canvas'),gl=canvas.getContext('webgl2'),ext=gl.getExtension('WEBGL_debug_renderer_info');return {title:document.title,renderer:ext?gl.getParameter(ext.UNMASKED_RENDERER_WEBGL):gl.getParameter(gl.RENDERER),vendor:ext?gl.getParameter(ext.UNMASKED_VENDOR_WEBGL):gl.getParameter(gl.VENDOR),size:[canvas.width,canvas.height],hud:document.querySelector('#hud')?.innerText,buttons:[...document.querySelectorAll('button')].map(b=>b.textContent),selects:[...document.querySelectorAll('select')].map(s=>({value:s.value,options:[...s.options].map(o=>({value:o.value,text:o.text}))})),text:document.body.innerText};});
 await page.screenshot({path:path.join(out,'upstream-default.png')});
 fs.writeFileSync(path.join(out,'initial-browser.json'),JSON.stringify({meta,messages},null,2));
 console.log(JSON.stringify({meta,messages},null,2));await browser.close();
})().catch(e=>{console.error(e);process.exit(1)});
