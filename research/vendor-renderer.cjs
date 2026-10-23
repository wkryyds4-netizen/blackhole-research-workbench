// Copy only the render dependency graph, including type imports and shader includes.
const fs=require('fs'),path=require('path'),crypto=require('crypto');
const source=path.resolve(__dirname,'blackhole-sim-upstream'),dest=path.resolve(__dirname,'../simulator/renderer-src/vendor/blackhole-sim');
const seen=new Set(),manifest=[];
function copy(relative){
 const p=path.resolve(source,relative);if(!p.startsWith(source+path.sep))throw Error('Outside source');
 if(seen.has(p))return;seen.add(p);const body=fs.readFileSync(p,'utf8');
 const out=path.join(dest,relative);fs.mkdirSync(path.dirname(out),{recursive:true});fs.writeFileSync(out,body);
 manifest.push({path:relative.replaceAll('\\','/'),sha256:crypto.createHash('sha256').update(body).digest('hex')});
 const refs=[...body.matchAll(/(?:from\s+|import\s*)['"](\.[^'"]+)['"]/g),...body.matchAll(/#include\s+(\.\S+?);/g)];
 for(const m of refs){let dep=path.resolve(path.dirname(p),m[1]);if(!fs.existsSync(dep))dep+='.ts';copy(path.relative(source,dep));}
}
for(const file of ['src/render/blackHolePass.ts','src/render/pipeline.ts','src/render/starfield.ts','LICENSE'])copy(file);
fs.writeFileSync(path.join(dest,'PROVENANCE.json'),JSON.stringify({repository:'https://github.com/algometrix/blackhole-sim',commit:'0ee5541ec4ffcd9ca708f37c080bf00b7f3ef6d1',license:'MIT',files:manifest},null,2));
console.log(`Copied ${manifest.length} files`);
