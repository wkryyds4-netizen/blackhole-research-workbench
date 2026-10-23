"""Fetch small public audit documents only; never load pickle or large samples."""
from pathlib import Path
import hashlib
import json
from html.parser import HTMLParser
import requests

OUT = Path(__file__).resolve().parent / 'provenance-sources'
OUT.mkdir(exist_ok=True)
SOURCES = {
    'GWTC3-tutorial.html': 'https://colmtalbot.github.io/gwpopulation/examples/GWTC3.html',
    'injection-release.md': 'https://zenodo.org/records/7890398/files/o1%2Bo2%2Bo3-sensitivity-estimates.md?download=1',
    'vt_helper.html': 'https://ratesandpopulations.docs.ligo.org/gwpopulation_pipe/_modules/gwpopulation_pipe/vt_helper.html',
    'data_collection.html': 'https://ratesandpopulations.docs.ligo.org/gwpopulation_pipe/_modules/gwpopulation_pipe/data_collection.html',
    'GWTC-1-live.json': 'https://gwosc.org/eventapi/json/GWTC-1-confident/',
}

class PreParser(HTMLParser):
    def __init__(self):
        super().__init__(); self.depth=0; self.parts=[]
    def handle_starttag(self, tag, attrs):
        if tag=='pre': self.depth += 1
    def handle_endtag(self, tag):
        if tag=='pre': self.depth -= 1
    def handle_data(self, data):
        if self.depth: self.parts.append(data)

manifest=[]
for filename,url in SOURCES.items():
    try:
        response=requests.get(url,timeout=50); response.raise_for_status()
    except requests.RequestException as exc:
        manifest.append(dict(file=filename,url=url,error=type(exc).__name__,retrieved='2026-10-24'))
        print(filename,type(exc).__name__,flush=True)
        continue
    data=response.content
    (OUT/filename).write_bytes(data)
    manifest.append(dict(file=filename,url=url,final_url=response.url,bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),retrieved='2026-10-24'))
    if filename.endswith('.html'):
        parser=PreParser(); parser.feed(response.text)
        (OUT/(filename+'.txt')).write_text(''.join(parser.parts),encoding='utf-8')
    print(filename,len(data),flush=True)
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
