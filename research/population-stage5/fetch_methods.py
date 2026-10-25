"""Save primary catalog method papers for bounded text inspection."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import requests

BASE=Path(__file__).resolve().parent
SOURCES={'gwtc21-v2.pdf':'https://arxiv.org/pdf/2108.01045v2',
         'gwtc3-v3.pdf':'https://arxiv.org/pdf/2111.03606v3'}

def fetch(item):
    name,url=item;dest=BASE/'sources'/name
    if not dest.exists():
        r=requests.get(url,timeout=(15,45));r.raise_for_status()
        assert r.content.startswith(b'%PDF') and len(r.content)<25_000_000
        dest.write_bytes(r.content)
    data=dest.read_bytes()
    return dict(file=name,url=url,bytes=len(data),sha256=hashlib.sha256(data).hexdigest())

if __name__=='__main__':
    (BASE/'sources').mkdir(exist_ok=True)
    with ThreadPoolExecutor(max_workers=2) as pool:records=list(pool.map(fetch,SOURCES.items()))
    (BASE/'sources/manifest.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
    print(json.dumps(records,indent=2))
