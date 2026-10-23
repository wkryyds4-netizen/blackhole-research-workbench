"""Download public GWOSC catalogs and GWPopulation tutorial data only."""
from pathlib import Path
import hashlib
import json
import requests

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "data"
OUT.mkdir(exist_ok=True)
sources = {
    "GWTC-1-confident.json": "https://gwosc.org/eventapi/json/GWTC-1-confident/",
    "GWTC-2.1-confident.json": "https://gwosc.org/eventapi/json/GWTC-2.1-confident/",
    "GWTC-3-confident.json": "https://gwosc.org/eventapi/json/GWTC-3-confident/",
    "gwtc-3-samples.pkl": "https://drive.google.com/uc?export=download&id=10pevUCM3V2-D-bROFEMAcTJsX_9RzeM6",
    "gwtc-3-injections.pkl": "https://drive.google.com/uc?export=download&id=16gStLIjt65gWBkw-gNOVUqNbZ89q8CLF",
}
manifest = []
for name, url in sources.items():
    p = OUT / name
    try:
        if not p.exists():
            r = requests.get(url, timeout=90)
            r.raise_for_status()
            if name.endswith('.pkl') and r.content[:1] != b'\x80':
                raise ValueError('Download is not a binary pickle; refusing to save')
            if name.endswith('.json'):
                r.json()
            p.write_bytes(r.content)
        b = p.read_bytes()
        manifest.append(dict(file=name, url=url, bytes=len(b), sha256=hashlib.sha256(b).hexdigest(), retrieved='2026-10-22'))
        print(name, len(b), flush=True)
    except Exception as exc:
        print(name, type(exc).__name__, str(exc)[:160], flush=True)
(OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
