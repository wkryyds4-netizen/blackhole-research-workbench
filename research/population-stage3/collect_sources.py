"""Freeze bounded official metadata for the event provenance audit."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import requests

BASE = Path(__file__).resolve().parent
OUT = BASE / 'sources'
SOURCES = {
    'population-2111.03634v3.html': 'https://arxiv.org/html/2111.03634v3',
    'zenodo-6513631.json': 'https://zenodo.org/api/records/6513631',
    'zenodo-5546663.json': 'https://zenodo.org/api/records/5546663',
}

def fetch(item):
    name, url = item
    path = OUT / name
    if not path.exists():
        response = requests.get(url, timeout=(15, 45))
        response.raise_for_status()
        path.write_bytes(response.content)
    data = path.read_bytes()
    return dict(file=name, url=url, bytes=len(data),
                sha256=hashlib.sha256(data).hexdigest(),
                observed_utc=datetime.now(timezone.utc).isoformat())

if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=3) as pool:
        results = list(pool.map(fetch, SOURCES.items()))
    (OUT / 'manifest.json').write_text(json.dumps(results, indent=2), encoding='utf-8')
    print(json.dumps(results, indent=2))
