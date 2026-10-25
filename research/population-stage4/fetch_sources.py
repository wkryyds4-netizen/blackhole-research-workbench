"""Fetch bounded source history; retain source as data, never import it."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import requests

BASE = Path(__file__).resolve().parent
OUT = BASE / 'sources'
SOURCES = {
 'pipe-history.json': 'https://git.ligo.org/api/v4/projects/RatesAndPopulations%2Fgwpopulation_pipe/repository/commits?path=gwpopulation_pipe/data_collection.py&until=2024-06-15T00:00:00Z&per_page=5',
 'pipe-pypi.json': 'https://pypi.org/pypi/gwpopulation_pipe/json',
 'bilby-cosmological.py.txt': 'https://raw.githubusercontent.com/bilby-dev/bilby/v2.2.1/bilby/gw/prior.py',
 'pipe-april2024.py.txt': 'https://git.ligo.org/api/v4/projects/RatesAndPopulations%2Fgwpopulation_pipe/repository/files/gwpopulation_pipe%2Fdata_collection.py/raw?ref=5dbcb49d602393f23d91e4d2b9a4cbae52d843a8',
 'pipe-may2024.py.txt': 'https://git.ligo.org/api/v4/projects/RatesAndPopulations%2Fgwpopulation_pipe/repository/files/gwpopulation_pipe%2Fdata_collection.py/raw?ref=20651a892864f8bb7e23adea33a8a7cf379fdde8',
 'pesummary-cosmology.py.txt': 'https://git.ligo.org/api/v4/projects/lscsoft%2Fpesummary/repository/files/pesummary%2Fgw%2Fcosmology.py/raw?ref=master',
}

def fetch(item):
    name, url = item
    path = OUT / name
    rec = dict(file=name, url=url, observed_utc=datetime.now(timezone.utc).isoformat())
    try:
        if not path.exists():
            r = requests.get(url, timeout=(15,35))
            r.raise_for_status()
            assert len(r.content) < 5_000_000
            path.write_bytes(r.content)
        data = path.read_bytes()
        rec.update(bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
    except (requests.RequestException, AssertionError) as e:
        rec['error'] = str(e)
    return rec

if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=3) as pool:
        records = list(pool.map(fetch,SOURCES.items()))
    (OUT/'manifest.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
    print(json.dumps(records,indent=2))
