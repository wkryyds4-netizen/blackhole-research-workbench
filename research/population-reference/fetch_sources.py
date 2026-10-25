"""Bounded public documentary lookup and checksum-pinned precession wheel."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import requests

BASE=Path(__file__).resolve().parent


def main():
    folder=BASE/'sources';folder.mkdir(parents=True,exist_ok=True)
    manifest=folder/'manifest.json'
    previous=json.loads(manifest.read_text()) if manifest.exists() else []
    targets=[
        ('mould-gerosa-2110.05507v2.html','https://arxiv.org/html/2110.05507v2',True),
        ('gwtc3-populations-2111.03634v3.html','https://arxiv.org/html/2111.03634v3',True),
        ('precession-pypi-2.1.2.json','https://pypi.org/pypi/precession/2.1.2/json',True),
        ('ligo-generation-project','https://git.ligo.org/api/v4/projects/RatesAndPopulations%2Flvc-rates-and-pop',False),
        ('ligo-generation-commit','https://git.ligo.org/RatesAndPopulations/lvc-rates-and-pop/-/tree/28446706ad69ff815a2ddc0b1dc3cf813215899c',False),
        ('ligo-o3-document','https://dcc.ligo.org/LIGO-T2100113/public',False)]
    # Preserve additional documentary sources and retrieve them on a clean checkout.
    target_urls={t[1] for t in targets}
    for entry in previous:
        if entry.get('file') and not entry['file'].endswith('.whl') and entry['url'] not in target_urls:
            targets.append((entry['file'],entry['url'],True));target_urls.add(entry['url'])
    def fetch(item):
        name,url,save=item;path=folder/name
        if save and path.exists():
            content=path.read_bytes();sha=hashlib.sha256(content).hexdigest()
            old=next((e for e in previous if e.get('file')==name),{})
            if old.get('sha256'):assert sha==old['sha256'],f'Cached source changed: {name}'
            return dict(old, file=name,url=url,bytes=len(content),sha256=sha,status='cached')
        try:
            r=requests.get(url,timeout=(10,25),allow_redirects=True)
            row=dict(file=name if save else None,url=url,http_status=r.status_code,final_url=r.url,bytes=len(r.content),content_type=r.headers.get('Content-Type'))
            assert len(r.content)<20_000_000
            if r.ok and save:
                old=next((e for e in previous if e.get('file')==name),{})
                if old.get('sha256'):assert hashlib.sha256(r.content).hexdigest()==old['sha256'],f'Remote source changed: {name}'
                path.write_bytes(r.content);row['sha256']=hashlib.sha256(r.content).hexdigest()
            if not save:
                old=next((e for e in previous if e.get('url')==url),{})
                if old.get('interpretation'):row['interpretation']=old['interpretation']
            return row
        except requests.RequestException as exc:return dict(url=url,status='unavailable',error=type(exc).__name__)
    with ThreadPoolExecutor(max_workers=4) as pool:records=list(pool.map(fetch,targets))
    meta=json.loads((folder/'precession-pypi-2.1.2.json').read_text())
    assert meta['info']['version']=='2.1.2'
    entry=next(x for x in meta['urls'] if x['packagetype']=='bdist_wheel')
    assert entry['size']<5_000_000
    dest=folder/entry['filename'];assert dest.parent==folder
    if not dest.exists():
        r=requests.get(entry['url'],timeout=(10,30));r.raise_for_status()
        assert len(r.content)==entry['size'] and hashlib.sha256(r.content).hexdigest()==entry['digests']['sha256']
        dest.write_bytes(r.content)
    content=dest.read_bytes();assert len(content)==entry['size'] and hashlib.sha256(content).hexdigest()==entry['digests']['sha256']
    records.append(dict(file=dest.name,url=entry['url'],bytes=len(content),sha256=entry['digests']['sha256'],official_pypi_digest_verified=True))
    (folder/'manifest.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
    print(json.dumps(records,indent=2),flush=True)


if __name__=='__main__':main()
