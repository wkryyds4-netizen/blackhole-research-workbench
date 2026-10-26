"""Download one pinned official O3 mixture; never execute upstream content."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import hashlib
import json
import time
import requests

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
NAME = 'endo3_mixture-LIGO-T2100113-v12.hdf5'
SIZE = 410957336
MD5 = 'ebbae2546c847a8656022f4c067fb266'


def main():
    metadata_path = ROOT/'research/population-selection-audit/sources/zenodo-7890437.json'
    meta = json.loads(metadata_path.read_text(encoding='utf-8'))
    assert meta['id'] == 7890437
    entry = next(x for x in meta['files'] if x['key'] == NAME)
    assert entry['size'] == SIZE and entry['checksum'] == 'md5:' + MD5
    folder = BASE/'sources'
    folder.mkdir(exist_ok=True)
    dest = folder/NAME
    chunk_size = 4*1024*1024
    starts = list(range(0, SIZE, chunk_size))
    chunks = folder/(NAME + '.chunks')
    if not dest.exists():
        chunks.mkdir(exist_ok=True)
        def fetch(start):
            end = min(start+chunk_size, SIZE)-1
            part = chunks/f'{start}.part'
            if part.exists() and part.stat().st_size == end-start+1:
                return
            for attempt in range(3):
                try:
                    with requests.get(entry['links']['self'], headers={'Range':f'bytes={start}-{end}'},
                                      timeout=(15,40), stream=True) as response:
                        response.raise_for_status()
                        assert response.status_code == 206
                        assert response.headers.get('Content-Range') == f'bytes {start}-{end}/{SIZE}'
                        raw = b''.join(response.iter_content(65536))
                    assert len(raw) == end-start+1
                    part.write_bytes(raw)
                    return
                except (requests.RequestException, AssertionError):
                    if attempt == 2:
                        raise
                    time.sleep(1+attempt)
        with ThreadPoolExecutor(max_workers=4) as pool:
            jobs = [pool.submit(fetch, start) for start in starts]
            for n, job in enumerate(as_completed(jobs), 1):
                job.result()
                if n % 5 == 0 or n == len(jobs):
                    print(f'O3 mixture: {n}/{len(jobs)} chunks', flush=True)
        temporary = folder/(NAME + '.assembled')
        with temporary.open('wb') as out:
            for start in starts:
                out.write((chunks/f'{start}.part').read_bytes())
        with temporary.open('rb') as stream:
            assert hashlib.file_digest(stream, 'md5').hexdigest() == MD5
        assert temporary.stat().st_size == SIZE
        temporary.replace(dest)
    with dest.open('rb') as stream:
        assert hashlib.file_digest(stream, 'md5').hexdigest() == MD5
    with dest.open('rb') as stream:
        sha = hashlib.file_digest(stream, 'sha256').hexdigest()
    report = dict(file=NAME, url=entry['links']['self'], bytes=SIZE,
                  official_md5=MD5, sha256=sha,
                  fixed_record_metadata_sha256=hashlib.sha256(metadata_path.read_bytes()).hexdigest())
    (BASE/'download.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
