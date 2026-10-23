"""Summarize pinned notebook history and release file hashes, without execution."""
from pathlib import Path
import hashlib,json
BASE=Path(__file__).resolve().parent;S=BASE/'provenance-sources'
def read(name):return json.loads((S/name).read_text())
old=read('zenodo-5636816.json');new=read('zenodo-7890398.json')
a={f['key']:f['checksum'] for f in old['files'] if f['key'].endswith('.hdf5')}
b={f['key']:f['checksum'] for f in new['files'] if f['key'].endswith('.hdf5')}
history=read('notebook-history.json')
notebooks={}
for file in ['notebook-original.ipynb','notebook-current.ipynb']:
    notebook=read(file)
    text='\n'.join(''.join(c['source']) for c in notebook['cells'])
    notebooks[file]=dict(sha256=hashlib.sha256((S/file).read_bytes()).hexdigest(),
       source_cells_with_download_or_preprocess=[dict(cell_index=i,text=''.join(c['source'])) for i,c in enumerate(notebook['cells']) if any(x in ''.join(c['source']) for x in ['gdown','pre-processed','gwpopulation-pipe'])],
       includes_local_google_drive_ids=all(x in text for x in ['16gStLIjt65gWBkw-gNOVUqNbZ89q8CLF','10pevUCM3V2-D-bROFEMAcTJsX_9RzeM6']),
       has_original_hdf_filename='real+semianalytic' in text,
       has_threshold_arguments=any(x in text for x in ['ifar_threshold','snr_threshold','vt-ifar-threshold','vt-snr-threshold']))
result=dict(schema='blackhole-population-provenance-history/1',date='2026-10-24',
   notebook_path_history=[dict(sha=x['sha'],date=x['commit']['author']['date'],message=x['commit']['message'].splitlines()[0]) for x in history],
   original_pinned_url='https://github.com/ColmTalbot/gwpopulation/blob/2daa71e8c56362f2f64e2166a0055478557f239a/examples/GWTC3.ipynb',
   notebook_checks=notebooks,
   conclusion='The earliest commit returned for this path (2024 documentation restructure) and current notebook both point to the same Drive objects; neither contains the original HDF filename or detection-threshold invocation. This is not an exhaustive claim about private or unrelated repositories.',
   release_comparison=dict(records=[old['doi'],new['doi']],all_hdf_checksums_identical=a==b,hdf_checksums=b,
      source_metadata_files=['provenance-sources/zenodo-5636816.json','provenance-sources/zenodo-7890398.json']),
   snapshot_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in S.iterdir() if p.is_file()})
(BASE/'provenance-history.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(dict(original_and_current_same_drive_ids=all(d['includes_local_google_drive_ids'] for d in notebooks.values()),hdf_hashes_unchanged=a==b),indent=2))
