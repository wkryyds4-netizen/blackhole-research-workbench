"""HDF-capable environment: independently check the serialized interface rows."""
from pathlib import Path
import json
import hashlib
import numpy as np
import h5py

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]


def main():
    path=BASE/'gw190805-weighted-event.npz'
    with np.load(path,allow_pickle=False) as a:
        meta=json.loads(str(a['metadata']))
        source=ROOT/'research/population-stage5/source-data'/meta['source']['file']
        assert hashlib.sha256(source.read_bytes()).hexdigest()==meta['source']['sha256']
        with h5py.File(source,'r') as h:
            samples=h[meta['analysis_label']]['posterior_samples'][:]
            rows=a['row_ids']
            assert np.array_equal(rows,np.arange(len(samples)))
            checked=[]
            for name in meta['measure']['coordinates']:
                original='mass_1_source' if name=='mass_1' else name
                assert np.array_equal(a['sample_'+name],samples[original][rows])
                checked.append(name)
    result=dict(schema='blackhole-weighted-source-row-validation/1',status='pass',
                source_sha256=meta['source']['sha256'],n_rows=len(rows),
                exact_coordinates=checked,all_original_rows_in_original_order=True,
                cache_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    (BASE/'source-row-validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
