"""Independent full-row set intersection and negative controls for PE matching."""
from pathlib import Path
import hashlib
import json
import h5py
import numpy as np

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
KEYS = ['mass_1', 'mass_ratio', 'redshift', 'a_1', 'a_2', 'cos_tilt_1', 'cos_tilt_2']
EXPECTED = {
    'results/model.json': '6dcf3c9f4fd51b2007977d131c61f69c7d861afb986551ea82b70a7667a51083',
    'data/posteriors.npz': '3aec709780ac373e485682323c20cbe6ed0c14028c66b030bb73d848739489c4',
    'data/injections.npz': '00ce46283aaf5c1e3df1a8b479a2f3c5bf904e4aef6aa0cf641f93df54bf7a13',
}

def main():
    for name, sha in EXPECTED.items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == sha
    for rec in json.loads((BASE / 'sources/manifest.json').read_text()):
        assert hashlib.sha256((BASE / 'sources' / rec['file']).read_bytes()).hexdigest() == rec['sha256']
    report = json.loads((BASE / 'pe-pilot-match.json').read_text())
    download = json.loads((BASE / 'pilot-downloads.json').read_text())
    results = []
    with np.load(ROOT / 'data/posteriors.npz', allow_pickle=False) as a:
        for e in report['events']:
            index = e['original_index']
            local = list(zip(*(a[f'{index}_{k}'].tolist() for k in KEYS)))
            assert len(set(local)) == len(local)
            file = BASE / 'source-data' / e['source_file']
            dl = next(x for x in download if x['event'] == e['event'])
            assert 'md5:' + hashlib.md5(file.read_bytes()).hexdigest() == dl['checksum']
            with h5py.File(file, 'r') as h:
                for label in e['labels']:
                    d = h[label['label'] + '/posterior_samples'][:]
                    source = set(zip(*(d['mass_1_source' if k == 'mass_1' else k].tolist() for k in KEYS)))
                    matches = len(set(local) & source)
                    assert matches == label['exact_seven_coordinate_matches']
                    if label['full_match']:
                        other = 48 if index == 38 else 38
                        wrong = set(zip(*(a[f'{other}_{k}'].tolist() for k in KEYS)))
                        assert len(source & wrong) == 0
                        # Keep the identification coordinates, deliberately break one tilt.
                        altered = [tuple(list(row[:-1]) + [2.0]) for row in local]
                        assert len(source & set(altered)) == 0
                    results.append(dict(event=e['event'], label=label['label'],
                                        full_row_intersection_count=matches))
    result = dict(schema='blackhole-population-pe-pilot-validation/1', status='pass',
                  frozen_hashes=EXPECTED, all_source_manifest_hashes_verified=True,
                  official_pilot_md5_verified=True, independent_set_intersections=results,
                  negative_controls=['wrong event: zero full-row matches',
                                     'same mass/redshift but corrupted tilt: zero full-row matches'])
    (BASE / 'validation.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    main()
