"""Exact seven-coordinate PE provenance pilots; never execute prior strings.

Report stored tilt transformations as event-posterior diagnostics, not population
fractions or a substitute for consistent event AND injection transformations.
"""
from pathlib import Path
import hashlib
import json
import h5py
import numpy as np

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
PARAMS = ['mass_1', 'mass_ratio', 'redshift', 'a_1', 'a_2', 'cos_tilt_1', 'cos_tilt_2']
SOURCE = dict(zip(PARAMS, ['mass_1_source'] + PARAMS[1:]))

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def metadata(group, path):
    if path not in group:
        return None
    if isinstance(group[path], h5py.Group):
        return {key: metadata(group[path], key) for key in group[path]}
    value = group[path][()].tolist()
    if not isinstance(value, list):
        value = [value]
    return [v.decode() if isinstance(v, bytes) else v for v in value]

def join(local, source, keys):
    left = np.rec.fromarrays([local[k] for k in keys])
    right = np.rec.fromarrays([source[SOURCE[k]] for k in keys])
    unique, first, counts = np.unique(right, return_index=True, return_counts=True)
    pos = np.searchsorted(unique, left)
    bounded = np.minimum(pos, len(unique) - 1)
    found = (pos < len(unique)) & (unique[bounded] == left)
    ambiguous = found & (counts[bounded] != 1)
    return first[bounded], found, ambiguous

def main():
    downloads = json.loads((BASE / 'pilot-downloads.json').read_text())
    audit = json.loads((ROOT / 'data/sample_audit.json').read_text())
    events, rows, plot_data = [], {}, {}
    with np.load(ROOT / 'data/posteriors.npz', allow_pickle=False) as archive:
        for entry in downloads:
            event = next(x for x in audit if x['name'] == entry['event'])
            index = event['index']
            local = {k: archive[f'{index}_{k}'] for k in PARAMS}
            path = BASE / 'source-data' / entry['file']
            assert digest(path) == entry['sha256']
            labels = []
            with h5py.File(path, 'r') as h:
                for name in h:
                    group = h[name]
                    if not isinstance(group, h5py.Group) or 'posterior_samples' not in group:
                        continue
                    source = group['posterior_samples'][:]
                    indices, found, ambiguous = join(local, source, ['mass_1', 'redshift'])
                    checked = found & ~ambiguous
                    columns = {}
                    for k in PARAMS:
                        equal = local[k] == source[SOURCE[k]][indices]
                        columns[k] = int(np.sum(checked & equal))
                        checked &= equal
                    full = bool(np.all(checked))
                    record = dict(label=name, source_rows=len(source), local_rows=len(indices),
                        join_matches=int(found.sum()), ambiguous_join_matches=int(ambiguous.sum()),
                        exact_seven_coordinate_matches=int(checked.sum()), full_match=full,
                        progressive_equal_counts=columns,
                        approximant=metadata(group, 'approximant'),
                        f_ref_hz=metadata(group, 'meta_data/meta_data/f_ref'),
                        config_reference_frequency=metadata(group, 'config_file/config/reference-frequency'),
                        cosmology=metadata(group, 'meta_data/meta_data/cosmology'))
                    if full:
                        independent, found2, ambiguous2 = join(local, source, ['a_1', 'redshift'])
                        assert found2.all() and not ambiguous2.any()
                        assert np.array_equal(independent, indices)
                        assert len(np.unique(indices)) == len(indices)
                        record['independent_a1_redshift_join_agrees'] = True
                        record['analytic_prior_strings'] = {k: metadata(group, 'priors/analytic/' + k)
                            for k in ['chirp_mass','mass_ratio','mass_1','mass_2',
                                      'luminosity_distance','a_1','a_2','tilt_1','tilt_2']}
                        record['backward_spin_evolution_metadata'] = metadata(group, 'meta_data/meta_data/backward_spin_evolution')
                        rows[f'{index}_original_row_indices'] = indices
                        diagnostics = {}
                        for component in [1, 2]:
                            k = f'cos_tilt_{component}'
                            inf_key = k + '_infinity_only_prec_avg'
                            if inf_key not in source.dtype.names:
                                continue
                            now, inf = local[k], source[inf_key][indices]
                            finite = np.isfinite(now) & np.isfinite(inf)
                            assert finite.all(), 'Non-finite tilt conversion requires an explicit missing-data policy'
                            plot_data[f'{index}_{k}'] = now
                            plot_data[f'{index}_{k}_infinity'] = inf
                            diagnostics[k] = dict(n=len(now),
                                reference_negative_fraction=float(np.mean(now < 0)),
                                infinity_negative_fraction=float(np.mean(inf < 0)),
                                sign_crossing_count=int(np.sum((now < 0) != (inf < 0))),
                                sign_crossing_fraction=float(np.mean((now < 0) != (inf < 0))),
                                delta_cosine_quantiles=np.quantile(inf-now, [.05,.5,.95]).tolist())
                        record['stored_infinity_tilt_diagnostic'] = diagnostics
                    labels.append(record)
            assert sum(x['full_match'] for x in labels) == 1
            events.append(dict(event=entry['event'], original_index=index, release_record=entry['record'],
                               source_file=entry['file'], source_sha256=entry['sha256'], labels=labels))
    result = dict(schema='blackhole-population-pe-pilot/1',
        local_npz_sha256=digest(ROOT / 'data/posteriors.npz'), events=events,
        interpretation='Exact finite-file sample membership, not recovery of historical commands or software version.',
        tilt_diagnostic_scope='Paired event posterior draws only; not population r, not 10-to-20 Hz evolution; no injection transformation applied.',
        prior_caution='Analytic strings can describe original sampling priors; cosmo releases reweight distance. Do not treat strings as final density.')
    np.savez_compressed(BASE / 'pilot-row-map.npz', **rows)
    np.savez_compressed(BASE / 'pilot-tilt-diagnostics.npz', **plot_data)
    (BASE / 'pe-pilot-match.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps([dict(event=e['event'], matches=[dict(label=l['label'], n=l['exact_seven_coordinate_matches'],
          f_ref=l['f_ref_hz'], tilt=l.get('stored_infinity_tilt_diagnostic')) for l in e['labels']]) for e in events], indent=2))

if __name__ == '__main__':
    main()
