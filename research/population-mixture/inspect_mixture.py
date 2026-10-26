"""Audit exact pure-BBH membership in the fixed official O3 mixture."""
from pathlib import Path
import hashlib
import json
import h5py
import numpy as np

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
OLD = ROOT/'research/population-selection-audit'


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def plain(value):
    if isinstance(value, bytes): return value.decode()
    if isinstance(value, np.ndarray): return [plain(x) for x in value]
    if isinstance(value, np.generic): return plain(value.item())
    return value


def main():
    source = json.loads((BASE/'download.json').read_text(encoding='utf-8'))
    path = BASE/'sources'/source['file']
    assert digest(path) == source['sha256']
    with h5py.File(path, 'r') as handle:
        group = handle['injections']
        attrs = {k:plain(v) for k,v in handle.attrs.items()}
        group_attrs = {k:plain(v) for k,v in group.attrs.items()}
        descriptors = {k:dict(shape=list(d.shape), dtype=str(d.dtype)) for k,d in group.items()}
        fields = ['mass1_source','mass2_source','redshift']
        target = np.rec.fromarrays([group[k][:].astype(float) for k in fields])
        assert len(np.unique(target)) == len(target)
        order = np.argsort(target)
        with np.load(OLD/'o3-bank.npz', allow_pickle=False) as old:
            wanted = np.rec.fromarrays([old['full_'+k].astype(float) for k in fields])
            positions = np.searchsorted(target[order], wanted)
            assert np.all(positions < len(target))
            indices = order[positions]
            assert np.array_equal(target[indices], wanted)
            assert len(np.unique(indices)) == len(indices) == 284672
            checked = {}
            for k in group:
                key = 'full_'+k
                if key not in old.files or k in ['sampling_pdf','mixture_weight'] or k.endswith('_sampling_pdf'):
                    continue
                actual = group[k][:][indices]; expected = old[key]
                equal = np.array_equal(actual, expected, equal_nan=True)
                checked[k] = dict(exact=bool(equal))
                if np.issubdtype(actual.dtype, np.number):
                    finite = np.isfinite(actual) & np.isfinite(expected)
                    checked[k]['maximum_absolute_difference'] = float(np.max(np.abs(actual[finite]-expected[finite]))) if finite.any() else 0.
            assert all(x['exact'] for x in checked.values()), checked
        search = [k for k in group if k.startswith('ifar_')]
        union = np.any([group[k][:] > 1 for k in search], axis=0)
        matched = np.any([group[k][:] > 1 for k in search if k != 'ifar_cwb'], axis=0)
        z = group['redshift'][:]; m = group['mass1_source'][:]; n = group['mass2_source'][:]
        a = [np.sqrt(sum(group[f'spin{i}{ax}'][:].astype(float)**2 for ax in 'xyz')) for i in [1,2]]
        target_mass = (m >= 3)&(m <= 100)&(n >= 3)&(n <= m)&(z <= 2.3)
        bbh_domain = (m >= 2)&(m <= 100)&(n >= 2)&(z <= 1.9)&(a[0] <= .998)&(a[1] <= .998)
        regions = {}
        for label, mask in dict(all=np.ones(len(z),bool), target_mass=target_mass,
                                target_new_union=target_mass&~bbh_domain,
                                target_redshift_tail=target_mass&(z > 1.9)).items():
            regions[label] = dict(rows=int(mask.sum()), matched_filter=int((mask&matched).sum()),
                                 with_cwb=int((mask&union).sum()))
        report = dict(schema='blackhole-o3-mixture-source-audit/1', source=source,
                      attributes=attrs, group_attributes=group_attrs,
                      root_group_attributes_equal=(attrs == group_attrs), fields=descriptors,
                      pure_bbh_exact_matches=len(indices), exact_shared_fields=checked,
                      regions=regions,
                      mixture_weight_range=[float(group['mixture_weight'][:].min()),float(group['mixture_weight'][:].max())],
                      redshift_range=[float(z.min()),float(z.max())],
                      candidate_reference_fields=[k for k in list(group)+list(attrs) if any(v in k.lower() for v in ['freq','waveform','approximant','reference'])],
                      source_identity_scope='Exact membership/metadata only; no recovered waveform, reference or hidden-trial detection law.')
        assert attrs['n_accepted'] == len(z)
        assert attrs['n_accepted'] + attrs['n_rejected'] == attrs['total_generated']
    np.savez_compressed(BASE/'mixture-membership.npz', pure_bbh_to_mixture=indices)
    report['mapping_sha256'] = digest(BASE/'mixture-membership.npz')
    report['producer_sha256'] = digest(Path(__file__))
    (BASE/'provenance.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['attributes','pure_bbh_exact_matches','regions','candidate_reference_fields']}, indent=2))


if __name__ == '__main__': main()
