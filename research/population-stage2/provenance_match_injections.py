"""Compare candidate public HDF rows and explicit masks with safe local NPZ.

Matching is by exact (primary mass, redshift) with uniqueness checks, then independent
checks of all seven physical columns and the stored draw probability.
"""
from pathlib import Path
import hashlib,json
import h5py,numpy as np
BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
path=BASE/'source-data/o1+o2+o3_bbhpop_real+semianalytic-LIGO-T2100377-v2.hdf5'
assert hashlib.md5(path.read_bytes()).hexdigest()=='c21337a25b7318fdfd471cce549a547e','Official HDF MD5 mismatch'
local=dict(np.load(ROOT/'data/injections.npz',allow_pickle=False))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def stats(a):
    a=np.asarray(a);return dict(min=float(a.min()),median=float(np.median(a)),max=float(a.max()))
def compare(a,b):
    assert len(a)==len(b)
    scale=np.maximum(np.maximum(abs(a),abs(b)),1e-300)
    return dict(exact_equal=bool(np.array_equal(a,b)),max_abs=float(np.max(abs(a-b))),
                max_relative=float(np.max(abs(a-b)/scale)),allclose_rtol_1e_13_atol_1e_15=bool(np.allclose(a,b,rtol=1e-13,atol=1e-15)))
with h5py.File(path,'r') as h:
    g=h['injections']; raw={k:g[k][:] for k in g}
    attrs={k:float(v) for k,v in g.attrs.items()}
m=raw['mass1_source'];join=np.rec.fromarrays([m,raw['redshift']]);order=np.argsort(join)
local_join=np.rec.fromarrays([local['mass_1'],local['redshift']])
assert len(np.unique(join))==len(m),'Two-coordinate identification would be ambiguous'
positions=np.searchsorted(join[order],local_join)
assert np.all(positions<len(m))
indices=order[positions]
assert np.array_equal(join[indices],local_join),'Candidate does not exactly contain all local mass/redshift pairs'
assert len(np.unique(indices))==len(indices)
params=dict(mass_1=m,mass_ratio=raw['mass2_source']/m,redshift=raw['redshift'])
for i in [1,2]:
    a=np.sqrt(sum(raw[f'spin{i}{axis}']**2 for axis in 'xyz'))
    params[f'a_{i}']=a;params[f'cos_tilt_{i}']=raw[f'spin{i}z']/a
checks={k:compare(local[k],v[indices]) for k,v in params.items()}
assert all(c['allclose_rtol_1e_13_atol_1e_15'] for c in checks.values())
assert all(c['exact_equal'] for c in checks.values()),'Physical coordinates not exactly reproduced'
names=raw['name'];old=np.isin(names,[b'o1',b'o2'])
ifar_keys=sorted(k for k in raw if k.startswith('ifar_'))
far_keys=sorted(k for k in raw if k.startswith('far_'))
mask=np.any(np.array([raw[k]>1 for k in ifar_keys]),axis=0)|(old&(raw['optimal_snr_net']>10))
actual=np.zeros(len(m),dtype=bool);actual[indices]=True
mask_far=np.any(np.array([raw[k]<1 for k in far_keys]),axis=0)|(old&(raw['optimal_snr_net']>10))
prior=raw['sampling_pdf']*m*(2*np.pi*params['a_1']**2)*(2*np.pi*params['a_2']**2)
weights=raw['mixture_weight']
prior_checks={
 'jacobian_only':compare(local['prior'],prior[indices]),
 'jacobian_divided_by_mixture_weight':compare(local['prior'],(prior/weights)[indices]),
}
assert np.array_equal(mask,actual),'Explicit selection rule does not reproduce source rows'
assert prior_checks['jacobian_only']['exact_equal'],'Candidate Jacobian does not reproduce stored prior'
run_info={}
for run in np.unique(names):
    r=names==run; sel=r&actual
    run_info[run.decode()]=dict(n_accepted=int(r.sum()),n_retained=int(sel.sum()),
        mixture_weight_all=stats(weights[r]),mixture_weight_retained=stats(weights[sel]),
        n_distinct_mixture_weights=int(len(np.unique(weights[r]))),
        snr_retained=stats(raw['optimal_snr_net'][sel]))
np.savez_compressed(BASE/'provenance-injection-row-map.npz',original_row_indices=indices,
                    run=np.array([s.decode() for s in names[indices]],dtype='U2'),mixture_weight=weights[indices])
record=dict(schema='blackhole-population-injection-match/1',date='2026-10-24',
    candidate_source_file=str(path.relative_to(BASE)),source_sha256=digest(path),
    local_npz_sha256=digest(ROOT/'data/injections.npz'),official_md5=hashlib.md5(path.read_bytes()).hexdigest(),
    raw_hdf_attributes=attrs,
    metadata_comparison=dict(local_total_generated=int(local['total_generated']),
        raw_total_generated=attrs['total_generated'],
        int_truncation_exact_match=int(attrs['total_generated'])==int(local['total_generated']),
        relative_truncation=(attrs['total_generated']-int(local['total_generated']))/attrs['total_generated'],
        local_analysis_time_years=float(local['analysis_time']),
        converted_analysis_time_years=attrs['analysis_time_s']/365.25/24/60/60,
        time_exact_match=attrs['analysis_time_s']/365.25/24/60/60==float(local['analysis_time']),
        warning='The original total_generated is fractional. It is effective release metadata, not proof of a literal IID trial count.'),
    row_match=dict(method='Exact primary-mass/redshift pair join, all candidate pairs unique, followed by seven physical-coordinate checks',
        source_primary_mass_unique_count=len(np.unique(m)),source_primary_mass_redshift_pair_unique_count=len(np.unique(join)),
        n_exact_primary_mass_matches=len(indices),n_unique_original_rows=len(np.unique(indices)),
        original_order_preserved=bool(np.all(np.diff(indices)>0)),physical_columns=checks),
    selection_match=dict(explicit_rule='any of five IFAR fields > 1 yr OR (run in O1/O2 AND optimal_snr_net > 10)',
        ifar_keys=ifar_keys,far_keys=far_keys,
        local_retained_count=len(indices),rule_retained_count=int(mask.sum()),
        rule_selects_exact_same_rows=bool(np.array_equal(mask,actual)),
        rule_row_symmetric_difference=int(np.sum(mask!=actual)),
        far_less_than_one_rule_same_rows=bool(np.array_equal(mask_far,actual)),
        qualification='An exact empirical reconstruction on this finite HDF; historical preprocessing command, package version and uniqueness of threshold remain unproven.'),
    prior_match=prior_checks,
    prior_formula='sampling_pdf * mass1_source * (2*pi*a1**2) * (2*pi*a2**2)',
    conditional_likelihood_bound=dict(n_events=69,
       mixture_weight_min=float(weights[indices].min()),mixture_weight_max=float(weights[indices].max()),
       max_log_likelihood_shape_change=69*float(np.log(weights[indices].max()/weights[indices].min())),
       derivation='For nonnegative importance weights, xi_weighted/xi_unweighted lies between min(w) and max(w). Thus the range of the conditional log-likelihood correction is at most N*log(max(w)/min(w)).',
       total_generated_truncation_effect='A parameter-independent multiplicative selection normalization; cancels exactly in the normalized conditional-on-N posterior. Not necessarily irrelevant to rate inference.',
       limitation='Bounds only the recovered mixture-weight omission, not event prior errors, selection-policy differences or Monte Carlo sampling uncertainty.'),
    run_summary=run_info,
    iid_status='Not established: source combines disjoint runs with run-dependent mixture weights; original run design counts/configuration absent from root attributes.',
    output_row_map='provenance-injection-row-map.npz',
    changes_to_original_inputs='none')
(BASE/'provenance-injection-match.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
print(json.dumps(record,indent=2))
