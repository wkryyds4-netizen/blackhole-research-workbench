"""Deterministic control using recovered public injection mixture weights."""
from pathlib import Path
import json

import numpy as np

from population_mc_bootstrap import ROOT, KEYS, prepare_weights, digest
from population_functionals import aligned_cdf, functional_coefficients, log_likelihood, normalize, summarize, transform


def evaluate(data, injection, mixture_weight, total, points):
    f = np.linspace(0, 1, points)
    lower = float(aligned_cdf(0, 1.2))
    r = np.linspace(lower, .5, points)
    rows = {}
    weighted = (injection*mixture_weight[:, None, None]).sum(axis=0)/total
    for ki, key in enumerate(KEYS):
        sigma = [.3, .6, 1.2][int(key[0])]
        rec = data['scenarios'][key]
        amended = {**rec, 'selection': weighted[ki].tolist()+rec['selection'][2:]}
        original_xi = np.asarray(rec['selection'][:2])
        fp = (.5-r)/(.5-float(aligned_cdf(0, sigma)))
        row = dict(selection_relative_changes=(weighted[ki]/original_xi-1).tolist())
        for scope, x, fg in [('native_r', f, f), ('common_r', r, fp)]:
            old_logs = log_likelihood(rec, fg)
            new_logs = log_likelihood(amended, fg)
            pold, pnew = normalize(x, old_logs), normalize(x, new_logs)
            if scope == 'native_r':
                yold, pyold = transform(x, pold, *functional_coefficients(sigma)['r'])
                ynew, pynew = transform(x, pnew, *functional_coefficients(sigma)['r'])
            else:
                yold, pyold, ynew, pynew = x, pold, x, pnew
            before, after = summarize(yold, pyold), summarize(ynew, pynew)
            row[scope] = dict(before_quantiles=before['quantiles'], after_quantiles=after['quantiles'],
                              quantile_changes=(np.array(after['quantiles'])-before['quantiles']).tolist(),
                              total_variation=float(.5*np.trapezoid(abs(pynew-pyold), yold)),
                              log_likelihood_correction_range=float(np.ptp(new_logs-old_logs)))
        rows[key] = row
    return rows


def main():
    out = ROOT/'research/population-stage2'
    data = json.loads((ROOT/'results/model.json').read_text(encoding='utf-8'))
    source = json.loads((out/'provenance-injection-match.json').read_text())
    with np.load(out/'provenance-injection-row-map.npz', allow_pickle=False) as rowmap:
        mix = rowmap['mixture_weight']
    _, injection, total, error = prepare_weights(data)
    assert len(injection) == len(mix) == source['row_match']['n_unique_original_rows']
    assert digest(ROOT/'data/injections.npz') == source['local_npz_sha256']
    assert mix.min() > 0
    rows = evaluate(data, injection, mix, total, 4001)
    fine = evaluate(data, injection, mix, total, 8001)
    grid_error = max(float(np.max(abs(np.array(rows[k][s]['quantile_changes'])-fine[k][s]['quantile_changes'])))
                     for k in KEYS for s in ['native_r', 'common_r'])
    bound = len(data['events'])*np.log(mix.max()/mix.min())
    assert all(v[s]['log_likelihood_correction_range'] <= bound+1e-11 for v in rows.values() for s in ['native_r', 'common_r'])
    assert grid_error < 1e-9
    # A parameter-independent efficiency rescaling must cancel upon normalization.
    rec = data['scenarios'][data['default']]
    ratio = total/source['raw_hdf_attributes']['total_generated']
    rescaled = {**rec, 'selection': [v*ratio for v in rec['selection'][:2]]+rec['selection'][2:]}
    f = np.linspace(0, 1, 4001)
    p0 = normalize(f, log_likelihood(rec, f))
    p1 = normalize(f, log_likelihood(rescaled, f))
    total_variation_scale = float(.5*np.trapezoid(abs(p0-p1), f))
    assert total_variation_scale < 1e-12
    result = dict(schema='recovered-mixture-weight-control/1', date='2026-10-24', points=4001, scenarios=rows,
                  correction='multiply source mixture_weight into injection importance weights; original inputs unchanged',
                  max_raw_moment_relative_error=error,
                  max_common_r_median_abs_change=max(abs(v['common_r']['quantile_changes'][1]) for v in rows.values()),
                  max_native_r_median_abs_change=max(abs(v['native_r']['quantile_changes'][1]) for v in rows.values()),
                  max_4001_8001_correction_difference=grid_error,
                  log_likelihood_correction_range_bound=float(bound),
                  noninteger_normalization_total_variation=total_variation_scale,
                  input_sha256={str(p.relative_to(ROOT)).replace('\\', '/'): digest(p) for p in
                                [ROOT/'results/model.json', out/'provenance-injection-match.json', out/'provenance-injection-row-map.npz']},
                  script_sha256=digest(Path(__file__)),
                  limitation='This isolates recovered near-unity mixture weights, not all selection error or true injection sampling covariance.')
    (out/'mixture-weight-control.json').write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps({k: v for k, v in result.items() if k not in ['scenarios', 'input_sha256']}, indent=2))


if __name__ == '__main__':
    main()
